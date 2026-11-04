#!/usr/bin/env python3
"""Reconcile raw records, censor coding, method retention and real-data overlap."""
import csv
import hashlib
import itertools
import json
import math
from collections import Counter, defaultdict
from pathlib import Path
import re
import sys

out = Path(sys.argv[1]).resolve()
config = json.loads((out / 'config.json').read_text())
root = Path.cwd()

def rows(path):
    with path.open(newline='', encoding='utf-8-sig') as stream:
        return list(csv.DictReader(stream))

def digest(path):
    value = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1048576), b''): value.update(block)
    return value.hexdigest()

def number(value):
    try:
        parsed = float(value)
        return None if math.isnan(parsed) else parsed
    except (ValueError, TypeError): return None

def finite(value): return value is not None and math.isfinite(value)

def numeric_equal(actual, expected, tolerance=1e-10):
    a = number(actual)
    if expected is None: return a is None
    if a is None: return False
    return a == expected or math.isclose(a, expected, rel_tol=tolerance, abs_tol=tolerance)

def normalize(value):
    value = value.strip().upper()
    return value or None

def collection(value): return value.strip() or None

def key(row): return (row['species'], row['state'], None if row['method_collected'] == 'NA' else row['method_collected'])

def result_key(row): return (row['element'], row['state'], None if row['method_collected'] == 'NA' else row['method_collected'])

def expected_observation(row):
    original = row['qvalue']; effective = original; provenance = 'qvalue'
    first = row['summary'].split(',', 1)[0].strip()
    original_number = number(first); reported = number(original)
    suffix = r'\.1{4,5}$'
    if re.search(suffix, first) and finite(original_number) and finite(reported) and abs(round(original_number, 2)-reported) < 1e-8:
        effective = first; provenance = 'summary_upper_limit'
    text = effective.strip(); value = number(text)
    left = finite(value) and value < 0
    right = bool(re.search(suffix, text))
    lower_limit = abs(value) if left else None
    upper_limit = number(re.sub(suffix, '', text)) if right else None
    if left and right: status = 'conflicting_qualifiers'
    elif not text: status = 'missing_value'
    elif not finite(value): status = 'invalid_value'
    elif left and (not finite(lower_limit) or lower_limit <= 0): status = 'invalid_limit'
    elif right and (not finite(upper_limit) or upper_limit <= 0): status = 'invalid_limit'
    elif not left and not right and value <= 0: status = 'nonpositive_detection'
    else: status = 'usable'
    if status != 'usable': left = right = False; lower_limit = upper_limit = None
    factor = {'ppm':1000, 'mg/kg':1000, 'ppb':1, 'ug/kg':1, 'µg/kg':1, 'μg/kg':1,
              'pct':1e7, '%':1e7, 'percent':1e7}.get(row['unit'].strip().lower())
    detected = value if status == 'usable' and not left and not right else None
    if status == 'usable' and factor is None: status = 'unsupported_unit'
    converted = lambda value: value*factor if value is not None and factor is not None else None
    detected, lower_limit, upper_limit = map(converted, (detected, lower_limit, upper_limit))
    eligible = ((left and finite(lower_limit) and lower_limit > 0) or
                (right and finite(upper_limit) and upper_limit > 0) or
                (not left and not right and finite(detected) and detected > 0))
    code = normalize(row['amethod'])
    info = ('method_not_supplied' if code is None else 'method_unspecified' if code == 'UNKNOWN'
            else 'review_semiquantitative' if code in {'ES_SQ','ES_H2O_P'} else 'method_not_reviewed')
    return dict(status=status, left=left, right=right, detected=detected, left_limit=lower_limit,
                right_limit=upper_limit, eligible=eligible, method=code, method_information=info,
                semi=code in {'ES_SQ','ES_H2O_P'}, provenance=provenance, factor=factor)

for manifest in ['input_manifest.csv','source_manifest.csv']:
    for row in rows(out / manifest):
        assert digest(Path(row['file'])) == row['sha256'], ('changed input/source', row['file'])
raw_selected = {}; main_selected = {}; extraction_checks = []
for medium in config['media']:
    source = root / '../data' / ('ngdb' + medium)
    selected_metadata = []
    with (source / 'main.csv').open(newline='', encoding='utf-8-sig') as stream:
        for row in csv.DictReader(stream):
            if row['country'].strip().upper() in {'US','USA','UNITED STATES','U.S.'} and row['state'].strip().upper() == config['state']:
                selected_metadata.append(row)
    extracted_meta = rows(out / 'inputs' / medium / 'main_selected.csv')
    assert selected_metadata == extracted_meta, ('metadata selection mismatch', medium)
    meta = {row['lab_id'].strip(): row for row in selected_metadata}; main_selected[medium] = meta
    extracted = rows(out / 'inputs' / medium / 'chemistry_selected.csv')
    cursor = 0; source_rows = 0
    with (source / 'bestvalue.csv').open(newline='', encoding='utf-8-sig') as stream:
        for ordinal, raw in enumerate(csv.DictReader(stream), 1):
            source_rows += 1
            if raw['lab_id'].strip() not in meta: continue
            assert cursor < len(extracted)
            expected = dict(raw, source_record_number=str(ordinal))
            assert expected == extracted[cursor], ('raw field loss/change', medium, ordinal)
            cursor += 1
    assert cursor == len(extracted)
    raw_selected[medium] = {row['source_record_number']: row for row in extracted}
    extraction_checks.append(dict(medium=medium,n_source_records=source_rows,
        n_selected_metadata=len(meta),n_selected_records=cursor,all_raw_fields_reconciled=True))
    print(medium, 'raw records reconciled', cursor, flush=True)

checks = []; overlap_checks = []; status_totals = Counter()
for stage, medium in [('pilot','soil'),('expansion','soil'),('expansion','sed')]:
    folder = out / 'results' / stage / medium
    assert (folder / 'COMPLETED.txt').is_file(), ('incomplete phase', stage, medium)
    meta = main_selected[medium]
    selected = {ordinal: row for ordinal, row in raw_selected[medium].items()
        if stage != 'pilot' or (row['species'] in config['pilot_elements'] and
            meta[row['lab_id'].strip()]['method_collected'].strip() == config['pilot_collection'])}
    observed = rows(folder / 'observations.csv')
    assert len(observed) == len(selected)
    assert len({row['source_record_number'] for row in observed}) == len(observed)
    groups = defaultdict(list)
    for row in observed:
        raw = selected[row['source_record_number']]
        for field, value in raw.items():
            target = row['lab_id_raw'] if field == 'lab_id' else row[field]
            assert target == value, ('raw field changed', stage, medium, field, row['source_record_number'])
        assert row['lab_id'] == raw['lab_id'].strip()
        source_meta = meta[row['lab_id']]
        assert row['state'] == config['state']
        assert (None if row['method_collected'] == 'NA' else row['method_collected']) == collection(source_meta['method_collected'])
        predicted = expected_observation(raw)
        assert row['validation_status'] == predicted['status'], ('validation', medium, raw['rec_id'], row['validation_status'], predicted)
        assert (row['left_cen'] == 'TRUE') == predicted['left']
        assert (row['right_cen'] == 'TRUE') == predicted['right']
        for field, value in [('detected_val', predicted['detected']),('left_limit',predicted['left_limit']),
                             ('right_limit',predicted['right_limit']),('unit_factor',predicted['factor'])]:
            assert numeric_equal(row[field], value), ('parsed numeric field', medium, raw['rec_id'], field)
        assert row['censoring_source'] == predicted['provenance']
        assert row['method_review_status'] == predicted['method_information']
        assert row['review_semiquantitative'] == ('TRUE' if predicted['semi'] else 'NA')
        x, y = number(row['x_m']), number(row['y_m'])
        predicted['spatial'] = predicted['eligible'] and finite(x) and finite(y)
        predicted['x'], predicted['y'] = x, y
        groups[key(row)].append(predicted)
        status_totals[(stage, medium, predicted['status'])] += 1
    def result_rows(name):
        path = folder / (name + '.csv')
        if not path.exists():
            assert not observed, ('missing nonempty result table', stage, medium, name)
            return []
        return rows(path)
    estimates = result_rows('estimates'); support = result_rows('method_overlap_support')
    assert len(estimates) == len(groups) == len(support)
    for table in [estimates, support]:
        for row in table:
            values = groups[result_key(row)]; eligible = [value for value in values if value['eligible']]
            expected = {'n_input':len(values), 'n_total':len(eligible),
                'n_left':sum(value['left'] for value in eligible), 'n_right':sum(value['right'] for value in eligible),
                'n_usable_chemistry':len(eligible), 'n_excluded_chemistry':len(values)-len(eligible),
                'n_spatial':sum(value['spatial'] for value in values),
                'n_unlocated_usable':sum(value['eligible'] and not value['spatial'] for value in values),
                'n_semiquantitative_input':sum(value['semi'] for value in values),
                'n_semiquantitative_used':sum(value['semi'] for value in eligible),
                'n_method_unspecified_used':sum(value['method'] is None or value['method']=='UNKNOWN' for value in eligible),
                'n_method_not_reviewed_used':sum(value['method_information']=='method_not_reviewed' for value in eligible)}
            for field, value in expected.items():
                if field in row: assert numeric_equal(row[field], value), ('count',stage,medium,result_key(row),field,row[field],value)
    coverage = result_rows('method_coverage')
    for row in coverage:
        code = None if row['method_code'] == 'NA' else row['method_code']
        values = [value for value in groups[result_key(row)] if value['method']==code]
        assert int(row['n_input']) == len(values)
        assert int(row['n_spatial']) == sum(value['spatial'] for value in values)
        if int(row['n_spatial']) == 0: assert number(row['occupied_cells_mean']) is None
    pairs = result_rows('method_overlap')
    expected_pair_count = 0
    for group, values in groups.items():
        named = sorted({value['method'] for value in values if value['method'] not in {None,'UNKNOWN'}})
        named_pairs = list(itertools.combinations(named, 2));expected_pair_count += len(named_pairs)*len(config['cell_sizes_m'])
        if not named_pairs: continue
        points = [value for value in values if value['spatial']]
        anchor_x = min((value['x'] for value in points), default=None)
        anchor_y = min((value['y'] for value in points), default=None)
        targets = {(float(row['cell_size']),frozenset((row['method_a'],row['method_b']))):row for row in pairs if result_key(row)==group}
        for size in config['cell_sizes_m']:
            history = {pair: [] for pair in named_pairs}
            for ox in range(config['n_offsets']):
                for oy in range(config['n_offsets']):
                    cells = {code:set() for code in named}
                    for point in points:
                        if point['method'] not in cells: continue
                        cells[point['method']].add((math.floor((point['x']-anchor_x)/size+ox/config['n_offsets']),
                            math.floor((point['y']-anchor_y)/size+oy/config['n_offsets'])))
                    for pair in named_pairs:
                        a,b = cells[pair[0]], cells[pair[1]]
                        if a and b:
                            shared = len(a & b); union = len(a | b)
                            history[pair].append((shared,union,shared/union))
            for pair, trials in history.items():
                target = targets[(float(size),frozenset(pair))]
                assert int(target['n_origins_evaluated']) == len(trials)
                if trials:
                    for index, metric in enumerate(['shared_cells','union_cells','jaccard']):
                        values_metric = [trial[index] for trial in trials]
                        for suffix, expected in [('min',min(values_metric)),('mean',sum(values_metric)/len(values_metric)),('max',max(values_metric))]:
                            assert numeric_equal(target[metric+'_'+suffix],expected,1e-8), ('overlap',stage,medium,group,pair,size,metric,suffix)
                    assert numeric_equal(target['fraction_origins_with_overlap'],sum(trial[0]>0 for trial in trials)/len(trials),1e-8)
                else:
                    assert number(target['jaccard_mean']) is None and number(target['shared_cells_mean']) is None
                overlap_checks.append(dict(stage=stage,medium=medium,element=group[0],
                    collection=group[2] or '',method_a=pair[0],method_b=pair[1],cell_size=size,
                    n_evaluated=len(trials),independent_overlap=True))
    assert len(pairs) == expected_pair_count
    for row in result_rows('declustering'):
        assert number(row['ci_low']) is None and number(row['ci_high']) is None
    checks.append(dict(stage=stage,medium=medium,n_records=len(observed),n_groups=len(groups),
        n_pair_scale_comparisons=len(pairs),raw_fields_and_censoring_reconciled=True))
    print(stage, medium, 'observations/counts/overlap verified', len(observed), flush=True)

verification = out / 'verification';verification.mkdir(exist_ok=True)
def write(path, records):
    with path.open('w',newline='',encoding='utf-8') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(records[0]));writer.writeheader();writer.writerows(records)
write(verification/'source_reconciliation.csv', extraction_checks)
write(verification/'observation_reconciliation.csv', checks)
if overlap_checks: write(verification/'independent_overlap.csv', overlap_checks)
write(verification/'validation_status_counts.csv',[dict(stage=key[0],medium=key[1],status=key[2],n_records=value)
    for key,value in sorted(status_totals.items())])
for manifest in ['input_manifest.csv','source_manifest.csv']:
    for row in rows(out/manifest):assert digest(Path(row['file']))==row['sha256']
(verification/'PYTHON_VERIFIED.txt').write_text('Full-source raw field extraction, parsed observations/counts, method preservation and all geographic overlap comparisons reconcile. No spatial confidence validation is claimed.\n')
print('Independent streaming/numerical verification complete.',flush=True)
