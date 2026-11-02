#!/usr/bin/env python3
"""Freeze a fresh, method-preserving South Carolina validation input set."""
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path
import shutil
import subprocess
import sys

ROOT = Path.cwd()
OUT = Path(sys.argv[1] if len(sys.argv) > 1 else 'output/usgs_method_validation_2026-11-02').resolve()
if OUT.exists():
    raise SystemExit('Use a new output directory; previous runs are preserved')
OUT.mkdir(parents=True)

def fingerprints(path):
    sha = hashlib.sha256(); md5 = hashlib.md5()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            sha.update(block); md5.update(block)
    return {'file': str(path.resolve()), 'bytes': path.stat().st_size,
            'sha256': sha.hexdigest(), 'md5': md5.hexdigest()}

def write_csv(path, rows, fields):
    with path.open('w', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader(); writer.writerows(rows)

protocol = '''# Frozen observational validation protocol - November 2, 2026

Scope is declared before chemistry values are read: explicitly U.S.-identified
South Carolina (SC) sample metadata records, separately for soil and sediment.
Country/state labels are used as recorded after trimming/case normalization;
no jurisdiction is inferred from coordinates. Missing/blank source information
stays missing. No physical-site identity or homogeneous population is asserted.

Pilot: soil Au and Al in records explicitly labelled single/grab. These elements
were used in an earlier censoring diagnostic, so this is a deliberately selected
pipeline check, not independent statistical-performance validation. Expansion:
every nonblank reported species in the SC soil and sediment records, separately.
Analysis groups are state and recorded method_collected; blank collection labels
remain an NA bucket and literal unknown remains explicit. Common labels do not
prove comparable support. Depth, size fraction and collection protocol are
largely absent in these records; no values or eligibility are inferred.

Use the existing per-state Lambert azimuthal equal-area coordinates and recorded
CRS definitions. Match IDs and source coordinates; independently reproject finite
source coordinates and require agreement within 0.00001 m. Grid widths are fixed
at 5, 10 and 20 km with 10 offsets per axis (100 placements), as in the earlier
exploratory run. They are diagnostic scales, not optimized geological spacings.

Preserve every chemistry field: rec_id, lab_id, species, unit, qvalue, amethod,
summary, parameter and mthd_parameter, plus its source data-record ordinal.
Read identifiers and qvalue as text. Keep invalid/unsupported/all-left groups;
no missing sample-by-element results, methods, limits or class bins are created.

Run the checked geoCensor 1.1.0.9003 installed build. Retain unweighted,
method-stratified, declustered and geographic-overlap tables. Model means/ordinary
limits are conditional on recorded coding and unvalidated model/independence
assumptions; they are not validated geologic-population or spatial intervals.
Weighted routine limits stay NA. No p-values, method-bias estimates, automatic
method choice, spatial covariance calibration or reporting-class conversion.
Spatial confidence limits remain experimental.

Acceptance requires exact source-record reconciliation, method retention,
consistent input/usable/coordinate counts, explicit missingness and independent
checks of fits/weights/CDF bounds/overlap. The direct-fit/weighted-CDF panel is
fixed to Au, Al, Ag, Cu, Pb and Zn; CDF thresholds are 10^(-2), 10^(-1.5), ...,
10^9 ug/kg. Every geographic-overlap comparison is checked independently. Processing checkpoint reuse
requires unchanged inputs, source, protocol and runtime signatures. Preserve
all earlier outputs; no source data or prior study is overwritten.
'''
(OUT / 'PROTOCOL.md').write_text(protocol)
config = {'state': 'SC', 'media': ['soil', 'sed'], 'pilot_elements': ['Au', 'Al'],
          'pilot_collection': 'single/grab', 'cell_sizes_m': [5000, 10000, 20000],
          'n_offsets': 10, 'group_cols': ['state', 'method_collected'],
          'package_version': '1.1.0.9003',
          'library_dir': str((ROOT / 'output/batch_method_overlap_2026-11-02/geoCensor.Rcheck').resolve()),
          'git_head': subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()}
(OUT / 'config.json').write_text(json.dumps(config, indent=2) + '\n')
with (OUT / 'config.csv').open('w', newline='') as stream:
    writer = csv.writer(stream); writer.writerow(['key', 'value'])
    for key, value in config.items():
        writer.writerow([key, ','.join(map(str, value)) if isinstance(value, list) else value])

research = [ROOT / 'scripts/analysis' / name for name in [
    'prepare_usgs_method_validation.py', 'run_usgs_method_validation.R',
    'verify_usgs_method_validation.py', 'verify_usgs_method_validation.R']]
package_manifest = ROOT / 'output/batch_method_overlap_2026-11-02/package_source_manifest.csv'
source_rows = []
with package_manifest.open() as stream:
    for entry in csv.DictReader(stream):
        path = ROOT / entry['file']; stamp = fingerprints(path)
        if stamp['sha256'] != entry['sha256']:
            raise SystemExit('Package source differs from the checked build: ' + entry['file'])
        destination = OUT / 'source_snapshot/package' / entry['file']
        destination.parent.mkdir(parents=True, exist_ok=True); shutil.copy2(path, destination)
        source_rows.append(stamp)
for path in research:
    stamp = fingerprints(path); source_rows.append(stamp)
    destination = OUT / 'source_snapshot/research' / path.name
    destination.parent.mkdir(parents=True, exist_ok=True); shutil.copy2(path, destination)
write_csv(OUT / 'source_manifest.csv', source_rows, ['file', 'bytes', 'sha256', 'md5'])

inputs = []
for medium in config['media']:
    inputs.extend([ROOT / '../data' / ('ngdb' + medium) / name for name in ['main.csv', 'bestvalue.csv']])
    inputs.extend([ROOT / 'output/usgs_all_elements_2026-10-31' / medium / name
                   for name in ['metadata_projected.csv', 'projection_definitions.csv']])
input_rows = [fingerprints(path) for path in inputs]
write_csv(OUT / 'input_manifest.csv', input_rows, ['file', 'bytes', 'sha256', 'md5'])

us = {'US', 'USA', 'UNITED STATES', 'U.S.'}
for medium in config['media']:
    folder = OUT / 'inputs' / medium; folder.mkdir(parents=True)
    source = ROOT / '../data' / ('ngdb' + medium)
    metadata_scope = {}; selected = []; n_main = 0
    with (source / 'main.csv').open(newline='', encoding='utf-8-sig') as stream:
        reader = csv.DictReader(stream); metadata_fields = reader.fieldnames
        for row in reader:
            n_main += 1; key = row['lab_id'].strip()
            if not key or key in metadata_scope:
                raise SystemExit('Missing/duplicate normalized metadata identifier')
            is_us = row['country'].strip().upper() in us
            in_scope = is_us and row['state'].strip().upper() == config['state']
            metadata_scope[key] = (is_us, in_scope)
            if in_scope: selected.append(row)
    write_csv(folder / 'main_selected.csv', selected, metadata_fields)
    inventory = Counter(); species_counts = Counter(); methods = Counter()
    with (source / 'bestvalue.csv').open(newline='', encoding='utf-8-sig') as stream, \
            (folder / 'chemistry_selected.csv').open('w', newline='', encoding='utf-8') as destination:
        reader = csv.DictReader(stream); fields = reader.fieldnames
        required = {'rec_id', 'lab_id', 'species', 'unit', 'qvalue', 'amethod', 'summary', 'parameter', 'mthd_parameter'}
        if set(fields) != required:
            raise SystemExit('Unexpected chemistry schema; review rather than dropping columns')
        writer = csv.DictWriter(destination, fieldnames=fields + ['source_record_number']); writer.writeheader()
        for number, row in enumerate(reader, start=1):
            inventory['n_input_chemistry'] += 1
            key = row['lab_id'].strip(); scope = metadata_scope.get(key)
            if not key: inventory['n_blank_chemistry_id'] += 1
            if scope is None:
                inventory['n_unmatched_chemistry'] += 1
            elif scope[0]:
                inventory['n_explicit_us_chemistry'] += 1
            else:
                inventory['n_outside_us_scope_chemistry'] += 1
            if scope is None or not scope[1]: continue
            inventory['n_selected_chemistry'] += 1
            species_counts[row['species']] += 1
            methods[(row['species'], row['amethod'], row['unit'])] += 1
            writer.writerow(dict(row, source_record_number=number))
            if inventory['n_selected_chemistry'] % 10000 == 0:
                print(medium, 'selected', inventory['n_selected_chemistry'], flush=True)
    overview = dict(inventory, medium=medium, n_input_metadata=n_main, n_selected_metadata=len(selected),
                    n_explicit_us_metadata=sum(value[0] for value in metadata_scope.values()),
                    n_reported_species=len(species_counts))
    overview = {key: overview.get(key, 0) for key in ['medium', 'n_input_metadata', 'n_explicit_us_metadata',
        'n_selected_metadata', 'n_input_chemistry', 'n_explicit_us_chemistry', 'n_outside_us_scope_chemistry',
        'n_unmatched_chemistry', 'n_blank_chemistry_id', 'n_selected_chemistry', 'n_reported_species']}
    write_csv(folder / 'inventory.csv', [overview], list(overview))
    write_csv(folder / 'species_counts.csv', [{'species': key, 'n_records': value}
        for key, value in sorted(species_counts.items())], ['species', 'n_records'])
    write_csv(folder / 'method_counts.csv', [{'species': key[0], 'amethod_raw': key[1], 'unit_raw': key[2],
        'n_records': value} for key, value in sorted(methods.items())], ['species', 'amethod_raw', 'unit_raw', 'n_records'])
    print(json.dumps(overview), flush=True)
for stamp in input_rows:
    if fingerprints(Path(stamp['file'])) != stamp:
        raise SystemExit('Input changed during preparation: ' + stamp['file'])
for stamp in source_rows:
    if fingerprints(Path(stamp['file'])) != stamp:
        raise SystemExit('Source changed during preparation: ' + stamp['file'])
(OUT / 'PREPARED.txt').write_text('Method-preserving extraction complete; full input/source fingerprints unchanged.\n')
print('Prepared', OUT, flush=True)
