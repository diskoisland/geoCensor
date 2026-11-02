#!/usr/bin/env python3
"""Read-only schema/availability audit; never fit or test a USGS concentration."""
import argparse
import configparser
import csv
import hashlib
import io
import json
import re
import subprocess
import sys
from collections import Counter
from pathlib import Path

TABLES = ('main', 'bestvalue', 'chemistry')


def read_bounded(path, *, header=False, timeout=5):
    """A cloud read failure means unavailable, never an empty data table."""
    script = ('import sys,csv,json;f=open(sys.argv[1],encoding="utf-8-sig");'
              'print(json.dumps(next(csv.reader(f))))' if header else
              'import sys;print(open(sys.argv[1],encoding="utf-8-sig").read(),end="")')
    try:
        result = subprocess.run([sys.executable, '-c', script, str(path)],
                                capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        return 'read_timeout', None
    if result.returncode:
        return 'read_error', None
    return 'readable', json.loads(result.stdout) if header else result.stdout


def schema_fields(text):
    parser = configparser.ConfigParser(interpolation=None)
    parser.read_string(text)
    result = {}
    for table in TABLES:
        section = table + '.csv'
        fields = []
        for key, value in parser.items(section):
            match = re.fullmatch(r'col(\d+)', key)
            if match:
                fields.append((int(match[1]), value.split()[0]))
        fields.sort()
        assert [n for n, _ in fields] == list(range(1, len(fields)+1))
        result[table] = [name for _, name in fields]
    return result


def write_csv(path, rows, fields):
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def fingerprint(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(root, out, timeout):
    out.mkdir(parents=True, exist_ok=False)
    (out/'source_snapshot').mkdir()
    source = Path(__file__)
    (out/'source_snapshot'/source.name).write_bytes(source.read_bytes())
    availability, columns, requirements, agreement, fingerprints = [], [], [], [], []
    metadata_summaries, metadata_coverage, locations = [], [], []
    for medium in ('soil', 'sed'):
        folder = root/('ngdb'+medium)
        texts = {}
        for name in ('schema.ini', 'datadxny.csv'):
            path = folder/name
            status, content = read_bounded(path, timeout=timeout)
            assert status == 'readable', (path, status)
            texts[name] = content
            target = out/(medium+'_'+name)
            target.write_text(content)
            fingerprints.append(dict(file=str(path), sha256=fingerprint(path)))
            availability.append(dict(medium=medium, file=str(path), status=status,
                                     bytes=path.stat().st_size, evidence='full documentation read'))
        schema = schema_fields(texts['schema.ini'])
        dictionary = list(csv.DictReader(io.StringIO(texts['datadxny.csv'])))
        for table in TABLES:
            declared = [r for r in dictionary if r['tblname']==table and r['label']]
            assert [r['label'] for r in declared] == schema[table], (medium, table)
            path = folder/(table+'.csv')
            status, header = read_bounded(path, header=True, timeout=timeout)
            availability.append(dict(medium=medium, file=str(path), status=status,
                                     bytes=path.stat().st_size,
                                     evidence='CSV header read' if header else 'record content unavailable; not zero rows'))
            if header is not None:
                assert header == schema[table], (medium, table, 'actual/schema mismatch')
            agreement.append(dict(medium=medium, table=table, n_declared_fields=len(declared),
                                  dictionary_matches_schema=True,
                                  actual_header_matches='TRUE' if header is not None else 'NA'))
            columns.extend(dict(medium=medium, table=table, field=r['label'],
                                description=r['description'].replace('\n',' '),
                                evidence='bundled dictionary and schema agree') for r in declared)
        def req(name, table, fields, status, meaning):
            assert all(f in schema[table] for f in fields)
            requirements.append(dict(medium=medium, requirement=name, table=table,
                candidate_fields=';'.join(fields), documentation_status=status,
                record_completeness='NA', scientific_eligibility='NA', interpretation=meaning))
        req('element_method_unit','bestvalue',['species','amethod','unit'],'declared',
            'Method/unit labels do not establish equivalent digestion, reporting design or exchangeability.')
        req('observed_censoring_encoding','chemistry',['qvalue','qualifier'],'declared',
            'A censored record encodes one observed bound; exact values do not encode both design limits.')
        req('best_value_selection','bestvalue',['summary','rec_id'],'declared',
            'Summary lists related analyses; rec_id is table-local, not a documented cross-table join key.')
        req('complete_lower_reporting_design','bestvalue',[],'not_explicitly_declared',
            'No dedicated lower reporting-limit field for every record, including exacts; method-specific external lookup not established.')
        req('complete_upper_reporting_design','bestvalue',[],'not_explicitly_declared',
            'No dedicated upper reporting-limit field for every record; an unbounded endpoint needs documented design justification.')
        req('reporting_limit_provenance','chemistry',[],'not_explicitly_declared',
            'No dedicated record-level limit-source field; possible external method documentation has not been audited.')
        req('laboratory_identity','main',['lab_id','job_id'],'not_established',
            'lab_id identifies a submitted sample and job_id a batch; neither is documented as an analytical laboratory identity.')
        req('sample_support','main',['primary_class','mesh_size','depth','prep','collection_protocol'],'declared',
            'Labels and nonmissing entries require domain review; missing labels cannot be interpreted as common support.')
        req('site_identity','main',['lab_id','field_id','previous_lab_id','previous_job_id'],'not_established',
            'Sample and historical IDs aid review but do not establish a physical-site/repeat map; field_id can be truncated.')
        req('coordinates_and_precision','main',['lat_wgs84','long_wgs84','coord_prec_est'],'declared',
            'Identical coordinates do not prove identical physical support; rounding/precision needs review.')
        req('survey_domain_and_dates','main',['project_name','date_collect','date_submitted'],'declared',
            'Project/date labels alone do not establish representative sampling or common reporting limits.')
        req('exchangeability_under_null','main',[],'not_established',
            'Requires substantive geology, mean/trend and sampling-design review; no software assertion or matching labels prove it.')
        req('balanced_site_vectors','main',[],'not_established',
            'Prototype permits 4–128 distinct sites and 1–4 ordered records/site; no thinning or arbitrary repeat pairing is authorized.')
        req('physical_distance_scales','main',['lat_wgs84','long_wgs84'],'not_established',
            'Simulation scales .05/.15/.4 are unit-square values; field scales and metric projection need prespecification.')
        req('across_domain_multiplicity','bestvalue',[],'not_established',
            'Single-domain validation does not justify automatically testing all elements, states or method subsets.')
        main_path = folder/'main.csv'
        main_status = next(r['status'] for r in availability if r['file']==str(main_path))
        # Header success materializes the file before this ordinary complete read.
        if main_status == 'readable':
            counts=Counter(); coverage=Counter(); coordinate_counts=Counter()
            support=['lab_id','field_id','job_id','date_collect','date_submitted','project_name',
                     'coord_prec_est','depth','mesh_size','prep','collection_protocol','previous_lab_id']
            if medium=='soil': support.append('horizon')
            with main_path.open(encoding='utf-8-sig',newline='') as handle:
                for row in csv.DictReader(handle):
                    counts['input_rows']+=1
                    if row['country'].strip().upper() not in ('US','USA','UNITED STATES','U.S.'):continue
                    counts['us_rows']+=1
                    for field in support:
                        if row[field].strip():coverage[field]+=1
                    try:
                        lat,lon=float(row['lat_wgs84']),float(row['long_wgs84'])
                        valid=(-90<=lat<=90 and -180<=lon<=180)
                    except ValueError:valid=False
                    if valid:coordinate_counts[(lat,lon)]+=1;counts['valid_coordinate_rows']+=1
            n=counts['us_rows']
            metadata_summaries.append(dict(medium=medium,status='available',**counts,
                exact_coordinate_groups=len(coordinate_counts),
                repeated_coordinate_groups=sum(v>1 for v in coordinate_counts.values()),
                records_at_repeated_coordinates=sum(v for v in coordinate_counts.values() if v>1),
                established_physical_sites='NA',eligible_tests='NA',p_value='NA'))
            metadata_coverage.extend(dict(medium=medium,field=f,n_us_records=n,
                n_nonblank=coverage[f],fraction_nonblank=coverage[f]/n if n else 'NA',
                interpretation='Nonblank raw labels, not verified usable/comparable values') for f in support)
            locations.extend(dict(medium=medium,records_at_coordinate=k,n_coordinate_groups=v,
                interpretation='Exact coordinate matches only; not established physical sites')
                for k,v in sorted(Counter(coordinate_counts.values()).items()))
            fingerprints.append(dict(file=str(main_path),sha256=fingerprint(main_path)))
        else:
            metadata_summaries.append(dict(medium=medium,status=main_status,
                input_rows='NA',us_rows='NA',valid_coordinate_rows='NA',exact_coordinate_groups='NA',
                repeated_coordinate_groups='NA',records_at_repeated_coordinates='NA',
                established_physical_sites='NA',eligible_tests='NA',p_value='NA'))
    write_csv(out/'input_availability.csv',availability,['medium','file','status','bytes','evidence'])
    write_csv(out/'schema_columns.csv',columns,['medium','table','field','description','evidence'])
    write_csv(out/'schema_agreement.csv',agreement,['medium','table','n_declared_fields','dictionary_matches_schema','actual_header_matches'])
    write_csv(out/'eligibility_requirements.csv',requirements,['medium','requirement','table','candidate_fields','documentation_status','record_completeness','scientific_eligibility','interpretation'])
    write_csv(out/'metadata_summary.csv',metadata_summaries,['medium','status','input_rows','us_rows','valid_coordinate_rows','exact_coordinate_groups','repeated_coordinate_groups','records_at_repeated_coordinates','established_physical_sites','eligible_tests','p_value'])
    write_csv(out/'metadata_field_coverage.csv',metadata_coverage,['medium','field','n_us_records','n_nonblank','fraction_nonblank','interpretation'])
    write_csv(out/'coordinate_multiplicity.csv',locations,['medium','records_at_coordinate','n_coordinate_groups','interpretation'])
    fingerprints.append(dict(file=str(source),sha256=fingerprint(source)))
    write_csv(out/'source_input_manifest.csv',fingerprints,['file','sha256'])
    (out/'runtime.json').write_text(json.dumps(dict(python=sys.version,argv=sys.argv),indent=2)+'\n')
    (out/'COMPLETED.txt').write_text('Documentation/availability audit complete; record-level eligibility remains unestablished.\nNo p-value calculated; unavailable is NA, not zero sites or no spatial structure.\n')
    print('Audit complete:',len(columns),'documented fields;',len(requirements),'requirements;',
          sum(r['status']=='readable' and r['evidence']=='CSV header read' for r in availability),'CSV headers readable.')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-root',type=Path,default=Path('../data'))
    parser.add_argument('--out',type=Path,default=Path('output/usgs_spatial_order_eligibility_2026-11-01'))
    parser.add_argument('--timeout',type=float,default=10)
    args=parser.parse_args()
    run(args.data_root,args.out,args.timeout)
