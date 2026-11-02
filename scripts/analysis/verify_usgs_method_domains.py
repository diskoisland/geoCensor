#!/usr/bin/env python3
"""Independently reconcile every domain count and selected metadata edge cases."""
import csv
import hashlib
from collections import Counter
from pathlib import Path
import sys

out=Path(sys.argv[1] if len(sys.argv)>1 else 'output/usgs_method_domains_2026-11-02')
checks=[];panels=[]
def csv_rows(path):
 with path.open(encoding='utf-8-sig',newline='') as f:yield from csv.DictReader(f)
def key(row):return tuple(row[n] for n in ('species','amethod','unit','project_name'))
def digest(path,kind='md5'):
 h=hashlib.new(kind)
 with path.open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()
for medium in ('soil','sed'):
 folder=out/medium
 assert (folder/'COMPLETED.txt').exists(),medium
 overview=next(csv_rows(folder/'overview.csv'))
 domains=list(csv_rows(folder/'method_project_domains.csv'))
 expected={key(r):int(r['n_records']) for r in domains}
 assert len(expected)==len(domains)
 assert all(r['p_value']=='NA' and r['reporting_design_status']=='unestablished' and
            r['physical_site_map_status']=='unestablished' for r in domains)
 assert overview['n_eligible_domains']=='NA' and overview['p_value']=='NA'
 fields=[c[len('n_blank_'):] for c in domains[0] if c.startswith('n_blank_') and c!='n_blank_state']
 selected={}
 def choose(label,rows,metric='n_records'):
  if rows:
   r=max(rows,key=lambda x:int(x[metric]));selected.setdefault(key(r),[]).append(label)
 choose('largest_observed_group',domains)
 choose('recorded_project',[r for r in domains if r['project_name']])
 choose('blank_project',[r for r in domains if not r['project_name']])
 choose('blank_method',[r for r in domains if not r['amethod']])
 choose('blank_unit',[r for r in domains if not r['unit']])
 choose('multiple_size_labels',domains,'n_distinct_nonblank_mesh_size')
 choose('literal_unknown_collection',domains,'n_literal_unknown_method_collected')
 choose('multiple_state_labels',domains,'n_nonblank_states')
 choose('reused_coordinates',domains,'n_repeated_coordinate_groups')
 choose('single_record',[r for r in domains if int(r['n_records'])==1])
 choose('all_collection_dates_recorded',[r for r in domains if r['n_blank_date_collect']=='0'])
 if medium=='soil':choose('multiple_horizon_zones',domains,'n_distinct_nonblank_sample_zone')
 expected_panel={key(r):r for r in domains if key(r) in selected}
 data={k:dict(count=0,ids=set(),coordinates=Counter(),states=set(),blank_state=0,resubmissions=0,
                blanks=Counter(),unknowns=Counter(),levels={f:set() for f in fields}) for k in selected}
 meta={};country_counts=Counter()
 source=Path('../data')/('ngdb'+medium)
 for row in csv_rows(source/'main.csv'):
  sample=row['lab_id'].strip();assert sample and sample not in meta
  country=row['country'].strip().upper();is_us=country in ('US','USA','UNITED STATES','U.S.')
  country_counts[country]+=1
  vals={f:row[f].strip() for f in fields}
  try:
   lat,lon=float(row['lat_wgs84']),float(row['long_wgs84'])
   coordinate=(lat,lon) if -90<=lat<=90 and -180<=lon<=180 else None
  except ValueError:coordinate=None
  meta[sample]=(is_us,row['project_name'].strip(),coordinate,row['state'].strip(),
                row['previous_lab_id'].strip(),vals)
 assert len(meta)==int(overview['n_input_metadata'])
 assert sum(v[0] for v in meta.values())==int(overview['n_us_metadata'])
 observed=Counter();totals=Counter();species_counts=Counter()
 for row in csv_rows(source/'bestvalue.csv'):
  totals['input']+=1
  sample=row['lab_id'].strip();m=meta.get(sample)
  if not sample:totals['missing_id']+=1
  if m is None:totals['unmatched']+=1;continue
  if not m[0]:totals['outside_us_scope']+=1;continue
  totals['retained']+=1
  k=(row['species'].strip(),row['amethod'].strip(),row['unit'].strip(),m[1])
  observed[k]+=1;species_counts[k[0]]+=1
  if k not in data:continue
  p=data[k];p['count']+=1;p['ids'].add(sample)
  if m[2] is not None:p['coordinates'][m[2]]+=1
  if m[3]:p['states'].add(m[3])
  else:p['blank_state']+=1
  p['resubmissions']+=bool(m[4])
  for f,v in m[5].items():
   p['blanks'][f]+=not bool(v);p['unknowns'][f]+=v.lower()=='unknown'
   if v:p['levels'][f].add(v)
 assert dict(observed)==expected,(medium,'domain counts differ')
 for field,count in [('n_input_chemistry','input'),('n_missing_chemistry_ids','missing_id'),
                     ('n_unmatched_chemistry','unmatched'),('n_matched_non_us_chemistry','outside_us_scope'),
                     ('n_us_chemistry','retained')]:assert int(overview[field])==totals[count],(medium,field)
 assert len(species_counts)==int(overview['n_species'])
 assert {r['species']:int(r['n_records']) for r in csv_rows(folder/'species_summary.csv')}==dict(species_counts)
 prior=next(csv_rows(out/('prior_'+medium+'_inventory.csv')))
 for field,old in [('n_input_metadata','input_metadata_rows'),('n_us_metadata','us_metadata_rows'),
                   ('n_input_chemistry','input_chemistry_rows'),('n_us_chemistry','us_chemistry_rows'),
                   ('n_species','n_species')]:assert int(overview[field])==int(prior[old]),(medium,field)
 for k,p in data.items():
  r=expected_panel[k];loc=p['coordinates']
  metrics={'n_records':p['count'],'n_sample_ids':len(p['ids']),
   'n_valid_coordinate_records':sum(loc.values()),'n_exact_coordinate_groups':len(loc),
   'n_repeated_coordinate_groups':sum(v>1 for v in loc.values()),
   'n_records_at_repeated_coordinates':sum(v for v in loc.values() if v>1),
   'n_nonblank_states':len(p['states']),'n_blank_state':p['blank_state'],'n_resubmission_links':p['resubmissions']}
  for f in fields:
   metrics['n_blank_'+f]=p['blanks'][f]
   metrics['n_literal_unknown_'+f]=p['unknowns'][f]
   metrics['n_distinct_nonblank_'+f]=len(p['levels'][f])
  assert all(int(r[f])==v for f,v in metrics.items()),(medium,k,'metadata profile differs')
  panels.append(dict(medium=medium,species=k[0],amethod=k[1],unit=k[2],project_name=k[3],
                     selection_reasons=';'.join(selected[k]),n_records=p['count'],metrics_checked=len(metrics),passed=True))
 for row in csv_rows(folder/'source_input_manifest.csv'):assert digest(Path(row['file']))==row['md5']
 checks.append(dict(medium=medium,n_domain_counts_checked=len(expected),n_records_checked=totals['input'],
                    n_us_records=totals['retained'],n_species=len(species_counts),
                    all_domain_counts_match=True,prior_run_totals_match=True,
                    unknown_country_metadata_rows=country_counts[''],p_values_remain_NA=True))
 print(medium,'all',len(expected),'group counts and',len(data),'metadata edge cases match.',flush=True)
for name,rows in [('independent_verification.csv',checks),('verification_panel.csv',panels)]:
 with (out/name).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
(out/'verify_usgs_method_domains.py').write_bytes(Path(__file__).read_bytes())
(out/'VERIFIED.txt').write_text('Every method/unit/project count independently reconstructed from source CSVs.\nAll prior all-element run input/U.S./species totals reconcile.\nMetadata-only edge-case panel verifies support profiles and coordinate reuse.\nSource/input fingerprints match; all domain p-values and eligible-domain counts remain NA.\n')
