#!/usr/bin/env python3
"""Report inventories; link observed method codes to bundled descriptions only."""
from pathlib import Path
from collections import defaultdict,Counter
import csv,hashlib,re,sys
out=Path(sys.argv[1] if len(sys.argv)>1 else 'output/usgs_method_domains_2026-11-02')
assert (out/'VERIFIED.txt').exists()
def read(path):
 with path.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def write(name,rows):
 with (out/name).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
registry=[];method_totals=[];overview=[];source_manifest=[]
for medium in ('soil','sed'):
 path=Path('../data')/('ngdb'+medium)/('ngdb'+medium+'.met')
 text=path.read_text(encoding='utf-8-sig')
 definitions=defaultdict(list)
 blocks=re.split(r'(?m)^      Attribute:\s*$',text)
 for block in blocks:
  if not re.search(r'(?m)^        Attribute_Label:\s*amethod\s*$',block):continue
  for enum in re.split(r'(?m)^\s*Enumerated_Domain:\s*$',block)[1:]:
   code=re.search(r'(?m)^\s*Enumerated_Domain_Value:\s*(\S+)\s*$',enum)
   definition=re.search(r'Enumerated_Domain_Value_Definition:\s*(.*?)\n\s*Enumerated_Domain_Value_Definition_Source:',enum,re.S)
   source=re.search(r'(?m)^\s*Enumerated_Domain_Value_Definition_Source:\s*(\S+)',enum)
   if code and definition:
    definitions[code[1]].append((' '.join(definition[1].split()),source[1] if source else ''))
 methods=read(out/medium/'methods_by_species.csv')
 counts=Counter()
 for r in methods:counts[r['amethod']]+=int(r['n_records'])
 semi_records=0;documented_records=0
 for code,n in sorted(counts.items()):
  entries=definitions.get(code,[])
  descriptions=sorted({d for d,_ in entries})
  sources=sorted({s for _,s in entries})
  status='not_documented_in_bundled_amethod_list' if not entries else ('conflicting_descriptions' if len(descriptions)>1 else 'documented_description')
  explicit_semi=bool(len(descriptions)==1 and re.search(r'semi[- ]?quantitative',descriptions[0],re.I))
  if explicit_semi:semi_records+=n
  if status=='documented_description':documented_records+=n
  registry.append(dict(medium=medium,amethod=code,n_records=n,documentation_status=status,
   n_dictionary_occurrences=len(entries),n_distinct_descriptions=len(descriptions),
   bundled_description=' | '.join(descriptions),definition_source=' | '.join(sources),
   documented_semiquantitative='TRUE' if explicit_semi else 'NA',
   class_bounds_status='not_established',reporting_design_status='not_established'))
  method_totals.append(dict(medium=medium,amethod=code,n_records=n,documentation_status=status,
   semiquantitative_review_flag=explicit_semi))
 base=read(out/medium/'overview.csv')[0]
 base.update(n_records_documented_method_description=documented_records,
  n_records_documented_semiquantitative=semi_records,
  semiquantitative_record_fraction=semi_records/int(base['n_us_chemistry']))
 overview.append(base)
 source_manifest.append(dict(file=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
 for name in ['methods_by_species.csv','overview.csv']:
  p=out/medium/name;source_manifest.append(dict(file=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
write('method_documentation_registry.csv',registry)
write('method_record_totals.csv',method_totals)
write('overview.csv',overview)
source_manifest.append(dict(file=str(Path(__file__)),sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()))
write('method_documentation_manifest.csv',source_manifest)
(out/'report_usgs_method_domains.py').write_bytes(Path(__file__).read_bytes())
summary_rows=[]
for r in overview:
 summary_rows.append(f"| {r['medium']} | {int(r['n_us_metadata']):,} | {int(r['n_us_chemistry']):,} | {r['n_species']} | {int(r['n_method_unit_project_groups']):,} | {int(r['n_records_blank_project']):,} | {int(r['n_records_documented_semiquantitative']):,} ({100*float(r['semiquantitative_record_fraction']):.2f}%) |")
report='''# USGS method/project inventory — November 2, 2026

Both best-value chemistry tables and both main sample-metadata tables are now
readable. Longer bounded reads restored the cloud-only chemistry files after
short access probes timed out. Earlier timeout snapshots remain preserved.
This audit inventories identifiers/methods/support only; no concentration or
qualifier values were loaded for grouping, fitted or used to select domains.

| Medium | U.S. sample metadata | U.S. best-value records | Species | Method/unit/project groups | Records without project | Documented semiquantitative records |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
'''+ '\n'.join(summary_rows)+'''

These are observed sample/element result records, not independent samples or
sites. Species include reported properties/auxiliary measurements. No missing
sample-by-element results were constructed. Coordinate groups remain exact
position matches, not established physical sites. All observed U.S.-matched
species, method/unit/project groups and blank labels are retained.

## What the grouping does and does not establish

Groups use raw trimmed species, analytical-method, unit and project labels.
A common project/method label does not establish common reporting limits,
physical support, geological population or exchangeability. The per-group table
reports support/date/geological label gaps and mixtures without choosing a
preferred domain, merging methods, filling labels or imposing a new threshold.
Blank project groups are missing identifiers, not inferred shared projects.
Nonblank dates/support labels are not automatically valid or comparable values.
The token "unknown" is counted separately; other nonblank labels still require
review. Methods and concentration units have not been harmonized here.

`n_matched_non_us_chemistry` means outside the explicitly U.S.-identified scope.
It includes other geographic labels such as oceans and the Gulf of Mexico;
it is not a claim that every such location is foreign. No country is inferred
from coordinates. In these two metadata tables, country labels are nonblank.
All chemistry records match a metadata identifier and no chemistry IDs are blank.

The soil audit yields 10,812 groups and sediment 20,701: 31,513 descriptive
groups, not 31,513 validated statistical domains. The reporting-design and
physical-site-map statuses are unestablished, and every group p-value and the
number of eligible domains remain NA. Bounds encoded for censored records do
not establish both reporting-design endpoints for exact records. No spatial
range, concentration estimate, prospectivity score or confidence interval was
computed. Confidence limits remain experimental.

## Semiquantitative method review is a concrete next requirement

The bundled soil and sediment metadata explicitly define ES_SQ as
semiquantitative emission spectrography. The method documentation registry
records those definitions and their sources. It marks a semiquantitative
review flag only where the supplied description says so; other methods are
not thereby certified quantitative. Absent/conflicting descriptions remain
unresolved. Counts include all observed records under that code, not just
uncensored results or exact-fit-eligible observations.

These flags do not supply numeric class bounds. A method family label does
not identify a unique step table, instrument/procedure variant or laboratory/
era-specific reporting design. The public parser currently handles censoring
encodings; it does not decode semiquantitative reporting classes into finite
intervals. This audit does not silently change their concentrations, convert
every ES_SQ record into a guessed interval or recompute previous estimates.
A reported number lacking a censoring qualifier need not represent an exact
continuous measurement under a semiquantitative reporting procedure.

USGS's [AGDB method metadata](https://pubs.usgs.gov/ds/759/contents/AGDB2_DS759_metadata/AGDB2_MD.faq.html)
also describes ES_SQ variants. The [2020 USGS Alaska methods evaluation](https://pubs.usgs.gov/of/2020/1038/ofr20201038.pdf)
reports that comparison scatter can exceed nominal semiquantitative reporting
ranges in that study. We consulted online metadata and abstract/index material,
not the complete PDF. That is a reason for method review, not a quantified NGDB
bias estimate or authority to turn reporting bins into confidence limits.

## Validation and retained evidence

The independent Python verifier reconstructs every method/unit/project count
from all 14,247,668 input best-value records. All 13,792,077 retained U.S.
records, 80/81 species and metadata totals reconcile with the prior all-element
run. Fourteen metadata edge-case groups independently check support profiles,
state mixtures and coordinate reuse. Selection used metadata/count cases,
not concentration outcomes. All input/source fingerprints match.

The first R attempt stopped before domain outputs on a bare logical-column
selection in data.table. The original source/signature and failure note are
preserved in development_snapshot; the explicit row condition passed. No
statistical response was used to tune the audit. Completed-mode resume and
artifact fingerprints provide additional integrity checks.

- [Combined totals](overview.csv)
- [Independent all-group reconciliation](independent_verification.csv)
- [Metadata edge-case checks](verification_panel.csv)
- [Method definitions and unresolved mappings](method_documentation_registry.csv)
- [Method record counts](method_record_totals.csv)
- [Soil groups](soil/method_project_domains.csv) and [sediment groups](sed/method_project_domains.csv)
- [Soil species summary](soil/species_summary.csv) and [sediment species summary](sed/species_summary.csv)
- [Frozen audit protocol](PROTOCOL.md)

Next: use this inventory to establish traceable method/reporting-class sources
and physical-support/site mappings. Verify applicability before using any
numeric class or detection-limit table. Keep unknown bounds and eligibility
missing. This work is separate from the published v1.1.0 API and release tag.
'''
(out/'README.md').write_text(report)
print('Method descriptions linked; report and registry written.')
for r in overview:print(r['medium'],'semiquantitative review records',r['n_records_documented_semiquantitative'],'share',r['semiquantitative_record_fraction'])
