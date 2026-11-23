#!/usr/bin/env python3
"""Reproduce CN4461 evidence audit and three explicitly provisional pilot designs.

Requires Biopython 1.85. Coordinates in exported tables are 1-based inclusive.
No sequence is fetched, mutated, ordered, or assigned a new barcode by this script.
"""
import copy
import csv
import gzip
import hashlib
import json
from difflib import SequenceMatcher
from pathlib import Path
from zipfile import ZipFile
import xml.etree.ElementTree as ET

from Bio import SeqIO
from Bio import __version__ as biopython_version
from Bio.Restriction import MluI, SacI
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import MeltingTemp

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'designs/CN4461_AAV_MPRA'
REF = BASE / 'reference'
OUT = BASE / 'verification'
PILOT = OUT / 'pilot_review_only'
for p in (OUT, PILOT):
    p.mkdir(parents=True, exist_ok=True)


def sha(s):
    return hashlib.sha256(s.encode() if isinstance(s, str) else s).hexdigest()


def save_json(path, data):
    path.write_text(json.dumps(data, indent=2) + '\n')


def write_tsv(path, rows):
    with path.open('w') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t')
        w.writeheader()
        w.writerows(rows)


def hits(s, motif, circular=False):
    hay = s + s[:len(motif)-1] if circular else s
    return [i+1 for i in range(len(s)) if hay.startswith(motif, i)]


def longest_shared_run(query, subject):
    """Exact substring comparison on both orientations, not a folding model."""
    subjects = [subject, str(Seq(subject).reverse_complement())]
    for size in range(len(query), 0, -1):
        for offset in range(len(query)-size+1):
            word = query[offset:offset+size]
            if any(word in target for target in subjects):
                return dict(length=size, sequence=word)
    return dict(length=0, sequence='')


def fasta(path, records):
    SeqIO.write([SeqRecord(Seq(s), id=n, description='') for n, s in records], path, 'fasta')


def supplementary_row():
    ns = {'s': 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
    path = next(REF.glob('*Table_S1*.xlsx'))
    with ZipFile(path) as z:
        strings = [''.join(x.itertext()) for x in ET.fromstring(z.read('xl/sharedStrings.xml'))]
        sheet = ET.fromstring(z.read('xl/worksheets/sheet1.xml'))
        for row in sheet.findall('.//s:row', ns):
            cells = {}
            for c in row:
                v = c.find('s:v', ns)
                value = v.text if v is not None else ''
                if c.get('t') == 's':
                    value = strings[int(value)]
                cells[c.get('r').rstrip('0123456789')] = value
            if cells.get('B') == 'AiE0779m_3xC2':
                return {'sheet': 'sheet1', 'row': int(row.get('r')), 'cells': cells}
    raise ValueError('Published enhancer row missing')


a = SeqIO.read(REF / 'addgene-plasmid-236600-sequence-472831.gbk', 'genbank')
ad = SeqIO.read(REF / 'addgene-plasmid-236600-sequence-472831.dna', 'snapgene')
b = SeqIO.read(REF / 'addgene-plasmid-191706-sequence-379462.gbk', 'genbank')
bd = SeqIO.read(REF / 'addgene-plasmid-191706-sequence-379462.dna', 'snapgene')
s, t = str(a.seq).upper(), str(b.seq).upper()
assert s == str(ad.seq).upper() and t == str(bd.seq).upper()
assert len(s) == 4820 and len(t) == 4808
assert s[:1827] + s[1839:] == t
assert hits(s, 'ACGCGT', True) == [139]
assert hits(s, 'GAGCTC', True) == [727]
assert MluI.search(Seq(s), linear=False) == [140]
assert SacI.search(Seq(s), linear=False) == [732]
published = supplementary_row()
assert published['cells']['I'].upper() == s[144:726]
assert s[144:338] == s[338:532] == s[532:726]
save_json(OUT / 'published_enhancer_evidence.json', published)
fasta(OUT / 'published_enhancer_and_repeat.fasta', [('AiE0779m_3xC2_145_726', s[144:726]), ('C2_repeat_194bp', s[144:338])])

blocks = []
for tag, i, j, k, l in SequenceMatcher(None, s, t, autojunk=False).get_opcodes():
    blocks.append({'operation_CN4461_to_CN3044': tag, 'CN4461_start': i+1,
                   'CN4461_end': j, 'CN3044_start': k+1, 'CN3044_end': l,
                   'CN4461_difference': s[i:j] if tag != 'equal' else '',
                   'CN3044_difference': t[k:l] if tag != 'equal' else ''})
write_tsv(OUT / 'CN4461_CN3044_alignment.tsv', blocks)

# Exact WPRE-related match to the primary Choi CW3B accession, without
# treating an unannotated 247th terminal base as an established mutation.
choi = next(r for r in SeqIO.parse(REF / 'Choi2014_KJ411911_KJ411919.gb', 'genbank') if r.id == 'KJ411912.1')
assert s[1579:1825] == str(choi.seq[1300:1546]).upper()
assert s[1831:1839] == 'AGCAATTA'
assert s[1839:1854] == 'CGAGAGATCTTCGAC'

specs = [
    ('Left AAV2 ITR', 1, 130, '+ (map); inverted repeat', 'source_annotated', 'CN4461 and CN3044 native feature; alternate 130-bp ITR, not assumed wild-type 145 bp'),
    ('AiE0779m_3xC2', 145, 726, '+', 'confirmed_exact_published_match', 'Hunker Table S1 row 30: exact 582/582 match; three exact 194-bp repeats; CN3044 identical'),
    ('minBG candidate', 733, 785, '+', 'strong_inference_boundary_unconfirmed', 'Table S1 row 30 specifies Beta globin 53 bp; conserved 53 bp after SacI before BamHI; no independent base-resolved minBG annotation'),
    ('Promoter/linker protected interval', 733, 825, '+', 'confirmed_preserved_interval', 'Entire conserved sequence between SacI and SYFP2 start retained, irrespective of precise promoter boundary'),
    ('SYFP2 CDS', 826, 1542, '+', 'source_annotated_translation_checked', 'Native joined CDS is contiguous 717 nt; translation checked separately; following TAA stop retained'),
    ('SYFP2 stop', 1543, 1545, '+', 'confirmed_sequence', 'TAA immediately after annotated CDS'),
    ('WPRE3 aligned core', 1580, 1825, '+', 'confirmed_exact_core_terminal_boundary_unresolved', '246/246 exact match to Choi KJ411912.1 bases 1301-1546; paper calls WPRE3 247 bp; terminal inclusion unresolved'),
    ('WPRE/linker protected interval', 1546, 1831, '+', 'confirmed_preserved_interval', 'All post-CDS bases through upstream barcode flank retained; includes possible WPRE3 terminal base 1826'),
    ('BarcodeS candidate', 1832, 1839, '+', 'strong_inference_published_readout_pattern', 'AGCAATTA; 8 bp immediately before published (BC)CGAGAGATCTTCGAC pattern; only CN3044/CN4461 difference is at this junction; depositor barcode-S name mapping not independently supplied'),
    ('Barcode junction protected interval', 1826, 1854, '+', 'confirmed_preserved_interval', 'Entire CTCGAG-AGCAATTA-CGAGAGATCTTCGAC architecture retained; alignment insertion is 12 bp, not a 12-bp barcode; terminal C overlaps annotated BGHpA'),
    ('BGHpA annotated signal', 1854, 1965, '+', 'source_annotated', 'Native bGH poly(A) signal; sequence conserved with CN3044 at 1842-1953; not asserted to delimit entire processing region'),
    ('Downstream polyA context', 1966, 2082, '+', 'confirmed_preserved_interval', 'All downstream sequence through start of right ITR retained'),
    ('Right AAV2 ITR long annotation', 2083, 2223, '+ (map); inverted repeat', 'source_annotated_boundary_ambiguity', '141-bp native annotation overlaps 130-bp annotation; these are two annotations of ONE right ITR'),
    ('Right AAV2 ITR short annotation', 2083, 2212, '+ (map); inverted repeat', 'source_annotated_boundary_ambiguity', 'Alternative 130-bp native annotation; physical repeat endpoint/function not resolved by map alone'),
]
features = []
for name, start, end, strand, confidence, evidence in specs:
    seq = s[start-1:end]
    features.append(dict(feature=name, reference_start=start, reference_end=end, length=len(seq), strand=strand,
                         pilot_500bp_start=start if end <= 144 else (145 if start == 145 else start-82),
                         pilot_500bp_end=end if end <= 144 else (644 if start == 145 else end-82),
                         confidence=confidence, evidence=evidence, sequence=seq, sequence_sha256=sha(seq)))
write_tsv(OUT / 'feature_evidence.tsv', features)
fasta(OUT / 'feature_evidence.fasta', [(f['feature'].replace(' ', '_'), f['sequence']) for f in features])
cds = s[825:1542]
assert len(cds) == 717 and '*' not in str(Seq(cds).translate()) and s[1542:1545] == 'TAA'
native_translation = next(f for f in a.features if f.qualifiers.get('label') == ['SYFP2']).qualifiers['translation'][0].replace(',', '')
assert str(Seq(cds).translate()) == native_translation

# Synthesis flanks include restored recognition sites. The five
# enhancer-facing bases absent from each retained top strand are supplied by
# the synthetic insert. Left flank is kept short to reduce ITR homology:
# 15 bp are on the retained left top strand, 25 bp on the right.
left, right = s[124:144], s[726:756]
left_overlap_len = len(left)-5
right_overlap_len = len(right)-5
arms = []
backbone = s[731:] + s[:139]  # actual retained top strand after139 and after731
for name, seq, start, end, duplex in [('left', left, 125, 144, left[:left_overlap_len]), ('right', right, 727, 756, right[5:])]:
    arms.append(dict(arm=name, reference_start=start, reference_end=end, sequence_5to3=seq,
                     length=len(seq), retained_top_strand_overlap=duplex, retained_top_strand_overlap_length=len(duplex),
                     reference_exact_forward_hits=hits(s, seq, True), reference_exact_reverse_hits=hits(s, str(Seq(seq).reverse_complement()), True),
                     overlap_Tm_NN_C=round(MeltingTemp.Tm_NN(duplex), 2),
                     longest_exact_run_elsewhere=longest_shared_run(seq, s[:start-1]+'N'*len(seq)+s[end:]),
                     Tm_assumptions='Biopython DNA_NN3 defaults: 25 nM each strand, 50 mM Na, no Mg correction; comparative estimate only'))
assert all(len(x['reference_exact_forward_hits']) == 1 and not x['reference_exact_reverse_hits'] for x in arms)
assert backbone.endswith(left[:left_overlap_len]) and backbone.startswith(right[5:])
save_json(OUT / 'homology_arms.json', arms)
fasta(OUT / 'homology_arms.fasta', [('common_5prime_arm_125_144', left), ('common_3prime_arm_727_756', right)])

source = ROOT / 'reports/phase5/motif_attribution/ranked_library_with_parent_attribution.tsv.gz'
with gzip.open(source, 'rt') as f:
    candidates = {int(r['candidate_rank']): r for r in csv.DictReader(f, delimiter='\t') if r['sequence_role'] == 'parent' and int(r['candidate_rank']) in [140, 25, 779]}
assert set(candidates) == {140, 25, 779}
manifest, qcs = [], []


def annotated_record(sequence, name, replacement=False):
    record = SeqRecord(Seq(sequence), id=name[:16], name=name[:16], description=name + '; COMPUTATIONAL REVIEW DESIGN, NOT LAB VALIDATED')
    record.annotations = dict(molecule_type='DNA', topology='circular')
    for f in a.features:
        if f.type == 'source':
            continue
        new = copy.deepcopy(f)
        if int(f.location.start) >= 726 and replacement:
            new.location = f.location + (-82)
        assert str(f.extract(a.seq)).upper() == str(new.extract(record.seq)).upper()
        record.features.append(new)
    for f in features:
        if f['confidence'].startswith('source_') or f['feature'] in ['SYFP2 CDS']:
            continue
        start = f['pilot_500bp_start'] if replacement else f['reference_start']
        end = f['pilot_500bp_end'] if replacement else f['reference_end']
        label = 'MG pilot insert' if replacement and f['feature'] == 'AiE0779m_3xC2' else f['feature']
        record.features.append(SeqFeature(FeatureLocation(start-1, end, strand=1), type='misc_feature',
            qualifiers={'label': [label], 'note': [f['confidence'] + '; ' + f['evidence'] if label != 'MG pilot insert' else 'Natural 500-bp parent; original enhancer replaced exactly; provisional selection']}))
    return record


SeqIO.write(annotated_record(s, 'CN4461_evidence'), OUT / 'CN4461_evidence_annotated.gbk', 'genbank')
for rank in [140, 25, 779]:
    r = candidates[rank]
    insert = r['sequence'].upper()
    assert len(insert) == 500 and set(insert) <= set('ACGT')
    assert sha(insert) == r['sequence_sha256']
    name = 'pilot_parent_rank' + str(rank)
    fragment = left + insert + right
    expected = s[:144] + insert + s[726:]
    # Independent directed overlap assembly on the actual restriction-cut top
    # strand, followed by rotation back to original reference coordinate 1.
    assembled_at_right_arm = backbone + fragment[left_overlap_len:]
    assert assembled_at_right_arm[-right_overlap_len:] == assembled_at_right_arm[:right_overlap_len]
    assembled_at_right_arm = assembled_at_right_arm[:-right_overlap_len]
    origin = len(s) - 731
    assembled = assembled_at_right_arm[origin:] + assembled_at_right_arm[:origin]
    assert assembled == expected
    assert len(assembled) == 4738 and assembled[:144] == s[:144] and assembled[644:] == s[726:]
    assert assembled.count(insert) == 1
    assert s[144:338] not in assembled  # no intact old 194-bp core remains
    protected = [f for f in features if f['reference_start'] > 726 or f['reference_end'] <= 144]
    checks = {f['feature']: assembled[f['pilot_500bp_start']-1:f['pilot_500bp_end']] == f['sequence'] for f in protected}
    assert all(checks.values())
    assembled_record = annotated_record(assembled, name, True)
    SeqIO.write(assembled_record, PILOT / (name + '.gbk'), 'genbank')
    fasta(PILOT / (name + '_plasmid.fasta'), [(name, assembled)])
    fasta(PILOT / (name + '_fragment.fasta'), [(name + '_550bp_review_only', fragment)])
    qc = dict(id=name, parent_sequence_id=r['sequence_id'], parent_rank=rank, peak_id=r['peak_id'],
              design_status='COMPUTATIONAL_SEQUENCE_CHECKS_PASS_PILOT_SELECTION_PROVISIONAL',
              intended_orientation='+ relative to reference SYFP2',
              plasmid_length=len(assembled), fragment_length=len(fragment),
              vector_genome_including_long_annotated_ITRs=2141, alternative_short_annotation_span=2130,
              all_non_enhancer_reference_bases_preserved=True, protected_features=checks,
              intact_old_C2_repeat_present=False, parent_copies=1,
              insert_internal_MluI=hits(insert, 'ACGCGT'), insert_internal_SacI=hits(insert, 'GAGCTC'),
              assembled_MluI=hits(assembled, 'ACGCGT', True), assembled_SacI=hits(assembled, 'GAGCTC', True),
              GC_fraction=round((insert.count('G')+insert.count('C'))/500, 4),
              maximum_homopolymer=max(len(x) for x in __import__('re').findall(r'A+|C+|G+|T+', insert)),
              insert_exact_full_arm_matches={x['arm']: hits(insert, x['sequence_5to3']) + hits(insert, str(Seq(x['sequence_5to3']).reverse_complement())) for x in arms},
              longest_insert_arm_shared_run={x['arm']: longest_shared_run(x['sequence_5to3'], insert) for x in arms},
              plasmid_sha256=sha(assembled), fragment_sha256=sha(fragment))
    qcs.append(qc)
    manifest.append(dict(id=name, status='REVIEW_ONLY_NOT_AN_ORDER', parent_sequence_id=r['sequence_id'],
                         candidate_rank=rank, peak_id=r['peak_id'], nearest_gene=r['nearest_gene'],
                         human_genomic_coordinates=f"{r['chrom']}:{r['start']}-{r['end']} (0-based half-open)",
                         enhancer_length=500, fragment_length=550, left_arm=left, enhancer=insert,
                         right_arm=right, synthesis_fragment_5to3=fragment, fragment_sha256=sha(fragment)))
write_tsv(PILOT / 'pilot_fragment_manifest_REVIEW_ONLY.tsv', manifest)
save_json(PILOT / 'simulation_QC.json', qcs)
save_json(OUT / 'verification_summary.json', dict(
    date='2026-11-10', status='EXACT_ENHANCER_REPLACEMENT_VERIFIED_WITH_ANNOTATION_CAVEATS',
    reference_length=4820, related_reference_length=4808, reference_sequence_sha256=sha(s),
    related_sequence_sha256=sha(t), enhancer_1based_inclusive=[145,726], published_enhancer_exact_match=True,
    repeat_length=194, repeat_copies=3, nominal_minBG_length_from_publication=53,
    barcode_candidate='AGCAATTA', barcode_candidate_interval=[1832,1839], alignment_gap_interval=[1828,1839],
    full_reference_protected_interval_identity=True, candidate_source=str(source.relative_to(ROOT)),
    candidate_source_sha256=sha(source.read_bytes()), pilot_parent_ranks=[140,25,779],
    MluI=dict(recognition=[139,144], top_cut_after=139, bottom_cut_after=143, overhang='5prime CGCG'),
    SacI=dict(recognition=[727,732], top_cut_after=731, bottom_cut_after=727, overhang='3prime AGCT'),
    retained_top_strand_length=len(backbone), removed_top_strand_interval=[140,731],
    removed_top_strand_length=592, final_enhancer_replacement_length=582,
    limitations=['Exact minimal-promoter identity boundary inferred from published 53-bp length and conserved junctions',
                 'WPRE3 core exact 246 bp; terminal boundary convention not independently annotated',
                 'BarcodeS 8-bp identity inferred using published readout pattern; request depositor barcode-name mapping',
                 'Overlapping 130/141-bp right-ITR annotations unresolved; both sequences preserved',
                 'Pilot candidates illustrative, not finalized; no unique barcodes assigned',
                 'Digital overlap reconstruction does not validate physical assembly, ITR folding, packaging, or enhancer activity'] ))
save_json(OUT / 'run_parameters.json', dict(biopython_version=biopython_version,
    coordinate_convention='1-based inclusive unless explicitly labeled otherwise',
    alignment='difflib.SequenceMatcher(autojunk=False), full length plus orientation; original origins already identical',
    pilot_selection='Illustrative parents from previous activity-oriented top list, original candidate ranks 140,25,779; not final selection',
    simulation='Resolve 15-bp left and 25-bp right directed overlaps of actual restriction-cut top strand and synthetic dsDNA top strand, close circle, rotate to reference origin',
    assembly_method='NEBuilder HiFi-type two-fragment homology assembly; backbone digested separately; synthetic insert NOT restriction digested',
    excluded_claims=['No thermodynamic folding/reaction kinetics model', 'No laboratory or packaging validation', 'No barcode redesign', 'No order']))
(REF / 'all_reference_checksums.sha256').write_text(''.join(sha(p.read_bytes())+'  '+p.name+'\n' for p in sorted(REF.iterdir()) if p.is_file() and p.name != 'all_reference_checksums.sha256'))

# A readable vector-cassette map; inferred annotations explicitly styled dashed.
def draw_map(path, modified=False):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    fig, ax = plt.subplots(figsize=(13, 4.3))
    ax.plot([1, 2223-(82 if modified else 0)], [0,0], color='#34495e', lw=2)
    items = [('Left ITR',1,130,0.38,False), ('MG parent 500 bp' if modified else 'AiE0779m 3xC2',145,644 if modified else 726,0.38,False),
             ('minBG*',733,785,-0.70,True),('SYFP2 →',826,1545,0.38,False),('WPRE3 core*',1580,1825,0.38,True),
             ('BarcodeS*',1832,1839,-0.70,True),('BGHpA',1854,1965,0.38,False),('Right ITR†',2083,2223,0.38,False)]
    for n, st, en, y, inferred in items:
        if modified and st>726:st-=82;en-=82
        ax.add_patch(Rectangle((st,-0.13),en-st+1,0.26,facecolor='#f4a340' if inferred else '#43978d',edgecolor='#253746',ls='--' if inferred else '-'))
        ax.annotate(f'{n}\n{st}–{en}',xy=((st+en)/2,0),xytext=((st+en)/2,y),ha='center',va='bottom' if y>0 else 'top',fontsize=8,arrowprops=dict(arrowstyle='-',color='#667788'))
    ax.set_title(('Simulated pilot cassette' if modified else 'CN4461 reference cassette')+' | 1-based inclusive coordinates',loc='left',fontsize=14,pad=22)
    ax.set_xlim(-20,2270);ax.set_ylim(-1.25,1.0);ax.set_yticks([]);ax.set_xlabel('Plasmid reference coordinate (bp); bacterial backbone omitted from this view')
    for side in ['top','right','left']:ax.spines[side].set_visible(False)
    fig.text(.07,.035,'* Inferred or partially resolved feature boundaries; see feature evidence table.  † Two overlapping right-ITR annotations; only one right ITR.\nAll sequence outside reference bases 145–726 is preserved exactly. This is a computational map, not laboratory validation.',fontsize=9)
    fig.tight_layout(rect=[0,.17,1,1]);fig.savefig(path);fig.savefig(path.with_suffix('.png'),dpi=160);plt.close(fig)


draw_map(OUT/'CN4461_reference_cassette.svg')
draw_map(OUT/'CN4461_pilot_cassette.svg', True)
files = sorted(p for p in OUT.rglob('*') if p.is_file() and p.name!='checksums.sha256')
(OUT/'checksums.sha256').write_text(''.join(sha(p.read_bytes())+'  '+str(p.relative_to(OUT))+'\n' for p in files))
print(json.dumps({'published_enhancer_match': True, 'pilot_simulations': len(qcs), 'plasmid_bp': 4738, 'fragment_bp': 550, 'arms_5to3': [left,right], 'out': str(OUT)}, indent=2))
