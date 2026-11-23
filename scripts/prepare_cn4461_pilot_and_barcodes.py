#!/usr/bin/env python3
"""Freeze selected cloning pilot and simulate a separate pooled barcode design.

No orders. Pooled sequences are engineering prototypes, not a library release.
"""
import copy
import csv
import hashlib
import json
import random
import re
from pathlib import Path

import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
from Bio.Restriction import EcoRV, PmlI

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT/'designs/CN4461_AAV_MPRA'
OLD = BASE/'verification'
ORDER = BASE/'pilot_order'
POOL = BASE/'pooled_barcode_design'
for p in (ORDER, POOL, POOL/'engineering_examples'):
    p.mkdir(parents=True, exist_ok=True)


def sha(s):
    return hashlib.sha256(s.encode() if isinstance(s,str) else s).hexdigest()


def dump(p, data):
    p.write_text(json.dumps(data, indent=2)+'\n')


def table(p, rows, delimiter=','):
    with p.open('w') as f:
        w=csv.DictWriter(f, fieldnames=list(rows[0]), delimiter=delimiter)
        w.writeheader(); w.writerows(rows)


def positions(s, motif):
    return [i+1 for i in range(len(s)) if (s+s[:len(motif)-1]).startswith(motif,i)]


def fasta(p, rows):
    SeqIO.write([SeqRecord(Seq(s), id=n, description='') for n,s in rows],p,'fasta')


refrec=SeqIO.read(BASE/'reference/addgene-plasmid-236600-sequence-472831.gbk','genbank')
ref=str(refrec.seq).upper()
arms=json.loads((OLD/'homology_arms.json').read_text())
left,right=[a['sequence_5to3'] for a in arms]
assert left=='GTTCCTGCGGCCGCACGCGT' and right=='GAGCTCGGGCTGGGCATAAAAGTCAGGGCA'
assert ref[124:144]==left and ref[726:756]==right
with (OLD/'pilot_review_only/pilot_fragment_manifest_REVIEW_ONLY.tsv').open() as f:
    old=list(csv.DictReader(f,delimiter='\t'))
assert [int(r['candidate_rank']) for r in old]==[140,25,779]
orders=[];pilot_audit=[]
for r in old:
    name='MG_P'+r['candidate_rank']+'_parent500'
    parent=r['enhancer']; fragment=left+parent+right
    assert fragment==r['synthesis_fragment_5to3'] and sha(fragment)==r['fragment_sha256']
    actual=SeqIO.read(OLD/'pilot_review_only'/(r['id']+'.gbk'),'genbank')
    target=ref[:144]+parent+ref[726:]
    assert str(actual.seq).upper()==target
    assert len(parent)==500 and len(fragment)==550 and len(target)==4738
    assert fragment[:20]==left and fragment[20:520]==parent and fragment[520:]==right
    orders.append(dict(fragment_id=name, format='linear double-stranded DNA',
        length_bp=550, sequence_5to3=fragment,
        vendor_note='Exact sequence; no added adapters or substitutions. '+('Internal SacI in enhancer: backbone-only digestion.' if 'GAGCTC' in parent else 'Backbone-only MluI/SacI digestion.'),
        scale='Vendor standard minimum; quote before purchase',
        nearest_gene=r['nearest_gene'], parent_sequence_id=r['parent_sequence_id'],
        sequence_sha256=sha(fragment)))
    fasta(ORDER/(name+'.fasta'),[(name,fragment)])
    actual.id=name;actual.name=name[:16]
    actual.description=name+'; selected initial cloning pilot; existing 8bp BarcodeS preserved; computational design'
    for feature in actual.features:
        if 'note' in feature.qualifiers:
            feature.qualifiers['note']=[x.replace('provisional selection','user-selected initial cloning pilot') for x in feature.qualifiers['note']]
    SeqIO.write(actual,ORDER/(name+'_expected_plasmid.gbk'),'genbank')
    pilot_audit.append(dict(fragment_id=name, selected_by_user=True,
        parent_sequence_id=r['parent_sequence_id'], parent_sequence_sha256=sha(parent),
        fragment_sequence_sha256=sha(fragment), expected_plasmid_sha256=sha(target),
        expected_plasmid_length=4738, ITR_span_long_annotation=2141, ITR_span_short_annotation=2130,
        left_junction=left+'|'+parent[:20], right_junction=parent[-20:]+'|'+right,
        parent_500bp_exact_match=True, only_change_to_CN4461='replace reference 145-726 with supplied 500bp parent',
        BarcodeS_1832_1839='AGCAATTA (unchanged; pilot coordinate 1750-1757)',
        internal_MluI=positions(parent,'ACGCGT'), internal_SacI=positions(parent,'GAGCTC'),
        GC_fraction=(parent.count('G')+parent.count('C'))/500,
        max_homopolymer=max(map(len,re.findall(r'A+|C+|G+|T+',parent))),
        feature_annotation_review='minBG exact identity endpoints; WPRE3 terminal convention; BarcodeS name mapping; overlapping right-ITR endpoints',
        laboratory_validation='not yet performed'))
table(ORDER/'pilot_order.csv',orders)
fasta(ORDER/'pilot_order.fasta',[(r['fragment_id'],r['sequence_5to3']) for r in orders])
dump(ORDER/'pilot_sequence_audit.json',pilot_audit)
dump(ORDER/'order_sheet_inputs.json',dict(date='2026-11-10',orders=orders,arms=arms,audit=pilot_audit,
    verification_source='https://www.addgene.org/236600/sequences/',
    reference_sequence_sha256=sha(ref)))

# Blueprint: 16bp predesigned tags, 10 tags per enhancer/control, 4000 capacity.
# No full enhancer catalog has been chosen or assigned in this codebook.
N=4000;seed=446120261110
rng=random.Random(seed)
codebook=[];encoded=np.zeros((N,16),dtype=np.uint8);rev=np.zeros_like(encoded)
alphabet='ACGT'; trans=str.maketrans('ACGT','TGCA')
forbidden=['ACGCGT','GAGCTC','CTCGAG','GATATC','CACGTG','GGCGCGCC','AATAAA','ATTAAA']
rejected=0
while len(codebook)<N:
    bc=''.join(rng.choice(alphabet) for _ in range(16))
    if not 7<=bc.count('G')+bc.count('C')<=9:rejected+=1;continue
    rc=bc.translate(trans)[::-1]
    if re.search(r'(A{3}|C{3}|G{3}|T{3})',bc):rejected+=1;continue
    if any(bc[i:i+2]*3==bc[i:i+6] for i in range(11)):rejected+=1;continue
    if any(x in bc or x in rc for x in forbidden):rejected+=1;continue
    # Keep exactly the pre-existing XhoI at the barcode's left constant flank;
    # reject new restriction sites or polyA-like words spanning either junction.
    junction=ref[1825:1831]+bc+ref[1839:1854]
    if any(len(positions(junction,x))>int(x=='CTCGAG') for x in forbidden[:6]):rejected+=1;continue
    if any(x in junction for x in ['AATAAA','ATTAAA']):rejected+=1;continue
    v=np.array([alphabet.index(x) for x in bc],dtype=np.uint8)
    vr=np.array([alphabet.index(x) for x in rc],dtype=np.uint8)
    if (v!=vr).sum()<5:rejected+=1;continue
    n=len(codebook)
    if n and (np.any((encoded[:n]!=v).sum(axis=1)<5) or np.any((rev[:n]!=v).sum(axis=1)<5)):
        rejected+=1;continue
    encoded[n]=v;rev[n]=vr;codebook.append(bc)

# Independently audit all pairwise distances, in memory-bounded blocks.
minimum=16;minimum_rc=16
for start in range(0,N,100):
    end=min(start+100,N)
    dist=(encoded[start:end,None,:]!=encoded[None,:,:]).sum(axis=2)
    for i in range(end-start):dist[i,start+i]=17
    minimum=min(minimum,int(dist.min()))
    dr=(encoded[start:end,None,:]!=rev[None,:,:]).sum(axis=2)
    minimum_rc=min(minimum_rc,int(dr.min()))
assert minimum>=5 and minimum_rc>=5 and len(set(codebook))==N
book=[dict(barcode_id=f'BC{i+1:04}',barcode_5to3=bc,length_bp=16,
           assignment_status='UNASSIGNED_CODEBOOK_ENGINEERING_ONLY',
           GC_fraction=(bc.count('G')+bc.count('C'))/16) for i,bc in enumerate(codebook)]
table(POOL/'barcode_codebook_4000.csv',book)
fasta(POOL/'barcode_codebook_4000.fasta',[(r['barcode_id'],r['barcode_5to3']) for r in book])
dump(POOL/'codebook_QC.json',dict(seed=seed,length=16,barcodes=N,barcodes_per_enhancer_proposed=10,
    capacity_enhancers_including_controls=400,minimum_pairwise_hamming=minimum,
    minimum_hamming_to_reverse_complements_including_self=minimum_rc,
    GC_allowed_counts=[7,8,9],max_homopolymer=2,perfect_dinucleotide_6mer_excluded=True,
    forbidden_motifs=forbidden,junctions_screened=True,assignment='No full-library assignment; first 30 used only in engineering prototypes',
    position_base_fractions=[{a:float((encoded[:,j]==i).mean()) for i,a in enumerate(alphabet)} for j in range(16)],
    candidates_rejected=rejected,decoder='Extract known flanks, exactly 16nt; unique expected whitelist exact match or <=1 substitution; reject ambiguous/indel/low-quality tags',
    limitations=['Not empirically neutral barcodes','Not a minimal-edit-distance code','Hamming screening does not correct indels','No final experimental catalog assigned']))

# Replacement strategy verified against the actual reference. XhoI is NOT
# unique; EcoRV/PmlI are. Re-synthesize the intervening cassette literally,
# changing only 1832-1839 from 8 to 16bp. All ITRs remain in retained backbone.
sites={name:positions(ref,motif) for name,motif in dict(XhoI='CTCGAG',EcoRV='GATATC',PmlI='CACGTG',AscI='GGCGCGCC').items()}
assert sites['XhoI']==[205,399,593,1826]
assert sites['EcoRV']==[1574] and sites['PmlI']==[2061]
assert EcoRV.search(Seq(ref),linear=False)==[1577] and PmlI.search(Seq(ref),linear=False)==[2064]
bl,br=ref[1549:1579],ref[2060:2080]
fasta(POOL/'barcode_cassette_arms.fasta',[('EcoRV_left_1550_1579',bl),('PmlI_right_2061_2080',br)])
expected=[];cassette_rows=[];prototype_qc=[];prototype_sequences=[]
for j,r in enumerate(old):
    parent=r['enhancer'];pseq=ref[:144]+parent+ref[726:]
    assert EcoRV.search(Seq(pseq),linear=False)==[1495] and PmlI.search(Seq(pseq),linear=False)==[1982]
    prec=SeqIO.read(ORDER/(orders[j]['fragment_id']+'_expected_plasmid.gbk'),'genbank')
    for k in range(10):
        item=book[j*10+k];bc=item['barcode_5to3'];name=orders[j]['fragment_id']+'_'+item['barcode_id']
        cassette=ref[1579:1831]+bc+ref[1839:2060]
        fragment=bl+cassette+br
        target=pseq[:1749]+bc+pseq[1757:]
        # Actual cut backbone: pilot cut after1494 and1981, resolve27/17bp overlaps.
        backbone=pseq[1981:]+pseq[:1494]
        assert backbone.endswith(bl[:27]) and backbone.startswith(br[3:])
        assembled=backbone+fragment[27:]
        assert assembled[-17:]==assembled[:17]
        assembled=assembled[:-17];origin=len(pseq)-1981
        assembled=assembled[origin:]+assembled[:origin]
        assert assembled==target and len(target)==4746 and len(fragment)==539
        assert target[:1749]==pseq[:1749] and target[1765:]==pseq[1757:]
        assert target[144:644]==parent
        assert EcoRV.search(Seq(target),linear=False)==[1495] and PmlI.search(Seq(target),linear=False)==[1990]
        # Preserve all native and evidence annotations outside old barcode.
        short_id='P'+r['candidate_rank']+'_'+item['barcode_id']
        rec=SeqRecord(Seq(target),id=short_id,name=short_id,description=name+'; ENGINEERING PROTOTYPE ONLY; not in initial pilot order')
        rec.annotations=dict(molecule_type='DNA',topology='circular')
        for feature in prec.features:
            lo,hi=int(feature.location.start),int(feature.location.end)
            if lo<1757 and hi>1749:continue # old tag or junction annotation replaced explicitly
            new=copy.deepcopy(feature)
            if lo>=1757:new.location=feature.location+8
            assert str(feature.extract(prec.seq)).upper()==str(new.extract(rec.seq)).upper()
            rec.features.append(new)
        rec.features.append(SeqFeature(FeatureLocation(1749,1765,strand=1),type='misc_feature',qualifiers={'label':[item['barcode_id']],'note':['Designed 16bp tag; requires experimental validation']}))
        SeqIO.write(rec,POOL/'engineering_examples'/(name+'.gbk'),'genbank')
        prototype_sequences.append((name,target))
        expected.append(dict(prototype_construct_id=name,enhancer_id=r['parent_sequence_id'],
            enhancer_sha256=sha(parent),barcode_id=item['barcode_id'],barcode_5to3=bc,
            barcode_replicate=k+1,status='ENGINEERING_EXAMPLE_NOT_INITIAL_ORDER',
            expected_plasmid_sha256=sha(target),expected_plasmid_bp=4746,ITR_span_long=2149,ITR_span_short=2138))
        cassette_rows.append(dict(prototype_construct_id=name,length_bp=539,
            cassette_with_arms_5to3=fragment,sha256=sha(fragment),status='ENGINEERING_ONLY_NOT_ORDER'))
        prototype_qc.append(dict(id=name,exact_directed_assembly=True,all_nonbarcode_pilot_bases_unchanged=True,
            original_parent_500bp_unchanged=True,protected_feature_sequences_unchanged=True,
            EcoRV_sites=positions(target,'GATATC'),PmlI_sites=positions(target,'CACGTG'),
            new_barcode_1based=[1750,1765],new_BGHpA_signal_1based=[1780,1891],
            right_ITR_long_1based=[2009,2149]))
table(POOL/'engineering_examples_expected_associations.tsv',expected,'\t')
fasta(POOL/'engineering_expected_plasmids.fasta',prototype_sequences)
table(POOL/'engineering_barcode_cassettes_NOT_ORDER.csv',cassette_rows)
dump(POOL/'engineering_simulation_QC.json',prototype_qc)
dump(POOL/'barcode_replacement_design.json',dict(reference_recognition_site_starts=sites,
    replacement_8bp=[1832,1839],new_barcode_length=16,
    reference_cassette_between_complete_sites=[1580,2060],reference_cassette_length=481,new_cassette_length=489,
    left_flank=dict(sequence=bl,interval=[1550,1579],length=30,actual_overlap_bp=27),
    right_flank=dict(sequence=br,interval=[2061,2080],length=20,actual_overlap_bp=17),
    EcoRV_cut_after=1576,PmlI_cut_after=2063,removed_top_strand=[1577,2063],removed_top_strand_length=487,
    engineering_fragment_bp=539,only_base_change='Replace 8bp candidate BarcodeS with new16bp tag; retain all flanks and cassette sequence',
    scope='Validated digitally on selected three parents; screen sites anew for every eventual pooled enhancer/variant',
    linearization='Backbone separately EcoRV/PmlI digested; cassette not digested; both original sites restored',
    reference_unmodified_with_new_tag_bp=4828,pilot_parent_with_new_tag_bp=4746))
print(json.dumps(dict(initial_pilot_fragments=len(orders),initial_pilot_fragment_bp=550,
    barcode_codebook=N,minimum_Hamming=minimum,minimum_Hamming_RC=minimum_rc,
    barcode_engineering_examples=len(expected),barcode_cassette_fragment_bp=539),indent=2))
