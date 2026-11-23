#!/usr/bin/env python3
"""Screen existing library inserts without inferring unavailable CN4461 DNA.

This produces no homology arms, assembly designs, or order-ready sequences.
"""
import hashlib
import json
from pathlib import Path
import re
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'designs/CN4461_AAV_MPRA'
SOURCE=ROOT/'reports/phase5/motif_attribution/ranked_library_with_parent_attribution.tsv.gz'
SITES={'MluI':'ACGCGT','SacI':'GAGCTC'}


def digest(s):return hashlib.sha256(s.encode('ascii')).hexdigest()


def main():
    out=OUT/'preflight';out.mkdir(parents=True,exist_ok=True)
    df=pd.read_csv(SOURCE,sep='\t')
    assert df.sequence_id.is_unique and len(df)==2884
    rows=[]
    for r in df.itertuples():
        seq=r.sequence
        assert isinstance(seq,str) and seq==seq.upper()
        assert digest(seq)==r.sequence_sha256
        flags=[]
        invalid=''.join(sorted(set(seq)-set('ACGT')))
        if invalid:flags.append('non_ACGT')
        gc=(seq.count('G')+seq.count('C'))/len(seq)
        if not .25<=gc<=.75:flags.append('GC_outside_review_range_25_75_percent')
        longest=max((len(m.group()) for m in re.finditer(r'A+|C+|G+|T+',seq)),default=0)
        if longest>=8:flags.append('homopolymer_ge_8')
        # Explicitly identify maximal perfect two-base repeats with overlapping starts.
        max_di=0
        for i in range(len(seq)-3):
            motif=seq[i:i+2]
            if motif[0]==motif[1]:continue
            j=i+2
            while seq[j:j+2]==motif:j+=2
            max_di=max(max_di,j-i)
        if max_di>=12:flags.append('perfect_dinucleotide_repeat_ge_12bp')
        row=dict(sequence_id=r.sequence_id,sequence_role=r.sequence_role,peak_id=r.peak_id,
                 insert_length_bp=len(seq),sequence_sha256=r.sequence_sha256,
                 GC_fraction=gc,non_ACGT_characters=invalid,longest_homopolymer_bp=longest,
                 longest_perfect_dinucleotide_repeat_bp=max_di,
                 review_flags=';'.join(flags),assembly_conflicts='NOT_ASSESSED_REFERENCE_MISSING')
        for enzyme,motif in SITES.items():
            positions=[m.start()+1 for m in re.finditer('(?='+motif+')',seq)]
            row[enzyme+'_internal_site_count']=len(positions)
            row[enzyme+'_site_start_1based']=';'.join(map(str,positions))
        rows.append(row)
    table=pd.DataFrame(rows)
    table.to_csv(out/'insert_sequence_screen.tsv',sep='\t',index=False)
    summary=dict(status='REFERENCE_SEQUENCE_REQUIRED_NO_CLONING_DESIGN_RELEASED',
        reference=dict(plasmid='CN4461',addgene_id=236600,
                       official_url='https://www.addgene.org/236600/sequences/',
                       nucleotide_sequence_available=False,
                       reason='Official download/copy requires login; no matching local GenBank/SnapGene file found'),
        source_library=str(SOURCE.relative_to(ROOT)),
        source_library_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        insert_sequences_screened=len(df),by_role=df.sequence_role.value_counts().to_dict(),
        all_source_sequence_checksums_match=True,
        sequence_lengths={str(k):int(v) for k,v in df.sequence.str.len().value_counts().items()},
        unique_insert_sequences=int(df.sequence.nunique()),
        non_ACGT_sequences=int(table.non_ACGT_characters.ne('').sum()),
        review_flagged_sequences=int(table.review_flags.ne('').sum()),
        internal_sites={enzyme:dict(recognition_sequence=motif,sequences_with_site=int(table[enzyme+'_internal_site_count'].gt(0).sum()),
                                   total_sites=int(table[enzyme+'_internal_site_count'].sum())) for enzyme,motif in SITES.items()},
        caveats=['Sequence flags are transparent review heuristics, not vendor acceptance criteria or automatic exclusions.',
                 'Internal MluI/SacI sites do not automatically prevent homology assembly with separately linearized backbone.',
                 'Backbone sites, enhancer boundaries, junctions, overlaps, protected features, vector lengths and circular assemblies are unverified.',
                 'No homology arms, final synthesis fragments, unique MPRA barcodes or order-ready vendor files were generated.'])
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    (out/'output_checksums.sha256').write_text(''.join(hashlib.sha256(p.read_bytes()).hexdigest()+'  '+p.name+'\n' for p in sorted(out.iterdir()) if p.name!='output_checksums.sha256'))
    print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
