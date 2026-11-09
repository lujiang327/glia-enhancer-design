#!/usr/bin/env python3
"""Create a local IGV session for a three-track, shared-scale demonstration."""
import argparse
from pathlib import Path
import xml.etree.ElementTree as ET
import json
import numpy as np
import pyBigWig


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--tracks-dir',type=Path,required=True)
    args=parser.parse_args()
    root=args.tracks_dir.resolve()
    chrom,start,end='chr13',27504892,27505392
    locus='chr13:27503893-27506392'
    files=[('MG pooled','MG_pooled.IPM_per_base.100bp.bw','0,0,0'),
           ('Astrocyte pooled','Astrocyte_pooled.IPM_per_base.100bp.bw','216,129,50'),
           ('Microglia pooled','Microglia_pooled.IPM_per_base.100bp.bw','40,120,165'),
           ('Ranked natural parent regions','MG_top1000_parents.hg38.bed','80,80,80')]
    session=ET.Element('Session',genome='hg38',locus=locus,version='3')
    resources=ET.SubElement(session,'Resources')
    panel=ET.SubElement(session,'Panel',name='DataPanel',height='420',width='1000')
    summaries=[]
    for name,filename,color in files:
        path=root/filename
        if not path.is_file():raise FileNotFoundError(path)
        ET.SubElement(resources,'Resource',name=name,path=str(path))
        track=ET.SubElement(panel,'Track',id=str(path),name=name,color=color,height='100' if filename.endswith('.bw') else '45',visible='true',autoScale='false')
        if filename.endswith('.bw'):
            ET.SubElement(track,'DataRange',baseline='0.0',drawBaseline='true',minimum='0.0',maximum='0.06',type='LINEAR')
            with pyBigWig.open(str(path)) as bw:
                values=np.array(bw.values(chrom,start,end))
                whole=np.array(bw.values(chrom,start-1000,end+1000))
                assert np.isfinite(values).all() and whole.max()<.06
                summaries.append(dict(track=name,parent_mean_density=float(values.mean()),window_max_density=float(whole.max())))
    bed_lines=(root/files[-1][1]).read_text().splitlines()
    assert any(line.split('\t')[:3]==[chrom,str(start),str(end)] for line in bed_lines)
    ET.indent(session)
    output=root/'IGV_demo_rank1.xml'
    ET.ElementTree(session).write(output,encoding='utf-8',xml_declaration=True)
    parsed=ET.parse(output)
    assert len(parsed.findall('./Resources/Resource'))==4
    assert len(parsed.findall('./Panel/Track/DataRange'))==3
    notes={'session_file':str(output),'initial_locus':locus,'candidate_interval_1based':'chr13:27504893-27505392',
        'candidate':'rank0001__retina_peak_044687','fixed_range':[0,.06],
        'signal':'observed endpoint insertion density; IPM/base; 100bp mean','observed_values':summaries,
        'validation':'XML and local data paths checked; GUI load not verified because Java window could not be controlled'}
    (root/'IGV_demo_rank1_notes.json').write_text(json.dumps(notes,indent=2)+'\n')
    print(json.dumps(notes,indent=2))

if __name__=='__main__':main()
