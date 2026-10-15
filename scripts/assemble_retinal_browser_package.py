#!/usr/bin/env python3
"""Validate all completed cell-type tracks and assemble a portable review folder."""
import argparse
import json
import os
from pathlib import Path
import shutil
import tarfile
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyBigWig
from prepare_browser_tracks import digest


def make_archive(directory,archive):
    if archive.exists(): raise FileExistsError(archive)
    def neutral_header(info):
        info.uid=info.gid=0;info.uname=info.gname='';info.mtime=0
        info.mode=0o755 if info.isdir() else 0o644
        return info
    with tarfile.open(archive,'w') as tar:
        for path in sorted(directory.iterdir()):
            tar.add(path,arcname=directory.name+'/'+path.name,filter=neutral_header)
    return dict(file=archive.name,bytes=archive.stat().st_size,sha256=digest(archive))


def finish_sharing(directory,report_dir):
    shutil.copy2(report_dir/'README.md',directory/'README.md')
    shutil.copy2(report_dir/'MG_vs_glial_offtargets_preview.png',directory/'MG_vs_glial_offtargets_preview.png')
    sizes=Path('data/raw/reference/hg38.canonical.chrom.sizes')
    shutil.copy2(sizes,directory/sizes.name)
    (directory/'checksums.sha256').write_text(''.join(digest(p)+'  '+p.name+'\n' for p in sorted(directory.iterdir()) if p.name!='checksums.sha256'))
    full=make_archive(directory,directory.parent/(directory.name+'.tar'))
    core=directory.parent/(directory.name+'_glial_core')
    core.mkdir(exist_ok=False)
    for group in ['MG','Astrocyte','Microglia']:
        name=group+'_pooled.IPM_per_base.100bp.bw'
        os.link(directory/name,core/name)
    for name in ['MG_top1000_parents.hg38.bed','candidate_index.tsv','hg38.canonical.chrom.sizes','MG_vs_glial_offtargets_preview.png']:
        shutil.copy2(directory/name,core/name)
    (core/'README.md').write_text('''# Compact glial browser review

This subset contains three pooled observed ATAC BigWigs (MG, astrocyte, microglia), all 1,000 natural candidate regions and a six-locus preview. Use hg38 in IGV and Group Autoscale for a shared y-axis. Signal is insertion density per million insertions per base, averaged over 100bp.

This is a first-look subset: neuronal and donor tracks are in the full 25-track package. It cannot rule out neuronal off-target activity. No track represents a synthetic variant, model prediction, or reporter experiment. Donor-blocked differential-accessibility statistics remain separate inferential evidence.
''')
    (core/'checksums.sha256').write_text(''.join(digest(p)+'  '+p.name+'\n' for p in sorted(core.iterdir()) if p.name!='checksums.sha256'))
    compact=make_archive(core,core.parent/(core.name+'.tar'))
    for archive in [full,compact]:
        (directory.parent/(archive['file']+'.sha256')).write_text(archive['sha256']+'  '+archive['file']+'\n')
    result=dict(full=full,glial_core=compact)
    (report_dir/'sharing_archives.json').write_text(json.dumps(result,indent=2)+'\n')
    return result


def verify_group(root,group,config,config_sha,sizes):
    directory=root/group
    for line in (directory/'output_checksums.sha256').read_text().splitlines():
        expected,name=line.split()
        path=directory/name
        if path.resolve().parent!=directory.resolve(): raise ValueError('Unexpected checksum path')
        if digest(path)!=expected: raise ValueError('Checksum mismatch: '+str(path))
    summary=json.loads((directory/'browser_track_summary.json').read_text())
    expected=config['groups'][group]
    if summary.get('status')!='OBSERVED_CELLTYPE_BROWSER_TRACKS_READY' or summary['group']!=group or summary['assembly']!='hg38' or summary['bin_bp']!=100:
        raise ValueError('Unexpected group/assembly/bin convention')
    if summary['config_sha256']!=config_sha: raise ValueError('Configuration hash mismatch')
    observed={record['donor']:{k:record[k] for k in ['fragments','cells','sha256']} for record in summary['donors']}
    if observed!=expected['donors']: raise ValueError('Donor provenance mismatch')
    fragments=sum(v['fragments'] for v in expected['donors'].values())
    if summary['fragments']!=fragments or summary['insertions']!=2*fragments:
        raise ValueError('Pooled count mismatch')
    if group=='MG' and summary['mg_reference_check']!='PASS_EXACT_ALL_GENOME_INTEGER_BIN_MATCH':
        raise ValueError('MG reference comparison did not pass')
    desired={group+'_pooled.IPM_per_base.100bp.bw':('pooled',2*fragments)}
    if expected['keep_donor_tracks']:
        desired.update({group+'_'+donor+'.IPM_per_base.100bp.bw':(donor,2*v['fragments']) for donor,v in expected['donors'].items()})
    if {t['file'] for t in summary['tracks']}!=set(desired): raise ValueError('Unexpected track inventory')
    records=[]
    for record in summary['tracks']:
        filename=record['file'];path=directory/filename;donor,n=desired[filename]
        if path.name!=filename or record['insertions']!=n or digest(path)!=record['sha256']:
            raise ValueError('Track identity/hash/count mismatch')
        with pyBigWig.open(str(path)) as bw:
            if bw.chroms()!=sizes or not np.isclose(bw.header()['sumData'],1e6,rtol=2e-5):
                raise ValueError('Chromosome or normalized integral mismatch')
        records.append(dict(cell_type=group,donor=donor,**record))
    return summary,records


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--tracks-root',type=Path,required=True)
    p.add_argument('--output-dir',type=Path,required=True)
    p.add_argument('--report-dir',type=Path,required=True)
    p.add_argument('--config',type=Path,default=Path('config/browser_celltype_tracks.json'))
    a=p.parse_args()
    config=json.loads(a.config.read_text());config_sha=digest(a.config)
    chrom_sizes=Path('data/raw/reference/hg38.canonical.chrom.sizes')
    if digest(chrom_sizes)!=config['chrom_sizes_sha256']: raise ValueError('Genome hash mismatch')
    sizes=dict((c,int(n)) for c,n in (line.split() for line in chrom_sizes.read_text().splitlines()))
    summaries=[];records=[]
    for group in config['group_order']:
        summary,tracks=verify_group(a.tracks_root,group,config,config_sha,sizes)
        summaries.append(summary);records.extend(tracks)
        print('PASS',group,len(tracks),'tracks',flush=True)
    if len(records)!=25: raise ValueError('Expected 25 total tracks')
    a.output_dir.mkdir(parents=True,exist_ok=False)
    a.report_dir.mkdir(parents=True,exist_ok=True)
    for record in records:
        source=a.tracks_root/record['cell_type']/record['file'];destination=a.output_dir/record['file']
        try: os.link(source,destination)
        except OSError: shutil.copy2(source,destination)
    candidates=pd.read_csv('reports/phase4/parent_scoring/parent_predictions.tsv.gz',sep='\t').sort_values('candidate_rank')
    bed=candidates[['chrom','start','end']].copy()
    bed['name']=['rank'+str(rank).zfill(4)+'__'+pid for rank,pid in zip(candidates.candidate_rank,candidates.peak_id)]
    bed['score']=0;bed['strand']='.'
    bed.sort_values(['chrom','start']).to_csv(a.output_dir/'MG_top1000_parents.hg38.bed',sep='\t',index=False,header=False)
    candidates[['candidate_rank','peak_id','chrom','start','end','genomic_class','nearest_gene','maximum_offtarget_cell_type']].to_csv(a.output_dir/'candidate_index.tsv',sep='\t',index=False)
    table=pd.DataFrame(records)
    table.to_csv(a.output_dir/'track_index.tsv',sep='\t',index=False)
    table.to_csv(a.report_dir/'track_index.tsv',sep='\t',index=False)
    (a.report_dir/'celltype_source_summaries.json').write_text(json.dumps(summaries,indent=2)+'\n')
    preview_ranks=[1,200,400,600,800,1000]
    fig,axes=plt.subplots(2,3,figsize=(13,7),sharey=True)
    handles={}
    try:
        for group in ['MG','Astrocyte','Microglia']:
            handles[group]=pyBigWig.open(str(a.output_dir/(group+'_pooled.IPM_per_base.100bp.bw')))
        for ax,rank in zip(axes.flat,preview_ranks):
            row=candidates[candidates.candidate_rank.eq(rank)].iloc[0]
            left=max(0,int((row.start-1000)//100*100));right=min(sizes[row.chrom],int(np.ceil((row.end+1000)/100)*100))
            x=np.arange(left,right)-(row.start+250)
            for (group,bw),color in zip(handles.items(),['black','#d88132','#2878a5']):
                values=np.array(bw.values(row.chrom,left,right))
                if not np.isfinite(values).all(): raise ValueError('Missing values in preview')
                ax.step(x,values,where='mid',color=color,lw=1.2,label=group)
            ax.axvspan(-250,250,color='grey',alpha=.12)
            ax.set_title(f'Rank {rank}: {row.peak_id}\n{row.chrom}:{row.start+1}-{row.end}',fontsize=9)
            ax.set_xlabel('Position relative to natural parent center (bp)',fontsize=9)
            ax.set_xlim(x.min(),x.max());ax.spines[['top','right']].set_visible(False)
        for ax in axes[:,0]: ax.set_ylabel('Observed insertion density\n(IPM/base; 100bp mean)',fontsize=9)
        h,l=axes[0,0].get_legend_handles_labels();fig.legend(h,l,ncol=3,loc='upper center',bbox_to_anchor=(.5,.945),frameon=False)
        fig.suptitle('Observed accessibility: MG versus related glial off-targets',fontsize=15,weight='bold',y=.99)
        fig.text(.04,.02,'Predefined ranks 1, 200, 400, 600, 800, 1000. All panels share the same y-scale; shading marks each natural 500bp parent.\nLower-depth off-target tracks are noisier. This is observed endpoint signal, not a mutant prediction or an enhancer assay.',fontsize=9)
        fig.tight_layout(rect=[.02,.1,.99,.86],h_pad=2)
        for ext in ['png','svg']:fig.savefig(a.report_dir/('MG_vs_glial_offtargets_preview.'+ext),dpi=180)
        plt.close(fig)
    finally:
        for bw in handles.values():bw.close()
    result=dict(status='PASS_ALL_25_BROWSER_TRACKS',source_array_job=64054545,cell_types=len(summaries),
        tracks=len(records),pooled_tracks=int(table.donor.eq('pooled').sum()),donor_tracks=int(table.donor.ne('pooled').sum()),
        bytes=int(table.bytes.sum()),candidate_regions=len(candidates),mg_reference_match='exact integer bins across genome',
        preview_ranks=preview_ranks,config_sha256=config_sha,
        off_target_sequence_models_trained=False,signal='observed native-locus insertion density; IPM/base; 100bp mean')
    for directory in [a.report_dir,a.output_dir]:
        (directory/'retinal_browser_summary.json').write_text(json.dumps(result,indent=2)+'\n')
    if (a.report_dir/'README.md').exists():
        result['sharing_archives']=finish_sharing(a.output_dir,a.report_dir)
    (a.report_dir/'output_checksums.sha256').write_text(''.join(digest(path)+'  '+path.name+'\n' for path in sorted(a.report_dir.iterdir()) if path.name!='output_checksums.sha256'))
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
