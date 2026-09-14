#!/usr/bin/env python3
"""Run one ordered local preparation stage, recording the exact command.

Use the project .venv interpreter. Each stage requires earlier stage outputs.
Stages never submit HPC jobs or run differential accessibility.
"""
import argparse
import datetime
import json
from pathlib import Path
import subprocess
import sys


def run(stage):
    config=json.loads(Path('config/preparation.json').read_text())
    scripts={'reference':'prepare_reference.py','shift':'check_fragment_shift.py',
             'tracks':'build_training_tracks.py','regions':'prepare_training_regions.py',
             'validate':'validate_training_inputs.py'}
    required={'reference':['data/raw/reference/upstream.md5sum.txt'],
              'shift':['reports/preparation/reference.json'],
              'tracks':['reports/preparation/fragment_shift.json'],
              'peaks':['reports/preparation/tracks.json'],
              'regions':['data/intermediate/training/peaks/muller_relaxed_peaks.narrowPeak'],
              'validate':['reports/preparation/regions.json']}
    for path in required[stage]:
        if not Path(path).is_file():raise FileNotFoundError(f'Complete prerequisite first: {path}')
    if stage=='peaks':
        Path('data/intermediate/training/peaks').mkdir(exist_ok=True)
        cmd=[str(Path(sys.executable).parent/'macs3'),'callpeak','-t','data/intermediate/training/muller.insertions.bed.gz',
             '-f',config['peak_format'],'-g',config['peak_effective_genome'],'--keep-dup','all','--nomodel',
             '--shift',str(config['peak_shift']),'--extsize',str(config['peak_extension']),
             '-p',str(config['peak_pvalue']),'--call-summits','-n','muller_relaxed',
             '--outdir','data/intermediate/training/peaks']
    else:cmd=[sys.executable,'scripts/'+scripts[stage]]
    Path('logs').mkdir(exist_ok=True)
    with open('logs/preparation_commands.jsonl','a') as f:f.write(json.dumps(dict(stage=stage,time=datetime.datetime.now(datetime.timezone.utc).isoformat(),argv=cmd))+'\n')
    print('Running',stage,flush=True);subprocess.run(cmd,check=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('stage',choices=['reference','shift','tracks','peaks','regions','validate']);run(parser.parse_args().stage)
