#!/usr/bin/env python3
"""Fail closed: discovery and design require documented model validation."""
import argparse
import hashlib
import json
from pathlib import Path

CHECKS = {
    'chromosome_split_integrity', 'training_stability_across_seeds',
    'held_out_count_performance', 'held_out_profile_performance',
    'biological_donor_generalization', 'bias_model_no_tf_leakage',
    'tf_model_no_residual_tn5_bias', 'motif_recovery_and_attribution',
}


def require_validated(path):
    record=json.loads(path.read_text())
    if record.get('status')!='VALIDATED' or any(record.get('checks',{}).get(k)!='pass' for k in CHECKS):
        raise ValueError('STOP: Müller glia ChromBPNet is not validated; phases 3–5 remain disabled.')
    checkpoint=Path(record.get('checkpoint') or '')
    if not checkpoint.is_file():raise ValueError('Missing model checkpoint')
    sha=hashlib.sha256()
    with checkpoint.open('rb') as f:
        for b in iter(lambda:f.read(2**20),b''):sha.update(b)
    if sha.hexdigest()!=record.get('checkpoint_sha256'):raise ValueError('Checkpoint checksum mismatch')
    evidence=record.get('evidence',[])
    if not evidence or any(not Path(p).is_file() for p in evidence):raise ValueError('Missing validation evidence')
    return record


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--record',type=Path,default=Path('config/model_validation.json'))
    args=p.parse_args()
    try: require_validated(args.record)
    except (ValueError,OSError,KeyError) as e:p.exit(2,str(e)+'\n')
    print('Model validation gate passed.')
