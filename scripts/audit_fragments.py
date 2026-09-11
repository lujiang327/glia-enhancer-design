#!/usr/bin/env python3
"""Audit pooled GEO fragments without mistaking reused 10x barcodes for cells.

Counts records once; column 5 is read support, not independent molecules.
No donor assignment is made from ambiguous barcode matches.
"""
import argparse
import collections
import csv
import gzip
import hashlib
import json
from pathlib import Path
import statistics


def audit(fragments, rna_barcodes, out):
    out.mkdir(parents=True, exist_ok=True)
    rna = collections.defaultdict(set)
    with gzip.open(rna_barcodes, 'rt') as handle:
        for line in handle:
            sample, barcode = line.strip().split('_', 1)
            rna[barcode].add(sample)
    counts = collections.Counter()
    nuclear = collections.Counter()
    chromosomes = collections.Counter()
    lengths = collections.Counter()
    support = collections.Counter()
    suffixes = collections.Counter()
    canonical = {f'chr{i}' for i in range(1, 23)} | {'chrX'}
    rows = 0
    with gzip.open(fragments, 'rt') as handle:
        for line in handle:
            if line.startswith('#'):
                continue
            fields = line.rstrip().split('\t')
            if len(fields) != 5:
                raise ValueError(f'Expected 5 fields at record {rows + 1}')
            chrom, start, end, barcode, reads = fields
            start, end, reads = int(start), int(end), int(reads)
            if start < 0 or end <= start or reads < 1:
                raise ValueError(f'Invalid fragment at record {rows + 1}')
            rows += 1
            counts[barcode] += 1
            chromosomes[chrom] += 1
            lengths[min(end - start, 2000)] += 1
            support[reads] += 1
            suffixes[barcode.rsplit('-', 1)[-1]] += 1
            if chrom in canonical:
                nuclear[barcode] += 1
            if rows % 10000000 == 0:
                print(f'Audited {rows:,} records', flush=True)
    ambiguity = collections.Counter()
    for barcode, n in counts.items():
        ambiguity[str(len(rna.get(barcode, set())))] += n
    with (out / 'barcode_audit.tsv').open('w') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow(['bare_barcode', 'fragment_records', 'chr1_22_X_records', 'matching_RNA_sample_prefixes', 'n_matching_RNA_samples'])
        for barcode, n in sorted(counts.items()):
            samples = sorted(rna.get(barcode, set()))
            writer.writerow([barcode, n, nuclear[barcode], ','.join(samples), len(samples)])
    for name, counter in [('chromosome_records', chromosomes), ('fragment_lengths', lengths), ('read_support', support)]:
        with (out / f'{name}.tsv').open('w') as handle:
            writer = csv.writer(handle, delimiter='\t')
            writer.writerow([name, 'records'])
            writer.writerows(sorted(counter.items()))
    sha = hashlib.sha256()
    with open(fragments, 'rb') as handle:
        for block in iter(lambda: handle.read(2**20), b''):
            sha.update(block)
    result = dict(fragment_file=str(fragments), sha256=sha.hexdigest(),
                  bytes=fragments.stat().st_size, fragment_records=rows,
                  read_support_sum=sum(k*v for k,v in support.items()),
                  chr1_22_X_records=sum(nuclear.values()),
                  distinct_bare_barcodes=len(counts),
                  median_records_per_bare_barcode=statistics.median(counts.values()),
                  median_chr1_22_X_records_per_bare_barcode=statistics.median(nuclear.get(b, 0) for b in counts),
                  barcode_suffix_record_counts=dict(suffixes),
                  records_by_number_of_matching_RNA_samples=dict(ambiguity),
                  barcodes_by_number_of_matching_RNA_samples=dict(collections.Counter(str(len(rna.get(b,set()))) for b in counts)),
                  biological_cell_count=None, usable_fragments_after_blacklist=None,
                  warning='Bare barcodes are not globally unique cell IDs. Chromosome-filtered records are provisional, not fully QC-filtered usable depth.')
    (out / 'fragment_audit.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fragments', type=Path, required=True)
    parser.add_argument('--rna-barcodes', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    audit(args.fragments, args.rna_barcodes, args.out)
