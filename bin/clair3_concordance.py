#!/usr/bin/env python3
"""
Compare LoFreq variants against Clair3 calls and produce a concordance TSV.
Only variants with AF >= params.min_report_af in LoFreq are considered.
Clair3 is used as corroboration at AF >= 25% only.
"""

import argparse
import gzip
import sys


def parse_vcf(vcf_path, min_af=0.0):
    """Parse a bgzipped VCF and return {(chrom, pos, ref, alt): info_dict}."""
    variants = {}
    opener = gzip.open if vcf_path.endswith('.gz') else open
    with opener(vcf_path, 'rt') as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 8:
                continue
            chrom, pos, _, ref, alt = parts[0], int(parts[1]), parts[2], parts[3], parts[4]
            info = {}
            for field in parts[7].split(';'):
                if '=' in field:
                    k, v = field.split('=', 1)
                    info[k] = v
                else:
                    info[field] = True
            af = float(info.get('AF', 0))
            if af >= min_af:
                variants[(chrom, pos, ref, alt)] = {'AF': af, 'DP': info.get('DP', '.'), 'info': info}
    return variants


def main():
    parser = argparse.ArgumentParser(description='LoFreq vs Clair3 concordance')
    parser.add_argument('--lofreq', required=True, help='LoFreq filtered VCF.gz')
    parser.add_argument('--clair3', required=True, help='Clair3 merge_output.vcf.gz')
    parser.add_argument('--output', required=True, help='Output concordance TSV')
    parser.add_argument('--min-af', type=float, default=0.01, help='Minimum AF for LoFreq variants')
    parser.add_argument('--clair3-min-af', type=float, default=0.25,
                        help='Minimum AF for Clair3 corroboration (default: 0.25)')
    parser.add_argument('--sample-id', default='sample', help='Sample ID for output')
    parser.add_argument('--genotype', default='', help='Genotype label for output')
    args = parser.parse_args()

    lofreq_vars = parse_vcf(args.lofreq, min_af=args.min_af)
    clair3_vars = parse_vcf(args.clair3, min_af=args.clair3_min_af)

    with open(args.output, 'w') as out:
        out.write('\t'.join([
            'sample_id', 'genotype', 'chrom', 'pos', 'ref', 'alt',
            'lofreq_af', 'lofreq_dp', 'clair3_af', 'clair3_dp', 'concordant'
        ]) + '\n')

        all_keys = set(lofreq_vars) | set(clair3_vars)
        for key in sorted(all_keys):
            chrom, pos, ref, alt = key
            lf = lofreq_vars.get(key)
            c3 = clair3_vars.get(key)
            if lf is None:
                continue  # only report LoFreq variants; Clair3 unique calls are not reported
            concordant = 'yes' if c3 is not None else 'no'
            out.write('\t'.join([
                args.sample_id,
                args.genotype,
                chrom,
                str(pos),
                ref,
                alt,
                str(lf['AF']),
                str(lf['DP']),
                str(c3['AF']) if c3 else '.',
                str(c3['DP']) if c3 else '.',
                concordant,
            ]) + '\n')

    # Summary to stderr
    total_lf = len(lofreq_vars)
    concordant = sum(1 for k in lofreq_vars if k in clair3_vars)
    print(f"LoFreq variants: {total_lf}", file=sys.stderr)
    print(f"Concordant with Clair3 (AF>={args.clair3_min_af}): {concordant} ({100*concordant/total_lf:.1f}%)" if total_lf else "No LoFreq variants.", file=sys.stderr)


if __name__ == '__main__':
    main()
