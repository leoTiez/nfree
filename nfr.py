#!/usr/bin/env python3
import sys
import os
import argparse
import warnings
import pandas as pd
import numpy as np
import gtfparse
import matplotlib.pyplot as plt
from datahandler import reader


def parse_arguments(args):
    parser = argparse.ArgumentParser(
        description='Calculate the Nucleosome Free Region (NFR). Pass an MNase-seq bigwig (.bw/.bigwig) file which'
                    'represents the nucleosome distribution and a bed file (.bed) which indicates the transcripts. '
    )

    parser.add_argument('--bw', type=str, required=True,
                        help='Path to bigwig file with the MNase-seq data.')
    parser.add_argument('--annot', type=str, required=True,
                        help='Path to annotation file with coordinates for transcripts. Accept GTF and BED.')
    parser.add_argument('--out', type=str,
                        help='Path to output directory. Default is current directory.')
    parser.add_argument('--smooth_ws', type=int, default=50,
                        help='Size of smoothing window. Default is 50.')
    parser.add_argument('--mind_nfr', type=int, default=100,
                        help='Minimum distance from the TSS to NFR. Default is 50.')
    parser.add_argument('--maxd_total', type=int, default=500,
                        help='Maximum distance from the TSS for +1 and -1 nucleosome. Default is 500.')
    parser.add_argument('--sig_nfold', type=int, default=3,
                        help='Determines how much larger another peak must be to use it as +1/-1 nucleosome instead of'
                             'initially computed value.'
                        )
    parser.add_argument('--mind_nucl', type=int, default=150,
                        help='Minimum distance between +1 and -1 nucleosome. Default is 150.')
    parser.add_argument('--p1_tss', action='store_true', dest='is_p1_tss',
                        help='When flag is set, the +1 nucleosome is always the closest peak to the tss instead of '
                             'being based on the computation of the NFR'
                        )
    parser.add_argument('--scale', type=int, default=1,
                        help='Scale smoothed signal with that factor. Default is 1.'
                        )
    parser.add_argument('--verbosity', type=int, default=0,
                        help='Verbosity flag regulates output on console. '
                             '0: no output; '
                             '1: steps in computation; '
                             '2: print gene; '
                             '3: plot gene. '
                             'Default is 0.'
                        )

    parsed_args = parser.parse_args(args)
    return parsed_args


def calc_sig_better(candidates, mnase_data, sig_nfold=3, is_positive=True):
    index = 0 if is_positive else -1
    sig_better = mnase_data[candidates][index] * sig_nfold < mnase_data[candidates]
    if np.sum(sig_better) > 0:
        return candidates[sig_better][index]
    else:
        return candidates[index]


def main(args):
    pargs = parse_arguments(args)
    bw_path = pargs.bw
    annot_path = pargs.annot
    smooth_ws = pargs.smooth_ws
    sig_nfold = pargs.sig_nfold
    mind_nfr = pargs.mind_nfr
    is_p1_tss = pargs.is_p1_tss
    maxd_total = pargs.maxd_total
    scaling = pargs.scale
    mind_nucl = pargs.mind_nucl
    out_path = pargs.out
    verbosity = pargs.verbosity

    if verbosity > 0:
        print('Load bw file')
    mnase_file = reader.load_big_file(bw_path, is_abs_path=True)

    annot_type = annot_path.lower().split('.')[-1]
    if annot_type == 'gtf':
        if verbosity > 0:
            print('Load GTF file')
        annot_file = gtfparse.read_gtf(annot_path)
        annot_file = annot_file[annot_file['feature'] == 'transcript']
    elif annot_type == 'bed':
        if verbosity > 0:
            print('Load bed file')
        annot_file = reader.load_bam_bed_file(annot_path, is_abs_path=True)
    else:
        raise ValueError('Allow only GTF or bed files as annotation files.')

    if verbosity > 0:
        print('Compute NFRs')
    nfr_l = []
    nfr_depth_l = []
    gname_l = []
    if annot_type == 'gtf':
        gene_iterator = annot_file.iterrows()
    else:
        gene_iterator = annot_file.__iter__()

    for gene in gene_iterator:
        if annot_type == 'gtf':
            gene = gene[1]
            chromosome, start, end, direct, name = (
                gene['seqname'],
                int(gene['start']),
                int(gene['end']),
                gene['strand'],
                gene['gene_name']
            )
        else:
            chromosome, start, end, direct, name = gene[0], int(gene[1]), int(gene[2]), gene[5], gene[3]

        start = end if direct == '-' else start
        if verbosity > 1:
            print('Calculate NFR for %s at %s:%s with direction %s' % (name, chromosome, start, direct))
        # Load data to array
        max_bp = int(mnase_file.chroms(chromosome))
        mnase_data_raw = np.nan_to_num(mnase_file.values(
            chromosome,
            np.maximum(start - maxd_total, 0),
            np.minimum(start + maxd_total, max_bp)
        ))

        distance = np.abs(np.maximum(start - maxd_total, 0) - np.minimum(start + maxd_total, max_bp))
        # Smooth data
        mnase_data = np.convolve(mnase_data_raw, np.ones(smooth_ws)/float(smooth_ws), mode='same') * scaling

        # Compute gradient
        grad = np.gradient(mnase_data)
        grad = np.convolve(grad, np.ones(smooth_ws)/float(smooth_ws), mode='same') * scaling

        # Compute maxima and minima
        asign = np.sign(grad)
        descending = np.arange(distance)[(np.roll(asign, 1) - asign) > 0]
        if descending.size < 2:
            warnings.warn(
                'Could not find enough peaks for %s on chromosome %s. Skip.' % (name, chromosome),
                UserWarning
            )
            continue
        if not is_p1_tss:
            ascending = np.arange(distance)[(np.roll(asign, 1) - asign) < 0]
            # Compute position of nfr
            if direct == '+':
                nfr_centre_candidate = ascending[np.logical_and(
                    maxd_total + mind_nfr > ascending,
                    ascending > 1
                )]
            else:
                nfr_centre_candidate = ascending[np.logical_and(
                    maxd_total - mind_nfr < ascending,
                    ascending < maxd_total * 2
                )]
            try:
                nfr_centre = nfr_centre_candidate[np.argmin(mnase_data[nfr_centre_candidate])]
            except ValueError:
                warnings.warn(
                    'Could not find NFR centre for %s on chromsome %s. Skip.' % (name, chromosome),
                    UserWarning
                )
                continue

            # compute p1 and m1
            p1_candidates = descending[descending > nfr_centre]
            if p1_candidates.size == 0:
                p1_candidates = np.asarray([distance - 1])
            try:
                p1 = calc_sig_better(p1_candidates, mnase_data, sig_nfold)
            except IndexError:
                warnings.warn('Could not find +1 for %s on chromsome %s. Skip.' % (name, chromosome), UserWarning)
                continue
            m1_candidates = descending[np.logical_and(descending < nfr_centre, descending < p1 - mind_nucl)]
            if m1_candidates.size == 0:
                m1_candidates = np.asarray([0])
            try:
                m1 = calc_sig_better(m1_candidates, mnase_data, sig_nfold, is_positive=False)
            except IndexError:
                warnings.warn('Could not find -1 for %s on chromosome %s. Skip.' % (name, chromosome), UserWarning)
                continue

            depth = np.abs(np.maximum(mnase_data[p1], mnase_data[m1]) - mnase_data[nfr_centre])
        else:
            try:
                p1 = descending[np.argmin(np.abs(descending - maxd_total))]
            except IndexError:
                warnings.warn('Could not find +1 for %s on chromosome %s. Skip.' % (name, chromosome), UserWarning)
                continue
            # Compute position of nfr
            if direct == '+':
                m1_candidates = descending[descending <= p1 - mind_nucl]
                try:
                    m1 = calc_sig_better(m1_candidates, mnase_data, sig_nfold, is_positive=False)
                except IndexError:
                    warnings.warn('Could not find -1 for %s on chromosome %s. Skip.' % (name, chromosome), UserWarning)
                    continue
                depth = np.abs(np.maximum(mnase_data[m1], mnase_data[p1]) - np.min(mnase_data[m1:p1]))
            else:
                m1_candidates = descending[descending >= p1 + mind_nucl]
                try:
                    m1 = calc_sig_better(m1_candidates, mnase_data, sig_nfold)
                except IndexError:
                    warnings.warn('Could not find -1 for %s on chromosome %s. Skip.' % (name, chromosome), UserWarning)
                    continue
                depth = np.abs(np.maximum(mnase_data[m1], mnase_data[p1]) - np.min(mnase_data[p1:m1]))

        nfr = np.asarray([m1, p1])
        if verbosity > 2:
            x_position = np.arange(
                np.maximum(start - maxd_total, 0),
                np.minimum(start + maxd_total, max_bp)
            )
            plt.figure(figsize=(8, 7))
            plt.plot(x_position, mnase_data_raw, color='tab:green', alpha=.5)
            plt.plot(x_position, mnase_data, color='tab:red')
            plt.vlines(start, ymin=0, ymax=np.max(mnase_data), linestyle='--', color='black', linewidth=4)
            plt.scatter(x_position[nfr], mnase_data[nfr], color='red', s=9**2)
            plt.xlabel('Pos (bp)', fontsize=24)
            plt.ylabel('MNase Amplitude', fontsize=24)
            plt.xticks([start - maxd_total, start, start + maxd_total], [-maxd_total, 'TSS', maxd_total], fontsize=20)
            plt.yticks(fontsize=20)
            plt.title('NFR %s: %s' % (name, np.abs(nfr[0] - nfr[1])), fontsize=32)
            plt.show()

        nfr_depth_l.append(depth)
        nfr_l.append(np.abs(nfr[0] - nfr[1]))
        gname_l.append(name)

    df = pd.DataFrame()
    df['name'] = gname_l
    df['nfr'] = nfr_l
    df['depth'] = nfr_depth_l
    out_path = os.getcwd() if out_path is None else out_path
    out_path = '%s/nrf.tsv' % out_path
    if verbosity > 0:
        print('Save NFRs to file: %s' % out_path)
    df.to_csv('%s' % out_path, sep='\t')
    if verbosity > 0:
        print('Finished')


if __name__ == '__main__':
    main(args=sys.argv[1:])

