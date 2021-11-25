#!/usr/bin/env python

import argparse
import os
import time
import sys
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO

try:
    from urllib2 import URLError
except ImportError:
    from urllib.error import URLError


E_VALUE_THRESH = 0.001
NUM_WWW_ATTEMPTS = 3
INTER_TRY_DELAY = 1  # in seconds


def info(msg):
    sys.stdout.write(str(msg) + '\n')
    sys.stdout.flush()


def warn(msg):
    sys.stderr.write('WARNING! ' + str(msg) + '\n')
    sys.stderr.flush()


class BlastHit(object):
    def __init__(self, query, alignment, hsp):
        query_len = len(query.seq)
        self.query_id = query.id
        self.query_description = query.description
        self.title = alignment.title
        self.description = alignment.hit_def
        self.accession = alignment.accession
        self.acc_length = alignment.length
        self.bit_score = hsp.bits
        self.query_cover = 100.0 * (hsp.align_length - hsp.gaps) / query_len
        self.e_value = hsp.expect
        self.per_idy = 100.0 * hsp.identities / hsp.align_length

    def __str__(self):
        '''
        standard BLAST output is:
            Description
            Scientific Name
            Max Score
            Total Score
            Query Cover
            E value
            Per. Ident
            Acc. Len
            Accession
        '''
        return '\t'.join(map(str, [self.query_description, self.description, self.accession,
                                   int(self.bit_score), self.query_cover, self.e_value, '%.2f%%' % self.per_idy]))

    @staticmethod
    def get_header():
        return '\t'.join(['Query Description',
                          'Match Description', 'Accession', 'Score (bits)', 'Query Cover', 'E-value', 'Per. Ident'])


def blast_fasta(fpath, protein=False, num_top_hits=1):
    if not os.path.isfile(fpath):
        info('File does not exists! ' + fpath)
        return


    info(BlastHit.get_header())
    for record in SeqIO.parse(fpath, "fasta"):
        attempts = 0
        result_handle = None
        while attempts < 3:
            try:
                if protein:
                    result_handle = NCBIWWW.qblast("blastp", "nr", record.seq, hitlist_size=num_top_hits)
                else:
                    result_handle = NCBIWWW.qblast("blastn", "nr", record.seq, hitlist_size=num_top_hits)
                break
            except URLError:
                attempts += 1
                time.sleep(INTER_TRY_DELAY)
        if not result_handle:
            warn('Failed to retrieve BLAST results for Record: %s\n' % str(record.id))
        else:
            blast_records = NCBIXML.parse(result_handle)
            to_print = ''
            for blast_record in blast_records:
                for alignment in blast_record.alignments:
                    for hsp in alignment.hsps:
                        if hsp.expect < E_VALUE_THRESH:
                            to_print += str(BlastHit(record, alignment, hsp)) + '\n'
                        break
                break
            info(to_print)


def main():
    parser = argparse.ArgumentParser(description='Multi-fasta BLAST utility. Keeps only best hit per entry. '
                                                 'Works with both DNA and protein sequences')
    parser.add_argument(
        'inputs',
        metavar='FILE/DIR',
        type=str,
        nargs='+',
        help='paths to fasta files with sequences to BLAST'
    )
    parser.add_argument(
        '-p', '--proteins',
        action='store_true',
        default=False,
        help='Input sequences are proteins (default is DNA seqs)'
    )
    parser.add_argument(
        '-t', '--num-top-hits',
        default=1,
        help='Number of top BLAST hits per entry to save (default 1)'
    )

    args = parser.parse_args()
    for input_path in args.inputs:
        blast_fasta(input_path, args.proteins, args.num_top_hits)


if __name__ == "__main__":
    main()