#!/usr/bin/env python

import sys
import logging
from Bio import SeqIO

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s: %(message)s')

fasta = sys.argv[1]

seqs = dict()

with open(fasta, 'r') as f:
    for record in SeqIO.parse(f, 'fasta'):
        logging.info('Loaded sequence {}'.format(record.id))
        seqs[record.id] = record

sys.stderr.write('Finished loading sequences\n')

for record_id, record in sorted(seqs.items()):
    logging.info('Writing sequence {}'.format(record_id))
    SeqIO.write(record, sys.stdout, 'fasta')
