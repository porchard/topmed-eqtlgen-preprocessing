#!/usr/bin/env python
# coding: utf-8


import sys
import pandas as pd
import numpy as np
import variants
from Bio import SeqIO
import sys

EQTLGEN = sys.argv[1]
HG38_FASTA = sys.argv[2]
HG19_FASTA = sys.argv[3]
CHAIN = sys.argv[4]

# EQTLGEN = '/net/topmed10/working/porchard/rnaseq/data/eqtlgen/2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz'
# HG38_FASTA = '/net/topmed10/working/porchard/rnaseq/data/testfasta/Homo_sapiens_assembly38_noALT_noHLA_noDecoy_ERCC.fasta'
# HG19_FASTA = '/net/topmed10/working/porchard/rnaseq/data/fasta-hg19/hg19.fa'
# CHAIN = '/net/topmed10/working/porchard/rnaseq/data/chain/hg19ToHg38.over.chain.gz'

PICARD_JAR = '/sw/picard/picard.jar'


eqtlgen = pd.read_csv(EQTLGEN, sep='\t')
to_lift = eqtlgen[['SNPChr', 'SNPPos', 'SNP', 'AssessedAllele', 'OtherAllele']].drop_duplicates()
to_lift.columns = ['chrom', 'pos', 'id', 'ref', 'alt']
to_lift.chrom = 'chr' + to_lift.chrom.astype(str)


hg19 = dict()

with open(HG19_FASTA, 'r') as f:
    for record in SeqIO.parse(f, 'fasta'):
        sys.stderr.write('Loaded sequence {}\n'.format(record.id))
        hg19[record.id] = record

sys.stderr.write('Finished loading sequences\n')

to_lift['hg19_ref'] = [hg19[chrom][pos-1] for chrom, pos in zip(to_lift.chrom, to_lift.pos)]
to_lift['hg19_ref'] = to_lift['hg19_ref'].str.upper()
assert(all((to_lift.ref == to_lift.hg19_ref) | (to_lift.alt == to_lift.hg19_ref)))
assert(to_lift.id.value_counts().max() == 1)
to_lift['hg19_alt'] = [alt if ref == hg19_ref else ref for ref, alt, hg19_ref in zip(to_lift.ref, to_lift.alt, to_lift.hg19_ref)]
to_lift = to_lift[['chrom', 'pos', 'id', 'hg19_ref', 'hg19_alt']].rename(columns={'hg19_ref': 'ref', 'hg19_alt': 'alt'})


lifted, rejected = variants.lift_variants(to_lift, HG38_FASTA, CHAIN, PICARD_JAR)

lifted['topmed_id'] = lifted.chrom + '_' + lifted.pos.astype(str) + '_' + lifted.ref + '_' + lifted.alt
old_id_to_new_id = dict(zip(lifted.id, lifted.topmed_id))
topmed_id_to_chrom = dict(zip(lifted.topmed_id, lifted.chrom))
topmed_id_to_pos = dict(zip(lifted.topmed_id, lifted.pos.astype(int)))

# update
eqtlgen_hg38 = eqtlgen.drop(columns=['SNPChr', 'SNPPos', 'GeneChr', 'GenePos', 'OtherAllele'])
eqtlgen_hg38 = eqtlgen_hg38[~eqtlgen_hg38.SNP.isin(rejected.id)]
eqtlgen_hg38.SNP = eqtlgen_hg38.SNP.map(old_id_to_new_id)
eqtlgen_hg38['Zscore'] = eqtlgen_hg38.Zscore * np.where(eqtlgen_hg38.SNP.str.split('_', expand=True)[3] == eqtlgen_hg38.AssessedAllele, 1, -1)
eqtlgen_hg38 = eqtlgen_hg38.drop(columns=['AssessedAllele'])

COLUMNS = eqtlgen_hg38.columns.to_list()
eqtlgen_hg38['#chrom'] = eqtlgen_hg38.SNP.map(topmed_id_to_chrom)
eqtlgen_hg38['start'] = (eqtlgen_hg38.SNP.map(topmed_id_to_pos) - 1).astype(str)
eqtlgen_hg38['end'] = eqtlgen_hg38.SNP.map(topmed_id_to_pos).astype(str)
eqtlgen_hg38 = eqtlgen_hg38[['#chrom', 'start', 'end'] + COLUMNS]
eqtlgen_hg38.to_csv(sys.stdout, sep='\t', index=False)