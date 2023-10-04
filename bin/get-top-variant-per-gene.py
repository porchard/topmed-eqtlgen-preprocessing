#!/usr/bin/env python

import csv
import sys
import gzip
import pandas as pd
import logging

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s: %(message)s')

full_sumstats = sys.argv[1]
best_per_gene = {} # gene --> [gene, p, snp, 
COLUMNS = ['Gene', 'GeneSymbol', 'SNPChr', 'SNPPos', 'AssessedAllele', 'OtherAllele', 'Pvalue', 'FDR']
pval_index = COLUMNS.index('Pvalue')

line_count = 0
with gzip.open(full_sumstats, 'rt') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for line in reader:
        line_count += 1
        if line_count % 1e6 == 0:
            logging.info(f'Processed {round(line_count/1e6)}M lines')
        gene = line['Gene']
        if gene not in best_per_gene or float(best_per_gene[gene][pval_index]) > float(line['Pvalue']):
            best_per_gene[gene] = [line[i] for i in COLUMNS]

pd.DataFrame(list(best_per_gene.values()), columns=COLUMNS).to_csv(sys.stdout, sep='\t', index=False)

logging.info('Done')
