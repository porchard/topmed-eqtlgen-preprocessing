#!/usr/bin/env python

import tempfile
import subprocess
import gzip
import pandas as pd
import re
import os
import logging
import shutil

logging.basicConfig(level=logging.WARNING, format='%(asctime)s - %(levelname)s: %(message)s')


def _write_fake_vcf(variant_df, fh):
    """
    variant_df should have columns: chrom, pos, id, ref, alt
    """
    fh.write("##fileformat=VCFv4.2\n")
    fh.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    fh.write('\t'.join(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'fakesample']) + '\n')

    for chrom, pos, id, ref, alt in zip(variant_df.chrom, variant_df.pos, variant_df.id, variant_df.ref, variant_df.alt):
        fh.write(f'{chrom}\t{pos}\t{id}\t{ref}\t{alt}\t.\tPASS\t.\tGT\t./.\n')
    return True


def _make_sequence_dict(picard_jar, fasta):
    bn, suffix = re.match('(.*)\.(.*)', fasta).groups()
    dict_name = f'{bn}.dict'
    if not os.path.exists(dict_name):
        cmd = ['java', '-Xmx4g', '-Xms4g', '-jar', picard_jar, 'CreateSequenceDictionary', f'R={fasta}', f'O={dict_name}']
        logging.info('Generating sequence dictionary (command: {})'.format(' '.join(cmd)))
        try:
            sp = subprocess.run(cmd, capture_output=True, check=True)
        except subprocess.CalledProcessError as e:
            logging.error('Error generating sequence dictionary: {}'.format(e.stderr))
            raise e
    return True


def _lift_vcf(picard_jar, fasta, chain, input_vcf, output_vcf, reject_vcf):
    cmd = ['java', '-Xmx5g', '-Xms5g', '-jar', picard_jar, 'LiftoverVcf', '--CHAIN', chain, '--INPUT', input_vcf, '--OUTPUT', output_vcf, '--REJECT', reject_vcf, '--WARN_ON_MISSING_CONTIG', 'true', '--REFERENCE_SEQUENCE', fasta, '--RECOVER_SWAPPED_REF_ALT', 'true', '--WRITE_ORIGINAL_ALLELES', 'true', '--WRITE_ORIGINAL_POSITION', 'true']
    logging.info('Lifting VCF (command: {})'.format(' '.join(cmd)))
    try:
        sp = subprocess.run(cmd, capture_output=True, check=True)
    except subprocess.CalledProcessError as e:
        logging.error('Error lifting VCF: {}'.format(e.stderr))
        raise e
    return True


def _vcf_to_df(vcf):
    output = []
    with gzip.open(vcf, 'rt') as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            output.append(line.rstrip().split('\t'))
    if len(output) == 0:
        return None
    else:
        output = pd.DataFrame(output, columns=['chrom', 'pos', 'id', 'ref', 'alt', 'qual', 'filter', 'info', 'format', 'fakesample'])
        return output


def lift_variants(variants_df, fasta, chain, picard_jar):
    """
    Input:
    variants_df: dataframe with columns chrom, pos, id, ref, alt ("id" must be unique)
    fasta: path to fasta file of target genome
    chain: path to chain file
    picard_jar: path to picard.jar

    Returns:
    (lifted_df, rejected_df)
    """
    assert(variants_df.columns.to_list() == ['chrom', 'pos', 'id', 'ref', 'alt'])

    with tempfile.TemporaryDirectory() as tmpdir:
        newfasta = os.path.join(tmpdir, os.path.basename(fasta))
        logging.info(f'Copying fasta file to {newfasta}')
        shutil.copy2(fasta, tmpdir)
        
        _make_sequence_dict(picard_jar, newfasta)
        INPUT_VCF = os.path.join(tmpdir, 'in.vcf.gz')
        OUTPUT_VCF = os.path.join(tmpdir, 'out.vcf.gz')
        REJECT_VCF = os.path.join(tmpdir, 'rejected.vcf.gz')
        assert(variants_df.id.value_counts().max()==1)
        with gzip.open(INPUT_VCF, 'wt') as fh:
            _write_fake_vcf(variants_df, fh)
        tmp = _lift_vcf(picard_jar, newfasta, chain, INPUT_VCF, OUTPUT_VCF, REJECT_VCF)

        output = _vcf_to_df(OUTPUT_VCF)
        output['swapped_alleles'] = output['info'].str.contains('SwappedAlleles')
        output['reverse_complemented_alleles'] = output['info'].str.contains('ReverseComplementedAlleles')
        output = output[['chrom', 'pos', 'id', 'ref', 'alt', 'swapped_alleles', 'reverse_complemented_alleles']]
        rejected = _vcf_to_df(REJECT_VCF)
        rejected['rejected_reason'] = rejected['filter']
        rejected = rejected[['chrom', 'pos', 'id', 'ref', 'alt', 'rejected_reason']]
        
    return (output, rejected)