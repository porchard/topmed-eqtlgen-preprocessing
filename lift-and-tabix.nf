#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process sort_and_tabix {

    publishDir "${params.results}/tabixed"
    container 'docker.io/porchard/general:20220406125608'
    memory '200 GB'
    cpus 10

    input:
    path(eqtlgen) 
    path(hg38_fasta)
    path(hg38_fasta_dict)
    path(hg19_fasta)
    path(chain)

    output:
    path("eqtlgen.txt.gz")
    path("eqtlgen.txt.gz.tbi")

    """
    lift-eqtlgen.py $eqtlgen $hg38_fasta $hg19_fasta $chain > eqtlgen.unsorted.txt
    head -n 1 eqtlgen.unsorted.txt > eqtlgen.txt
    cat eqtlgen.unsorted.txt | grep -v "^#" | sort -k1,1 -k2n,2 --parallel=10 -S 20G >> eqtlgen.txt
    bgzip eqtlgen.txt
    tabix -p bed eqtlgen.txt.gz
    rm eqtlgen.unsorted.txt
    """
}




workflow {
    sort_and_tabix(Channel.fromPath(params.eqtlgen), Channel.fromPath(params.hg38_fasta), Channel.fromPath(params.hg38_fasta_dict), Channel.fromPath(params.hg19_fasta), Channel.fromPath(params.chain))
}