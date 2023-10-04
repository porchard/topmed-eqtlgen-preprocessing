ROOT=/net/topmed11/working/porchard/eqtlgen-preprocessing
DATA=$(ROOT)/data
WORK=$(ROOT)/work
BIN=$(ROOT)/bin

ANALYSIS=$(WORK)/$@
SIF=singularity exec --bind $(ROOT) $(ROOT)/general.sif

.PHONY: all

define NL


endef

data: eqtlgen fasta-hg19 fasta-hg38 chain

singularity:
	singularity pull general.sif docker://porchard/general:20220406125608

eqtlgen:
	mkdir -p $(DATA)/$@
	wget  https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz --directory-prefix $(DATA)/$@/
	wget https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz --directory-prefix $(DATA)/$@/
	wget https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/README_cis --directory-prefix $(DATA)/$@/
	wget https://molgenis26.gcc.rug.nl/downloads/eqtlgen/trans-eqtl/2018-09-04-trans-eQTLsFDR-CohortInfoRemoved-BonferroniAdded.txt.gz --directory-prefix $(DATA)/$@/
	wget https://molgenis26.gcc.rug.nl/downloads/eqtlgen/trans-eqtl/2018-09-04-trans-eQTLsFDR0.05-CohortInfoRemoved-BonferroniAdded.txt.gz  --directory-prefix $(DATA)/$@/

fasta-hg19:
	mkdir -p $(DATA)/fasta/hg19
	cd $(DATA)/fasta/hg19 && wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz && zcat hg19.fa.gz > hg19.unsorted.fa
	cd $(DATA)/fasta/hg19 && $(SIF) python $(BIN)/sort-fasta.py hg19.unsorted.fa > hg19.fa && rm hg19.unsorted.fa && $(SIF) samtools faidx hg19.fa

fasta-hg38:
	mkdir -p $(DATA)/fasta/hg38
	wget https://console.cloud.google.com/storage/browser/_details/gtex-resources/references/Homo_sapiens_assembly38_noALT_noHLA_noDecoy_ERCC.fasta --directory-prefix $(DATA)/fasta/hg38/
	wget https://console.cloud.google.com/storage/browser/_details/gtex-resources/references/Homo_sapiens_assembly38_noALT_noHLA_noDecoy_ERCC.dict --directory-prefix $(DATA)/fasta/hg38/

chain:
	mkdir -p $(DATA)/$@
	wget http://hgdownload.cse.ucsc.edu/goldenpath/hg38/liftOver/hg38ToHg19.over.chain.gz --directory-prefix $(DATA)/$@/

top-hit-per-gene:
	mkdir -p $(ANALYSIS)
	$(SIF) python $(BIN)/get-top-variant-per-gene.py $(DATA)/eqtlgen/2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz > $(ANALYSIS)/top-per-gene.txt

lift-and-tabix:
	mkdir -p $(ANALYSIS)
	cd $(ANALYSIS) && nohup nextflow run -resume --eqtlgen $(DATA)/eqtlgen/2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz --hg19_fasta $(DATA)/fasta/hg19/hg19.fa --hg38_fasta $(DATA)/fasta/hg38/Homo_sapiens_assembly38_noALT_noHLA_noDecoy_ERCC.fasta --hg38_fasta_dict $(DATA)/fasta/hg38/Homo_sapiens_assembly38_noALT_noHLA_noDecoy_ERCC.dict --chain $(DATA)/chain/hg38ToHg19.over.chain.gz --results $(ANALYSIS)/results $(ROOT)/lift-and-tabix.nf &
