sra_id_list=["SRR628582","SRR628583","SRR628584","SRR628585","SRR628586","SRR628587","SRR628588","SRR628589"]
list_chr=["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","MT","X","Y"]


### Rules
rule all:
	input: 
		"DE_analysis.txt", "DE_genes_comparison.txt","exons_analysis.txt"


### PART 1 : DOWNLOAD DATA
rule download_sra: # Download .sra files containing sequencing data of 8 samples in the sra_folder folder.
  output:
    expand("sra_folder/{SRAID}.sra",SRAID=sra_id_list)
    
  run:
    for k in range(len(sra_id_list)):
      shell("wget -O sra_folder/{SRAID}.sra https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-5/{SRAID}/{SRAID}.1".format(SRAID=sra_id_list[k]))



rule download_chr: # Download unzipped chromosomes of the GRCh38 version of the human genome in the chr_folder folder.
  output:
    expand("chr_folder/{CHR}.fa.gz",CHR=list_chr)
  
  run:
    for i in range(len(list_chr)):
      shell("wget -O chr_folder/{chr}.fa.gz ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.{chr}.fa.gz".format(chr=list_chr[i]))


rule download_annotations_genome: # Download genome annotations for the GRCh38 genome version.
	output:
		"Homo_sapiens.GRCh38.101.chr.gtf"
	shell:
		"""
		wget ftp://ftp.ensembl.org/pub/release-101/gtf/homo_sapiens/Homo_sapiens.GRCh38.101.chr.gtf.gz
		gunzip Homo_sapiens.GRCh38.101.chr.gtf.gz
		"""


rule convert_sra_fastq: # Create two fastq.gz files for each .sra file (use sratoolkit container). One file contains the genetic sequence in 5'3' and the other in 3'5'.
	input:
		"sra_folder/{SRAID}.sra"
	output:
		"fastq_folder/{SRAID}_1.fastq.gz","fastq_folder/{SRAID}_2.fastq.gz"
	singularity:
		"docker://evolbioinfo/sratoolkit:v2.10.8"
	shell:
		"fastq-dump --gzip --outdir fastq_folder/ --split-files {input}"


		
### PART 2 : MAPPING SEQUENCES		
rule unzip_genome: # Unzip human genome and put it in the ref folder.
	input:
		expand("chr_folder/{CHR}.fa.gz",CHR=list_chr)
	output:
		"ref/ref.fa"
	shell:
		"gunzip -c {input} > {output}"

		
rule indexation_genome: # Index human genome and create many output files (use STAR container).
	input:
		"ref/ref.fa"	
	output:
		"ref/chrLength.txt","ref/chrName.txt","ref/chrNameLength.txt","ref/chrStart.txt","ref/genomeParameters.txt","ref/Genome","ref/SA","ref/SAindex"
	threads: 16
	singularity: 
		"docker://evolbioinfo/star:v2.7.6a"
	shell:
		"STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir ref/ --genomeFastaFiles {input}"

		
rule mapping_FastQ_files: # Map samples of interest on the human genome (use STAR container). Creates BAM files in the bam_folder folder.
	input:
		fastq1="fastq_folder/{SRAID}_1.fastq.gz",fastq2="fastq_folder/{SRAID}_2.fastq.gz",chr_index="ref/chrLength.txt"
	output:
		"bam_folder/{SRAID}.bam"
	threads: 16
	singularity:
		"docker://evolbioinfo/star:v2.7.6a"
	shell:
		"""
		STAR --outSAMstrandField intronMotif --outFilterMismatchNmax 4 --outFilterMultimapNmax 10 --genomeDir ref --readFilesIn <(gunzip -c {input.fastq1}) <(gunzip -c {input.fastq2}) --runThreadN {threads} --outSAMunmapped None --outSAMtype BAM SortedByCoordinate --outStd BAM_SortedByCoordinate --genomeLoad NoSharedMemory --limitBAMsortRAM 7000000000 > {output}
		"""

rule index_bam_files: # Index BAM files created by the rule mapping_FastQ_files.
	input:
		"bam_folder/{SRAID}.bam"
	output:
		"bam_folder/{SRAID}.bam.bai"
	singularity:
		"docker://evolbioinfo/samtools:v1.11"
	shell:
		"samtools index {input} {output}"


rule gene_count: # Counts the number of reads mapped on each gene of the human genome, to get the level of expression of the genes (use the subread container). Creates two files, especially gene_output.counts, used in the analysis rules.
        input:
                annot="Homo_sapiens.GRCh38.101.chr.gtf",bam_files=expand("bam_folder/{SRAID}.bam",SRAID=sra_id_list),bai_files=expand("bam_folder/{SRAID}.bam.bai",SRAID=sra_id_list)
		# Pour avoir 8 files : bam_files="bam_folder/{SRAID}.bam"
        output:
                "gene_output.counts"
        threads: 16
        singularity:
                "docker://evolbioinfo/subread:v2.0.1"
        shell:
                "featureCounts -p -T {threads} -t gene -g gene_id -s 0 -a {input.annot} -o {output} {input.bam_files}"


rule exon_count: # Counts the number of reads mapped on the exons, to get the level of expression of the genes in case of alternative splicing (use the subread container). Creates two files, especially exon_output.counts, used in the analysis rules.
        input:
                annot="Homo_sapiens.GRCh38.101.chr.gtf",bam_files=expand("bam_folder/{SRAID}.bam",SRAID=sra_id_list),bai_files=expand("bam_folder/{SRAID}.bam.bai",SRAID=sra_id_list)
		# Pour avoir 8 files : bam_files="bam_folder/{SRAID}.bam"
        output:
                "exon_output.counts"
        threads:16
        singularity:
                "docker://evolbioinfo/subread:v2.0.1"
        shell:
                "featureCounts -p -T {threads} -t exon -g exon_id -s 0 -a {input.annot} -o {output} {input.bam_files}"


### PART 3 : STATISTICAL ANALYSIS
rule statsAnalysis: # Statistical analysis on the samples, using the counting of reads per gene done in the gene_count rule.
	input:
		"SraRunTable.txt","AnalyseGenes.R","gene_output.counts"
	output:
		"DE_analysis.txt", "DE_genes_comparison.txt"
	container:
		"docker://evolbioinfo/deseq2:v1.28.1"
	script:
		"AnalyseGenes.R"

rule statAnalysisExon: # Statistical analysis on the samples, using the counting of reads per exon done in the exon_count rule.
	input:
		"exon_output.counts"
	output:
		"exons_analysis.txt"
	container:
		"docker://evolbioinfo/deseq2:v1.28.1"
	script:
		"analyseDataExons.R"
	
