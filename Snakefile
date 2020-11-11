sra_id_list=["SRR628582","SRR628583","SRR628584","SRR628585","SRR628586","SRR628587","SRR628588","SRR628589"]
list_chr=["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","MT"]

### Chargement des donnees
# Donnees de sequencage (fichiers .sra)
#for k in range(len(sra_id_list)):
#        shell("wget -O {SRAID}.sra https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-5/{SRAID}/{SRAID}.1".format(SRAID=sra_id_list[k]))

# Chromosomes
#for i in range(len(list_chr)):
#        shell("wget -O {chr}.fa.gz ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.{chr}.fa.gz".format(chr=list_chr[i]))

# Annotations du genome
#shell("wget ftp://ftp.ensembl.org/pub/release-101/gtf/homo_sapiens/Homo_sapiens.GRCh38.101.chr.gtf.gz")
#shell("gunzip Homo_sapiens.GRCh38.101.chr.gtf.gz")


### Rules
rule all:
	input: # Met tout pour le moment, on mettra que les vrais inputs quand le workflow sera complet
		#expand("{SRAID}_1.fastq.gz",SRAID=sra_id_list),expand("{SRAID}_2.fastq.gz",SRAID=sra_id_list),"ref/ref.fa",
		#"chrLength.txt", "chrName.txt", "chrNameLength.txt","chrStart.txt","genomeParameters.txt","Genome","SA","SAindex",
		#expand("{SRAID}.bam",SRAID=sra_id_list),expand("{SRAID}.bam.bai",SRAID=sra_id_list),
		"gene_output.counts","exon_output.counts"


		
rule download_sra: #téléchargement des fichiers .sra
  output:
    expand("{SRAID}.sra",SRAID=sra_id_list)
    
  run:
    for k in range(len(sra_id_list)):
      shell("wget -O {SRAID}.sra https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-5/{SRAID}/{SRAID}.1".format(SRAID=sra_id_list[k]))



rule download_chr: #téléchargement des chromosomes non dézippés 
  output:
    expand("chromosomes/{CHR}.fa.gz",CHR=list_chr)
  
  run:
    for i in range(len(list_chr)):
      shell("wget -O {chr}.fa.gz ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.{chr}.fa.gz".format(chr=list_chr[i]))



rule convert_sra_fastq: # Cree deux fichiers fastq.gz pour chaque fichier .sra (utilisation du container sratoolkit)
	input:
		"{SRAID}.sra"
	output:
		"{SRAID}_1.fastq.gz","{SRAID}_2.fastq.gz"
	singularity:
		"docker://evolbioinfo/sratoolkit:v2.10.8"
	shell:
		"fastq-dump --gzip --split-files {input}"


rule unzip_genome: # Decompresser le genome et le mettre dans un repertoire ref
	input:
		expand("chromosomes/{CHR}.fa.gz",CHR=list_chr)
	output:
		"ref/ref.fa"
	shell:
		"gunzip -c {input} > {output}"

		
rule indexation_genome: # Indexe le genome et cree de nombreux fichiers de sortie (utilisation du container STAR)
	input:
		"ref/ref.fa"	
	output:
		"ref/chrLength.txt","ref/chrName.txt","ref/chrNameLength.txt","ref/chrStart.txt","ref/genomeParameters.txt","ref/Genome","ref/SA","ref/SAindex"
	threads: 16
	singularity: 
		"docker://evolbioinfo/star:v2.7.6a"
	shell:
		"STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir ref/ --genomeFastaFiles {input}"

		
rule mapping_FastQ_files: # Aligne les sequences d'interet sur le genome (utilisation du container STAR) --> cree des fichiers BAM
	input:
		fastq1="{SRAID}_1.fastq.gz",fastq2="{SRAID}_2.fastq.gz",chr_index="ref/chrLength.txt"
	output:
		"{SRAID}.bam"
	threads: 16
	singularity:
		"docker://evolbioinfo/star:v2.7.6a"
	shell:
		"""
		STAR --outSAMstrandField intronMotif --outFilterMismatchNmax 4 --outFilterMultimapNmax 10 --genomeDir ref --readFilesIn <(gunzip -c {input.fastq1}) <(gunzip -c {input.fastq2}) --runThreadN {threads} --outSAMunmapped None --outSAMtype BAM SortedByCoordinate --outStd BAM_SortedByCoordinate --genomeLoad NoSharedMemory --limitBAMsortRAM 7000000000 > {output}
		"""

rule index_bam_files: # Indexe les fichiers BAM crees par la règle mapping_FastQ_files
	input:
		"{SRAID}.bam"
	output:
		"{SRAID}.bam.bai"
	singularity:
		"docker://evolbioinfo/samtools:v1.11"
	shell:
		"samtools index {input} {output}"


rule gene_count:
        input:
                annot="Homo_sapiens.GRCh38.101.chr.gtf",bam_files=expand("{SRAID}.bam",SRAID=sra_id_list) #--> trouver comment l'ecrire pour eviter de mettre *.bam dans la commande shell
        output:
                "gene_output.counts"
        threads: 16
        singularity:
                "docker://evolbioinfo/subread:v2.0.1"
        shell:
                "featureCounts -T {threads} -t gene -g gene_id -s 0 -a {input.annot} -o {output} {input.bam_files}"


rule exon_count:
        input:
                annot="Homo_sapiens.GRCh38.101.chr.gtf",bam_files=expand("{SRAID}.bam",SRAID=sra_id_list) #--> trouver comment l'ecrire pour eviter de mettre *.bam dans la commande shell
        output:
                "exon_output.counts"
        threads:16
        singularity:
                "docker://evolbioinfo/subread:v2.0.1"
        shell:
                "featureCounts -T {threads} -t exon -g exon_id -s 0 -a {input.annot} -o {output} {input.bam_files}"
