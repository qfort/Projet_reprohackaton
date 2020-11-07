# Ajouter les containers

sra_id_list=["SRR628582","SRR628583","SRR628584","SRR628585","SRR628586","SRR628587","SRR628588","SRR628589"]
list_chr=["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","MT"]
id_bam=["1","2","3","4","5","6","7","8"] # COMBIEN DE FICHIERS BAM ON ATTEND ?

# Chargement des donnees
for k in range(len(sra_id_list)):
        shell("wget -O {SRAID}.sra https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-5/{SRAID}/{SRAID}.1".format(SRAID=sra_id_list[k]))

for i in range(len(list_chr)):
        shell("wget -O {chr}.fa.gz ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.{chr}.fa.gz".format(chr=list_chr[i]))

shell("wget ftp://ftp.ensembl.org/pub/release-101/gtf/homo_sapiens/Homo_sapiens.GRCh38.101.chr.gtf.gz")

# Rules
rule all:	# rule finale
	input: # Met tout pour le moment, on mettra que les vrais inputs quand le workflow sera complet
		expand("{SRAID}_1.fastq.gz",SRAID=sra_id_list),expand("{SRAID}_2.fastq.gz",SRAID=sra_id_list),"ref/ref.fa",
		"chrLength.txt", "chrName.txt", "chrNameLength.txt","chrStart.txt","genomeParameters.txt","Genome","SA","SAindex" 


rule convert_sra_fastq: #conversion du sra en fastq
	input:
		"{SRAID}.sra"
	output:
		"{SRAID}_1.fastq.gz","{SRAID}_2.fastq.gz"
	singularity:
		"docker://evolbioinfo/sratoolkit:v2.10.8"
	shell:
		 "vdb-config -i | fastq-dump --gzip --split-files {input}"


rule unzip_genome:#recuperer le genome et le mettre dans un repertoire ref
	input:
		expand("{CHR}.fa.gz",CHR=list_chr)
	output:
		"ref/ref.fa"
	shell:
		"gunzip -c {input} > {output}"

rule indexation_genome:
	input:
		"ref/ref.fa"	
	output:
		"chrLength.txt","chrName.txt","chrNameLength.txt","chrStart.txt","genomeParameters.txt","Genome","SA","SAindex"
	threads: 16
	singularity: 
		"docker://evolbioinfo/star:v2.7.6a"
	shell:
		"STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir ref/ --genomeFastaFiles {input}"

rule mapping_FastQ_files:
	input:
		fastq1="{SRAID}_1.fastq.gz",fastq2="{SRAID}_2.fastq.gz"
	output:
		"{id_bam}.bam"
	threads: 4
	shell:
		"STAR --outSAMstrandField intronMotif --outFilterMismatchNmax 4 --outFilterMultimapNmax 10 \
--genomeDir ref \
--readFilesIn <(gunzip -c {input.fastq1) <(gunzip -c {input.fastq2}) \ #Verifier la syntaxe
--runThreadN {threads} \
--outSAMunmapped None \
--outSAMtype BAM SortedByCoordinate \
--outStd BAM_SortedByCoordinate \
--genomeLoad NoSharedMemory \
--limitBAMsortRAM 1000000000 \
> {output}"

rule create_bam_file: # Pas sure : c’est pour case la commande samtools index *.bam
	input: "id_bam"
	output: "AUCUNE IDEE"
	shell: "samtools index {input}.bam"

