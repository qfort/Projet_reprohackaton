sra_id_list = ['SRR.......']#pour le moment on a un seul id SRA mais pour après

rule all:#rule final
	input:
		expand("{SRAID}.fastq.gz", SRAID=sra_id_list)

rule download_SRA_file:
	input:
		"{SRAID}"#le SRA id
	output:
		"{SRAID}.sra"#le fichier .sra
	shell:
		"wget -O ${SRAID}.sra https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-5/${SRAID}/${SRAID}.1"#le wget pour récupérer le fichier sra l'addresse change en fonction du sra id

rule convert_sra_fastq:#conversion du sra en fastq
	input:
		"{SRAID}.sra"
	output:
		"{SRAID}.fastq.gz"
	shell:
		"fastq-dump --gzip --split-files {SRAID}.sra --outdir $pwd"#La commande en utilisant le docker, ce n'est pas la même que celle du slide (celle du slide marche sur mon pc). Après je n'ai pas testé cette variante donc peut être à tester de façon individuel avant.
#j'ai pas trop vu comment utiliser un docker dans le snakefile.

