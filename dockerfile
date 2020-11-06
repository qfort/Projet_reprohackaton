# base image : Ubuntu
FROM ubuntu:16.04
RUN apt-get update --fix-missing \
&& apt-get install -y wget gzip gcc make libbz2-dev zlib1g zlib1g-dev liblzma5 liblzma-dev libncurses5 libncurses-dev bzip2 \
&& cd /usr/local \ 
# Installation de Star
&& cd /usr/local \
&& wget https://github.com/alexdobin/STAR/archive/2.7.6a.tar.gz \
&& tar -xzvf 2.7.6a.tar.gz \
&& cp STAR-2.7.6a/bin/Linux_x86_64_static/* /usr/local/bin \
&& rm -rf 2.7.6a.tar.gz STAR-2.7.6a \
#Installation de FeatureCount
&& wget -O featureCount.tar.gz https://sourceforge.net/projects/subread/files/subread-2.0.1/subread-2.0.1-source.tar.gz/download \
&& tar -zxvf featureCount.tar.gz \
&& rm -rf featureCount.tar.gz \
&& cd featureCount \
&& ./configure \
&& make \ 
&& make install \ 
&& cd /usr/local \ 
#Istallation de Samtools
&& wget -O samtools.tar.bz2 https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 \
&& tar -xjvf samtools.tar.bz2 \
&& rm -rf samtools.tar.bz2 \
&& cd samtools-1.9 \
&& ./configure \
&& make \
&& make install \
&& cd /usr/local \
#Chargement des donnees
&& for SRAID in "SRR628582" "SRR628583"; do wget -O ${SRAID}.sra https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-5/${SRAID}/${SRAID}.1;done \
&& for chr in "1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "MT"; do wget -O ${chr}.fa.gz ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.${chr}.fa.gz; done \
