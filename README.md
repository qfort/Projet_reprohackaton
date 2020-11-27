# Projet_reprohackaton
## Outils
 - Dépôt GitHub partagé entre tous les membres du groupe ;
 - Le gestionnaire de workflow Snakemake ;
 - Le container Docker ;
 - Le pipeline d'analyses a été exécuté sur une machine virtuelle de 16 CPUs, 64 Go de RAM et 920 Go d'espace de stockage de l'Institut Pasteur.


## Contexte
Nous avons à notre disposition deux articles. Le premier a été publié en 2013  par [Harbour et al.](https://pubmed.ncbi.nlm.nih.gov/23313955/) (Nat. Genet. 2013). Les chercheurs ont séquencé les ARN de patients atteints de mélanome uvéal, possédant le gène SF3B1 muté ou non. Bien que SF3B1 soit un facteur d'épissage, il n'a pas été mis en évidence de différence d'épissage entre les deux groupes de patients.
[Furney et al.](https://pubmed.ncbi.nlm.nih.gov/23861464/) (Cancer Discov. 2013) ont analysé les mêmes jeux de données et ont observé un épissage différentiel entre les deux groupes de patients.

## Objectifs
Nous cherchons à reproduire les résultats d’analyses de RNA-seq obtenus dans les deux précédentes études sur le mélanome uvéal, par l’élaboration d’un workflow contenant différentes étapes de récupération, de transformation et d’analyse de données. Nous étudierons les différents résultats d’analyse obtenus en sortie de ce workflow et nous les comparerons à ceux des publications, fournis en données supplémentaires. Ce workflow a donc pour but de valider ou non l’obtention de données reproductibles et dans le cas échéant, nous essaierons de déterminer les éventuelles sources de variabilité. 

Un objectif supplémentaire qui en découle est que notre pipeline d'analyses doit être reproductible.

## Données
Les données de séquençage sur lesquelles nous travaillons sont disponibles au format .sra sur le site du [NCBI](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP017413&o=acc_s%3Aa).

## Méthodes
Nous avons établi ce workflow, que nous avons implémenté en utilisant le gestionnaire de workflow Snakemake :
1. Récupérer les données de séquençage au format FastQ grâce à l'outil [SRA toolkit](https://github.com/ncbi/sra-tools) ;
2. Récupérer la version GRCh38 du génome humain sur le site du [NCBI](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.39) et indexer le génome avec l'outil [STAR](https://github.com/alexdobin/STAR) ;
3. Aligner les séquences d'intérêt sur le génome de référence, avec les annotations du génome préalablement obtenues, grâce à l'outil [STAR](https://github.com/alexdobin/STAR) ;
4. Compter le nombre reads alignés sur chaque gène du génome humain et sur chaque exon, avec l'outil SUBREAD ;
5. Analyse statistique avec le package [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) de R.

L'utilisation de Snakemake permet une garantie de la reproductibilité de notre workflow.

## Exécution du workflow
### Activation de l'environnement conda
`conda activate`

`conda install singularity=3.6.3`

### Se mettre dans un dossier avec au moins 80 Go d'espace disponible
`cd ../../mnt/mydatalocal`

### Cloner le dossier GitHub
`mkdir git_repository`

`cd git_repository`

`git clone https://github.com/qfort/Projet_reprohackaton`

`cd Projet_reprohackaton`

### Exécution du fichier sraConfig pour paramétrer sratoolkit v2.10.8
`snakemake --use-singularity -s sraConfig --cores 1`

Faire exit

### Lancement du workflow
`snakemake --use-singularity --cores 16 -s Snakefile`

### Visualisation des résultats
Pour la visualisation des données d'exons, lancer le script Python 'request_exon_Ensembl.py' pour récupérer des exons de la base de données Ensembl. Ces données, stockées dans le fichier texte 'exon_table.csv', sont utilisées dans le script R 'visualisationDataExons.R'.
