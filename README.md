# Projet_reprohackaton
## Outils
 - Dépôt GitHub partagé entre tous les membres du groupe ;
 - Le gestionnaire de workflow Snakemake ;
 - Le container Docker ;
 - Le pipeline d'analyses a été exécuté sur une machine virtuelle de l'institut Pasteur.


## Contexte
Nous avons à notre disposition deux articles. Le premier a été publié en 2013  par Harbour et al. (Nat. Genet. 2013)(https://pubmed.ncbi.nlm.nih.gov/23313955/). Les chercheurs ont séquencé les ARN de patients atteints de mélanome uvéal, possédant le gène SF3B1 muté ou non. Bien que SF3B1 soit un facteur d'épissage, il n'a pas été mis en évidence de différence d'épissage entre les deux groupes de patients.
Furney et al. (Cancer Discov. 2013) (https://pubmed.ncbi.nlm.nih.gov/23313955/) ont analysé les mêmes jeux de données et ont observé un épissage différentiel entre les deux groupes de patients.

## Objectifs
Nous cherchons à reproduire une partie des analyses de données décrites dans les deux articles et nous essayons de trouver si les gènes sont exprimés de manière différentielle entre les deux groupes de patients.
Un objectif supplémentaire qui en découle est que notre pipeline d'analyses doit être reproductible.

## Données
Les données sur lesquelles nous travaillons sont disponibles au format .sra sur https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP017413&o=acc_s%3Aa.

## Méthodes
Nous avons établi ce workflow, que nous avons implémenté en utilisant le gestionnaire de workflow Snakemake :
1. Récupérer les données de séquençage au format FastQ ;
2. Récupérer le génome humain de référence et indexer le génome avec l'outil STAR ;
3. Aligner les données de séquençage sur le génome, avec les annotations du génome préalablement obtenues ;
4. Compter les reads sur le fichier BAM obtenu avec l'outil SUBREAD ;
5. Faire l'analyse statistique avec le package DESeq de R.

L'utilisation de Snakemake permet une garantie de la reproductibilité de notre workflow.
