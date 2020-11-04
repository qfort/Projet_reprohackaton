# Projet_reprohackaton
## Outils
 - Dépôt GitHub partagé entre tous les membres du groupe.
 - Le gestionnaire de workflow Snakemake.
 - Le container Docker.
 - Le pipeline d'analyses a été exécuté sur une machine virtuelle de l'institut Pasteur.


## Contexte
Nous avons à notre disposition deux articles [a citer] parlant de [expliquer]

## Objectifs
Reproduire une partie des analyses de données décrites dans les deux articles et essayer de trouver des gènes exprimés de manière différentielle.
Un objectif supplémentaire qui en découle est que notre pipeline d'analyses doit être reproductible.

## Données
Les données sur lesquelles nous travaillons sont disponibles au format .sra sur https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-5/

## Méthodes
Nous avons établi ce workflow, que nous avons implémenté en utilisant le gestionnaire de workflow Snakemake :
	1. Récupérer les données de séquençage au format FastQ ;
	2. Récupérer le génome humain de référence et indexer le génome avec l'outil STAR ;
	3. Aligner les données de séquençage sur le génome, avec les annotations du génome préalablement obtenues ;
	4. Compter les reads sur le fichier BAM obtenu avec l'outil SUBREAD ;
	5. Faire l'analyse statistique avec le package DESeq de R.

L'utilisation de Snakemake permet une garantie de la reproductibilité de notre workflow.
