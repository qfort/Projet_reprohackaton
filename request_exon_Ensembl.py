#### Librairies importations ####
import pandas as pd
import requests

#### Genes list ####
# List of genes, with differenctial splicing due to the SF3B1 mutation,
# according to the literature.

list_gene = ["ABCC5", "UQCC1", "CRNDE", "GUSBP11", "ANKHD1", "ADAM12", "F8", "GAS8"]


#### Ensembl REST-API request ####
server = "http://rest.ensembl.org"
exon_df_total = pd.DataFrame()

# For each gene and each of its transcript, we keep exons id.
for g_name in list_gene:

    ext = "/lookup/symbol/homo_sapiens/" + g_name + "?expand=1"

    r = requests.get(server + ext, headers={"Content-Type": "application/json"})

    if not r.ok:
        r.raise_for_status()
        sys.exit()

    decoded = r.json()
    liste_exon = []
    for k in range(len(decoded["Transcript"])):
        exon = decoded["Transcript"][k]
        exon = exon["Exon"]
        for i in range(len(exon)):
            for k_exon, v_exon in exon[i].items():
                if k_exon == "id":
                    liste_exon.append(v_exon)
    exon_df = pd.DataFrame(liste_exon)
    exon_df_total = pd.concat([exon_df_total, exon_df], axis=1)
    exon_df_total = exon_df_total.rename(columns={0: g_name})

exon_df_total.to_csv("exon_table.csv")