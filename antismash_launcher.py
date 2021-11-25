import os.path

import pandas as pd
import requests
import subprocess


def download_genome_from_ncbi(genome_url, filename):
    print(genome_url)
    r = requests.get(genome_url)
    if r.status_code != 200:
        print("Error", r.status_code, filename)
        return
    with open(f"./production/genomes/{filename}.fna.gz", "wb") as f:
        f.write(r.content)


def run_antismash_processes(row):
    filename = " - ".join(row[0:3])
    filename = filename.replace("/", "_").replace("\\", "_")
    taxon = row[72]
    genefinding_tool = "glimmerhmm" if taxon == "fungi" else "prodigal"

    if not os.path.exists(f"./production/genomes/{filename}.fna.gz"):
        genome_url = row[3]
        download_genome_from_ncbi(genome_url, filename)
    subprocess.run(
        f"antismash --output-dir '/home/jediknight/Documents/Sirius/BGC/production/antismash_output/{filename}/' --genefinding-tool"
        f" {genefinding_tool} --taxon {taxon} '/home/jediknight/Documents/Sirius/BGC/production/genomes/{filename}.fna.gz'",
        shell=True)


genome_table = pd.read_csv("./production/ncbi_genomes_taxons.csv", header=None)
genome_table = genome_table[genome_table[72].notnull()]

genome_table.apply(run_antismash_processes, axis=1)
