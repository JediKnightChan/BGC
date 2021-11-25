from ete3 import NCBITaxa
import pandas as pd


ncbi = NCBITaxa()


def get_tax_id_for_name(name):
    tax_ids = ncbi.get_name_translator([name])
    if name not in tax_ids:
        return
    tax_id = tax_ids[name][0]
    return tax_id


def get_tax_id_for_simple_name(name):
    name = " ".join(name.split(" ")[:2])
    return get_tax_id_for_name(name)


def get_kingdom_of_species(name):
    tax_id = get_tax_id_for_name(name)
    if tax_id is None:
        tax_id = get_tax_id_for_simple_name(name)
        if tax_id is None:
            raise ValueError(name)

    lineage = ncbi.get_lineage(tax_id)
    names = ncbi.get_taxid_translator(lineage)
    ranks = {ncbi.get_rank([taxid])[taxid]: names[taxid] for taxid in lineage}

    if "kingdom" in ranks:
        return ranks["kingdom"]
    elif "superkingdom" in ranks:
        if ranks["superkingdom"] == "Bacteria":
            return ranks["superkingdom"]
    else:
        raise ValueError


def get_kingdom_for_rows(row):
    species = row[0]
    kingdom_casual = get_kingdom_of_species(species)
    if kingdom_casual == "Bacteria":
        antismash_taxon = kingdom_casual.lower()
    elif kingdom_casual == "Fungi":
        antismash_taxon = kingdom_casual.lower()
    else:
        print(f"! Got kingdom {kingdom_casual} for species {species}, not recognized")
        return
    print(antismash_taxon)
    return antismash_taxon


genome_table = pd.read_csv("./production/ncbi_genomes.csv", header=None)
genome_table["taxon"] = genome_table.apply(get_kingdom_for_rows, axis=1)

genome_table.to_csv("./production/ncbi_genomes_taxons.csv", header=None, index=None)
