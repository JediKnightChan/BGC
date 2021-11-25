import json
import re
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

print("Starting to process")

with open("./output/NC_015434.json", "r") as f:
    output = json.load(f)


def get_feature_subsequence(feature, sequence):
    pattern = r"\[(?P<start>\d+):(?P<end>\d+)](\((?P<chain>[+|-])\))*"
    m = re.match(pattern, feature["location"])
    if not m:
        raise ValueError("You passed wrong location format")
    start = int(m.group("start"))
    end = int(m.group("end"))

    need_reverse_complementary = m.group("chain") == "-"

    subseq = Seq(sequence[int(start):int(end)])
    if need_reverse_complementary:
        subseq = subseq.reverse_complement()
    return subseq


features = output["records"][0]["features"]
genome = output["records"][0]["seq"]["data"]

types_to_jsons = {}

for feature in features:
    type = feature["type"]
    if type not in types_to_jsons:
        types_to_jsons[type] = []
    types_to_jsons[type].append(feature)

print("All types of features:", list(types_to_jsons.keys()))

needed_types = ["CDS"]
for type, features in types_to_jsons.items():
    if type not in needed_types:
        continue
    seq_records = []
    for feature in features:
        qualifiers = feature["qualifiers"]
        seq_id = qualifiers["domain_id"][0] if "domain_id" in qualifiers else "unknown"

        # dna_sequence = get_feature_subsequence(feature, genome)
        # seq_record = SeqRecord(dna_sequence, id=seq_id)

        if "gene_kind" not in qualifiers:
            continue

        print(qualifiers["gene_kind"], qualifiers)

        seq = Seq(qualifiers["translation"][0])
        aa_sequence = SeqRecord(seq, id=seq_id)
        seq_records.append(aa_sequence)

    # with open(f"./fasta_proteins/{type}.fasta", "w") as output_handle:
    #     SeqIO.write(seq_records, output_handle, "fasta")

