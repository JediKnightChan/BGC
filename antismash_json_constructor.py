import os
import shutil
import shutil


antismash_results_dir = "./production/antismash_output/"
for node in os.listdir(antismash_results_dir):
    full_node_path = os.path.join(antismash_results_dir, node)
    if not os.path.isdir(full_node_path):
        continue
    json_result_file = os.path.join(full_node_path, f"{node}.json")
    if not os.path.exists(json_result_file):
        print("!", json_result_file, "does not exist")
        continue
    new_json_location = os.path.join(f"./production/jsons_output/{node}.json")
    shutil.copy(json_result_file, new_json_location)