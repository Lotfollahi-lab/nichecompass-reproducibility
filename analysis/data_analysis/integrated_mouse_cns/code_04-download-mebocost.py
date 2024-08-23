import pandas as pd
import os


def load_mebocost():
    mebocost_files = os.listdir("data/metabolite_enzyme_sensor_gps")
    mebocost_datasets = []
    for file in mebocost_files:
        dataset = pd.read_csv(os.path.join("data/metabolite_enzyme_sensor_gps", file), sep="\t")
        mebocost_datasets.append(dataset)
    return mebocost_files, mebocost_datasets


mebocost_files, mebocost_datasets = load_mebocost()
os.makedirs("mebocost-gene-programs", exist_ok=True)
for file, dataset in zip(mebocost_files, mebocost_datasets):
    file_name = file
    local_file = os.path.join("mebocost-gene-programs", file_name)
    dataset.to_csv(local_file, sep="\t")
