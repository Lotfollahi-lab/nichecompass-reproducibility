import pandas as pd
import os
import boto3

s3 = boto3.resource("s3")
artifacts_bucket = s3.Bucket("lotfollahi-wandb-artifacts")


def load_gene_orthologs():
    gene_orthologs = pd.read_csv("data/human_mouse_gene_orthologs.csv")
    return gene_orthologs

gene_orthologs = load_gene_orthologs()
os.makedirs("mouse-gene-orthologs", exist_ok=True)
file_name = "human_mouse_gene_orthologs.csv"
local_file = os.path.join("mouse-gene-orthologs", file_name)
gene_orthologs.to_csv(local_file)
