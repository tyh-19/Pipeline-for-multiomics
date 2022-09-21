#!/usr/bin/env python
import pandas as pd
import argparse
from collections import defaultdict
from tqdm import tqdm
import numpy as np

parser = argparse.ArgumentParser(description='Count promoter reactivity')
parser.add_argument('--input','-i',help="Input TPM matrix by transcripts",required=True)
parser.add_argument('--output','-o',help="Output TPM  matrix by promoters",required=True)
parser.add_argument('--promoters','-p',help="Transcript grouped by TSS",default="genome/promoter/tx2tss.10.txt")
args = parser.parse_args()

def main(args):
    #ENST00000002596.5|ENSG00000002587.9|HS3ST1|protein_coding|8031
    print("Load TPM matrix ...")
    matrix = pd.read_csv(args.input,sep="\t",index_col=0).fillna(0)
    print("Done .")
    n_genes = len(set(matrix.index.map(lambda x:x.split("|")[1])))
    gene_info = {}
    for tx in matrix.index:
        fields = tx.split("|")
        if fields[1] not in gene_info.keys():
            gene_info[fields[1]] = "|".join(fields[1:])
    gene_info = pd.Series(gene_info)
    print("Input matrix contains {} transcripts, {} corresponding genes".format(matrix.shape[0],n_genes))

    txUsed = set(matrix.index.map(lambda x:x.split("|")[0])) 
    print("Load promoter information ...")
    tx2tss = pd.read_csv(args.promoters,sep="\t")
    tx2tss = tx2tss[tx2tss["tx_name"].isin(txUsed)]
    tssGroupNumberByGene = tx2tss.groupby("gene.id").apply(lambda x:np.unique(x["tss.group"].values).shape[0])            
    #Only consider genes with more than one TSS group
    geneConsidered = set(tssGroupNumberByGene[tssGroupNumberByGene>1].index)
    print("{} genes with more than 1 TSS group are considered".format(len(geneConsidered)))
    tx2tss = tx2tss[tx2tss["gene.id"].isin(geneConsidered)]
    print("{} transcripts of these genes were considered".format(tx2tss.shape[0]))
    
    print("Aggregate expression by promoters ...")
    tx2tss = tx2tss.set_index("tx_name")
    gene_ids = gene_info.loc[tx2tss["gene.id"].values]
    gene_ids.index = tx2tss.index
    promoter_info = gene_ids + "|" +  tx2tss["tss.group"]
    matrix.index = matrix.index.map(lambda x:x.split("|")[0])
    matrix = matrix.loc[tx2tss.index,:]
    
    matrix["promoter_id"] = promoter_info.loc[matrix.index]
    matrix.to_csv("test2.txt",sep="\t")
    tpmByPromoter = matrix.groupby("promoter_id").apply(lambda x:x.sum(axis=0))
    del tpmByPromoter["promoter_id"]
    tpmByPromoter.to_csv(args.output,sep="\t")
    print("Done .")



if __name__ == "__main__":
    main(args)

