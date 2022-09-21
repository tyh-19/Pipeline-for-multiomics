#!/usr/bin/env python
import pandas as pd
import argparse
import os
from tqdm import tqdm
parser = argparse.ArgumentParser(description='Summarize Salmon Quantification Result')
parser.add_argument('--indir', '-i',help="Salmon output dir",required=True)
parser.add_argument('--output','-o',help="Output quantification matrix",required=True)
parser.add_argument('--quant','-q',help="Type of quantification used in salmon",default="TPM",choices=["EffectiveLength","TPM","NumReads"])
parser.add_argument('--transcript','-t',help="Transcript information",default="genome/tx-info.txt")
args = parser.parse_args()

def main(args):
    print("Load transcript info ...")
    txDict = {}
    with open(args.transcript) as txf:
        _ = next(txf)
        for line in txf:
            line = line.strip()
            txDict[line.split("\t")[0]] = "|".join(line.split("\t"))
    txInfo = pd.Series(txDict)
    print("Done .")
    print("Load salmon's {}".format(args.quant))
    records = []
    usedField = dict(zip(["EffectiveLength","TPM","NumReads"],[2,3,4]))[args.quant]
    for sample_id in tqdm(os.listdir(args.indir)):
        path = os.path.join(args.indir,sample_id,"quant.sf")
        with open(path) as f:
            _ = next(f)
            for line in f:
                fields = line.strip().split("\t")
                txId = fields[0].split("|")[0]
                records.append((txId, sample_id, float(fields[usedField])))
    print("Done .")
    print("Write to output ...")
    table = pd.DataFrame.from_records(records)
    table.columns = ["transcript_id","sample_id",args.quant]
    matrix = table.pivot(index="transcript_id",columns="sample_id",values=args.quant)
    matrix.index = txInfo.loc[matrix.index]
    matrix = matrix[matrix.mean(axis=1)>0]
    matrix.to_csv(args.output,sep="\t")
    print("All done.")
if __name__ == "__main__":
    main(args)

