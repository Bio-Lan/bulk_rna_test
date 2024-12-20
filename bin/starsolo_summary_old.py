#!/usr/bin/env python

import argparse
from collections import defaultdict

import pandas as pd

import utils
import parse_protocol

def parse_read_stats(read_stats):
    dtypes = defaultdict(lambda: "int")
    dtypes["CB"] = "object"
    df = pd.read_csv(
        read_stats, sep="\t", header=0, index_col=0, skiprows=[1], dtype=dtypes
    )
    df_bc = df.loc[:,["nUMIunique",'countedU','nGenesUnique']]
    df_bc.columns = ['UMI', 'read', 'gene']
    df_bc = df_bc.sort_values('UMI', ascending=False)

    df = df.loc[
        :, ["cbMatch", "cbPerfect", "genomeU", "genomeM", "exonic", "intronic", "exonicAS", "intronicAS", "countedU"]
    ]
    s = df.sum()
    # json does not recognize NumPy data types. TypeError: Object of type int64 is not JSON serializable
    valid = int(s["cbMatch"])
    perfect = int(s["cbPerfect"])
    corrected = valid - perfect
    genome_uniq = int(s["genomeU"])
    genome_multi = int(s["genomeM"])
    mapped = genome_uniq + genome_multi
    exonic = int(s["exonic"])
    intronic = int(s["intronic"])
    antisense = int(s["exonicAS"] + s["intronicAS"])
    intergenic = mapped - exonic - intronic - antisense
    counted_uniq = int(s["countedU"])
    data_dict = {
        "Corrected Barcodes": corrected / valid,
        "Reads Mapped To Unique Loci": genome_uniq / valid,
        "Reads Mapped To Multiple Loci": genome_multi / valid,
        "Reads Mapped Uniquely To Transcriptome": counted_uniq / valid,
        "Mapped Reads Assigned To Exonic Regions": exonic / mapped,
        "Mapped Reads Assigned To Intronic Regions": intronic / mapped,
        "Mapped Reads Assigned To Intergenic Regions": intergenic / mapped,
        "Mapped Reads Assigned Antisense To Gene": antisense / mapped,
    }
    for k in data_dict:
        data_dict[k] = utils.get_frac(data_dict[k])

    return df_bc, data_dict

def well_bctonum(df, file):
    barcodes = utils.read_one_col(file)
    n = 1
    data = {}
    for i in barcodes:
        data[i] = f"well{n}"
        n += 1
    df["BC"] = df.index
    df.index = df.index.map(data)
    return(df)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Starsolo summary")
    parser.add_argument("--read_stats", help="cellReadsStats file")
    parser.add_argument("--summary", help="summary file")
    parser.add_argument("--sample", help="sample name")
    parser.add_argument('--assets_dir', required=True)
    parser.add_argument('--protocol', required=True)
    parser.add_argument('--umi_cutoff', default=500, type=int,
        help='If the UMI number exceeds the threshold, it is considered a valid well and reported.'
    )
    parser.add_argument('--gene_cutoff', default=0, type=int,
        help='If the gene number exceeds the threshold, it is considered a valid well and reported.'
    )
    parser.add_argument('--read_cutoff', default=0, type=int,
        help='If the read number exceeds the threshold, it is considered a valid well and reported.'
    )
    args = parser.parse_args()


    df_well, data_dict = parse_read_stats(args.read_stats)

    if args.protocol == 'AccuraCode-V1':
        protocol_dict = parse_protocol.get_protocol_dict(args.assets_dir)
        whitelist_str = protocol_dict[args.protocol].get("well384", [])
        df_well = well_bctonum(df_well,whitelist_str)

    # out file
    read_stats_file = args.sample + ".bulk_rna.read.stats.json"
    summary_file = args.sample + ".bulk_rna.starsolo.stats.json"

    raw_count_file = args.sample + '.bulk_rna.counts.txt'
    marked_count_file = args.sample + '.bulk_rna.counts_report.txt'
    marked_count_json = args.sample + '.bulk_rna.counts_report.json'
    
    utils.write_json(data_dict, read_stats_file)

    # Detailed information per Well
    df_well.to_csv(raw_count_file, sep='\t')
    df_well_valid = df_well[ (df_well['UMI']>=args.umi_cutoff) & (df_well['read']>=args.read_cutoff) & (df_well['gene']>=args.gene_cutoff) ]
    if df_well_valid.shape[0]==0:
        df_well_valid = df_well[ (df_well['UMI']>0) & (df_well['read']>0) & (df_well['gene']>0) ]
    df_well_valid.to_csv(marked_count_file, sep='\t')
    utils.write_json(df_well_valid.to_dict('index'), marked_count_json)
    
    data_summary = utils.csv2dict(args.summary)
    stats = df_well_valid.describe()
    data_dict = {
        "Raw Reads" : int(data_summary['Number of Reads']),
        "Valid Reads" : utils.get_frac(data_summary["Reads With Valid Barcodes"]),
        "Median Reads per Well" : int(stats.loc['50%',"read"]),
        "Median UMI per Well" : int(stats.loc['50%',"UMI"]),
        "Median Genes per Well" : int(stats.loc['50%',"gene"]),
        "Mean Reads per Well" : int(stats.loc["mean","read"]),
        "Mean UMI per Well" : int(stats.loc["mean","UMI"]),
        "Mean Genes per Well" : int(stats.loc["mean","gene"])
    }
    # summary
    utils.write_json(data_dict, summary_file)