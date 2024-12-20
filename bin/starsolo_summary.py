#!/usr/bin/env python

import argparse
import os
import random
import json
import statistics
from collections import defaultdict

import pandas as pd
import pysam

import utils

class StarsoloSummary:
    def __init__(self, args):
        self.args = args

        # output file
        self.well_raw = args.sample + '.bulk_rna.counts.txt'
        self.well_filter = args.sample + '.bulk_rna.counts_report.txt'
        # multiqc report data
        multiqc_data = "./multiqc_data/"
        os.makedirs(f'./multiqc_data/', exist_ok=True)
        self.read_json = f"{multiqc_data}{args.sample}.bulk_rna.read.stats.json"
        self.summary_json = f"{multiqc_data}{args.sample}.bulk_rna.starsolo.stats.json"
        self.well_json = f"{multiqc_data}{args.sample}.bulk_rna.counts_report.json"
        self.median_gene_json = f"{multiqc_data}{args.sample}.bulk_rna.median_gene.json"

    def run(self):
        df_well, data_dict = self.parse_read_stats(self.args.read_stats)
        df_well.to_csv(self.well_raw, sep='\t')
        utils.write_json(data_dict, self.read_json)

        df_filter = df_well[ (df_well['UMI']>=args.umi_cutoff) & (df_well['read']>=args.read_cutoff) & (df_well['gene']>=args.gene_cutoff) ]
        if df_filter.shape[0]==0:
            df_filter = df_well[ (df_well['UMI']>0) & (df_well['read']>0) & (df_well['gene']>0) ]
        df_filter.to_csv(self.well_filter, sep='\t')
        utils.write_json(df_filter.to_dict('index'), self.well_json)

        summary_dict = self.parse_summary(df_filter)
        utils.write_json(summary_dict, self.summary_json)

        fraction_dcit = self.sub_gene(df_filter)
        with open(self.median_gene_json, "w") as f:
            f.write(json.dumps(fraction_dcit))

    def parse_read_stats(self,read_stats):
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
    
    def parse_summary(self,df):
        data_summary = utils.csv2dict(self.args.summary)
        stats = df.describe()
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
        return data_dict
    
    def sub_gene(self,df):
        a, cb_dict = self.get_records(self.args.bam)
        barcodes = df.index.to_list()
        barcodes = set(cb_dict[x] for x in barcodes)
        random.seed(0)
        random.shuffle(a)

        """get median gene for each fraction"""
        nread_fraction = {}
        n = len(a)
        for fraction in range(11):
            fraction /= 10.0
            nread = int(n * fraction)
            nread_fraction[nread] = fraction
        
        fraction_mg = {0.0: 0}
        cb_gx = defaultdict(set)
        for i, (cb, _, gx) in enumerate(a, start=1):
            if cb in barcodes:
                cb_gx[cb].add(gx)
            if i in nread_fraction:
                fraction = nread_fraction[i]
                fraction_mg[fraction] = int(statistics.median([len(x) for x in cb_gx.values()]))
        return fraction_mg
    
    def get_records(self,bam_file):
        a = []
        cb_int = {}
        ub_int = {}
        gx_int = {}
        n_cb = n_ub = n_gx = n_read = 0
        dup_align_read_names = set()
        with pysam.AlignmentFile(bam_file) as bam:
            for record in bam:
                cb = record.get_tag("CB")
                ub = record.get_tag("UB")
                gx = record.get_tag("GX")
                if all(x != "-" for x in (cb, ub, gx)):
                    if record.get_tag("NH") > 1:
                        if record.query_name in dup_align_read_names:
                            continue
                        else:
                            dup_align_read_names.add(record.query_name)
                    
                    # use int instead of str to avoid memory hog
                    if cb not in cb_int:
                        n_cb += 1
                        cb_int[cb] = n_cb
                    if ub not in ub_int:
                        n_ub += 1
                        ub_int[ub] = n_ub
                    if gx not in gx_int:
                        n_gx += 1
                        gx_int[gx] = n_gx
                    
                    a.append((cb_int[cb], ub_int[ub], gx_int[gx]))
                    n_read += 1
        return a, cb_int

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Starsolo summary")
    parser.add_argument("--read_stats", help="cellReadsStats file")
    parser.add_argument("--summary", help="summary file")
    parser.add_argument("--sample", help="sample name")
    parser.add_argument('--bam', required=True)
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

    runner = StarsoloSummary(args)
    runner.run()