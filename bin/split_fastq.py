#!/usr/bin/env python

import argparse
import sys
import os
from collections import Counter
from itertools import chain

import pandas as pd
import pysam
from xopen import xopen


import utils
import parse_protocol

logger = utils.get_logger(__name__)

def splitInf_to_dict(file,sample):
    """
    input:  raw_sample  well            sub_sample
            sampleX	    1-9,10,11	    sampleA
            sampleX	    56,64,85,21,12  sampleB
    output:{'SampleA': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11], 'sampleB': [56, 64, 85, 21, 12]}
    """
    df = pd.read_csv(file,sep="\t",header=0)
    if '_'.join(df.columns) != "raw_sample_well_sub_sample":
        sys.exit(f'Wrong file,header should be:raw_sample\\twell\\tsub_sample')
    if sample not in df["raw_sample"].values:
        sys.exit(f"{sample} doesn't in split csv")
    df = df[df["raw_sample"]==sample]
    if df["sub_sample"].nunique() != df.shape[0]:
        sys.exit("Pleas merge the same sub_sample row")
        
    split_dict = {}
    for index,row in df.iterrows():
        well_list = []
        sub_sample = row['sub_sample']
        wells = row['well'].split(",")
        for i in wells:
            if "-" in i:
                temp = i.split("-")
                if len(temp)>2:
                    sys.exit("More than two number between '-'")
                else:
                    for x in range(int(temp[0]),int(temp[1])+1):
                        well_list.append(int(x))
            else:
                well_list.append(int(i))
        split_dict[sub_sample] = well_list
    # Whether duplicate well exists
    wells_list = list(chain.from_iterable(split_dict.values()))
    for count in Counter(wells_list).values():
        if count > 1:
            sys.exit("Duplicate well exists")
    return split_dict

def get_all_bc(file,well_dict,well_split=False):
    """
    get raw barcode and mismatch barcode
    """
    barcodes = utils.read_one_col(file)
    all_bc ={}
    for i in well_dict.keys(): 
        all_bc[i]= {"map":{},"sample":{},"well":{}}
        for j in well_dict[i]:
            all_bc[i]["map"][barcodes[j-1]] = f'well{j}'
            all_bc[i]["sample"].update(parse_protocol.get_mismatch_dict([barcodes[j-1]], 1))
            if well_split:
                all_bc[i]["well"][f'well{j}']=parse_protocol.get_mismatch_dict([barcodes[j-1]], 1)
    return all_bc

class Split_Fastq:
    def __init__(self, args):
        self.args = args
        self.fq1_list = args.fq1.split(",")
        self.fq2_list = args.fq2.split(",")
        self.fq1_number = len(self.fq1_list)
        fq2_number = len(self.fq2_list)
        if self.fq1_number != fq2_number:
            sys.exit('fastq1 and fastq2 do not have same file number!')
        
        # pattern and whitlist
        if args.protocol == 'customized':
            pattern = args.pattern
            self.pattern_dict = parse_protocol.parse_pattern(pattern)
            self.whitelist_str = args.whitelist
        else:
            protocol = args.protocol
            protocol_dict = parse_protocol.get_protocol_dict(args.assets_dir)
            if protocol == 'AccuraCode-V1':
                if args.well == 96:
                    self.whitelist_str = protocol_dict[protocol].get("well96", [])
                else:
                    self.whitelist_str = protocol_dict[protocol].get("well384", [])
            self.pattern_dict = protocol_dict[protocol]["pattern_dict"]
        
        if len(self.pattern_dict["C"])!=1:
            sys.exit("Wrong pattern,only accept one barcode position!")
        if not self.whitelist_str:
            sys.exit("Whitelist not found!")
        if " " in self.whitelist_str:
            sys.exit("Only accept one whitelist")

    def run(self):
        raw_sample = self.args.sample
        # subsample wells
        split_dict = splitInf_to_dict(self.args.split_inf, raw_sample)

        # output file define
        out_dict = {}
        for i in split_dict.keys():
            os.makedirs(f'{raw_sample}/{i}', exist_ok=True)
            out_dict[i] = {"sample":{"out_R1":f'{raw_sample}/{i}/{i}_R1.fastq',"out_R2":f'{raw_sample}/{i}/{i}_R2.fastq'}}
            if self.args.split_to_well:
                out_dict[i]["well"] = {}
                os.makedirs(f'{raw_sample}/{i}/well', exist_ok=True)
                for j in split_dict[i]:
                    well_name = "well"+str(j)
                    out_dict[i]["well"][well_name] = {"out_R1":f'{raw_sample}/{i}/well/{well_name}_R1.fastq',"out_R2":f'{raw_sample}/{i}/well/{well_name}_R2.fastq'}
        
        # open output file
        fh_fq1 = {}
        fh_fq2 = {}
        for i in out_dict.keys():
            fh_fq1[i] = {"sample":{}}
            fh_fq1[i]["sample"] = xopen(out_dict[i]["sample"]["out_R1"],"w")
            fh_fq2[i] = {"sample":{}}
            fh_fq2[i]["sample"] = xopen(out_dict[i]["sample"]["out_R2"],"w")
            if self.args.split_to_well:
                fh_fq1[i]["well"]={}
                fh_fq2[i]["well"]={}
                for j in out_dict[i]['well'].keys():
                    fh_fq1[i]["well"][j]= xopen(out_dict[i]["well"][j]["out_R1"],"w")
                    fh_fq2[i]["well"][j]= xopen(out_dict[i]["well"][j]["out_R2"],"w")
        
        # fastq
        sub_bc = get_all_bc(self.whitelist_str,split_dict,self.args.split_to_well)
        for i in range(self.fq1_number):
            with pysam.FastxFile(self.fq1_list[i], persist=False) as fq1,pysam.FastxFile(self.fq2_list[i], persist=False) as fq2:
                for entry1, entry2 in zip(fq1, fq2):
                    header1, seq1, qual1 = entry1.name, entry1.sequence, entry1.quality
                    header2, seq2, qual2 = entry2.name, entry2.sequence, entry2.quality
                    temp_bc = seq1[self.pattern_dict["C"][0].start:self.pattern_dict["C"][0].stop]
                    for j in sub_bc.keys():
                        if temp_bc in sub_bc[j]["sample"].keys():
                            fh_fq1[j]["sample"].write(f'@{header1}\n{seq1}\n+\n{qual1}\n')
                            fh_fq2[j]["sample"].write(f'@{header2}\n{seq2}\n+\n{qual2}\n')
                            if self.args.split_to_well:
                                seq_bc = sub_bc[j]["sample"][temp_bc]
                                well_num = sub_bc[j]["map"][seq_bc]
                                fh_fq1[j]["well"][well_num].write(f'@{header1}\n{seq1}\n+\n{qual1}\n')
                                fh_fq2[j]["well"][well_num].write(f'@{header2}\n{seq2}\n+\n{qual2}\n')
        # close files
        for i in out_dict.keys():
            fh_fq1[i]["sample"].close()
            fh_fq2[i]["sample"].close()
            if self.args.split_to_well:
                for j in fh_fq1[i]["well"]:
                    fh_fq1[i]["well"][j].close()
                    fh_fq2[i]["well"][j].close()
        out_json = raw_sample + ".bulk_rna.well_bc.json"
        utils.write_json(sub_bc,out_json)
        out_json = raw_sample + ".bulk_rna.fastq_inf.json"
        utils.write_json(out_dict,out_json)
        
        logger.info(out_dict)
        logger.info("Analysis finish!")

if __name__ == "__main__":
    """
    Split fastq based on information provided by the user.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--sample', required=True)
    parser.add_argument('--fq1', required=True)
    parser.add_argument('--fq2', required=True)
    parser.add_argument('--split_inf', required=True)
    parser.add_argument('--split_to_well',action='store_true')
    parser.add_argument("--assets_dir", required=True)
    parser.add_argument("--protocol", required=True,default='AccuraCode-V1')
    parser.add_argument('--well', required=True,type=int,default=384)
    parser.add_argument("--pattern")
    parser.add_argument("--whitelist")
    parser.add_argument('--version', action='version', version='1.0')
    args = parser.parse_args()
    
    runner = Split_Fastq(args)
    runner.run()