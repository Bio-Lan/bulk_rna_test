#!/usr/bin/env python
import argparse
import gzip
import sys
import os

import pandas as pd
import pyfastx

import utils
import parse_protocol
logger = utils.get_logger(__name__)

class Split_Fastq:
    def __init__(self, args):
        self.args = args
        self.fq1_list = args.fq1.split(",")
        self.fq2_list = args.fq2.split(",")
        self.fq1_number = len(self.fq1_list)
        fq2_number = len(self.fq2_list)
        if self.fq1_number != fq2_number:
            sys.exit('fastq1 and fastq2 do not have same file number!')
        
        if args.protocol == 'customized':
            pattern = args.pattern
            self.pattern_dict = parse_protocol.parse_pattern(pattern)
        else:
            protocol = args.protocol
            protocol_dict = parse_protocol.get_protocol_dict(args.assets_dir)
            self.pattern_dict = protocol_dict[protocol]["pattern_dict"]

        if len(self.pattern_dict["C"])!=1:
            sys.exit("Wrong pattern,only accept one barcode position!")
    
    def run(self):
        raw_sample = self.args.sample
        BC_dict = self.get_bc_dict(self.args)
        
        # output file define
        out_dict = {}
        for i in BC_dict.keys():
            os.makedirs(f'{raw_sample}/{i}', exist_ok=True)
            out_dict[i] = {
                "read_num": { i:0 },
                "subsample":{"out_R1":f'{raw_sample}/{i}/{i}_R1.fastq.gz',"out_R2":f'{raw_sample}/{i}/{i}_R2.fastq.gz'},
                "subwell": {}
            }
            if self.args.split_to_well:
                os.makedirs(f'{raw_sample}/{i}/well', exist_ok=True)
                for j in BC_dict[i]["wellBC"].keys():
                    out_dict[i]["read_num"][j] = 0
                    out_dict[i]["subwell"][j] = {
                        "out_R1":f'{raw_sample}/{i}/well/{j}_R1.fastq.gz',
                        "out_R2":f'{raw_sample}/{i}/well/{j}_R2.fastq.gz'
                    }
        
        # open output file
        out_fq1 = {}
        out_fq2 = {}
        for i in out_dict.keys():
            out_fq1[i] = {"subsample":{}}
            out_fq1[i]["subsample"] = gzip.open(out_dict[i]["subsample"]["out_R1"],"wb")
            out_fq2[i] = {"subsample":{}}
            out_fq2[i]["subsample"] = gzip.open(out_dict[i]["subsample"]["out_R2"],"wb")
            if self.args.split_to_well:
                out_fq1[i]["subwell"]={}
                out_fq2[i]["subwell"]={}
                for j in out_dict[i]['subwell'].keys():
                    out_fq1[i]["subwell"][j]= gzip.open(out_dict[i]["subwell"][j]["out_R1"],"wb")
                    out_fq2[i]["subwell"][j]= gzip.open(out_dict[i]["subwell"][j]["out_R2"],"wb")
        
        for i in range(self.fq1_number):
            fq1 = pyfastx.Fastx(self.fq1_list[i])
            fq2 = pyfastx.Fastx(self.fq2_list[i])
            for R1,R2 in zip(fq1,fq2):
                name1, seq1, qual1 = R1[0], R1[1], R1[2]
                name2, seq2, qual2 = R2[0], R2[1], R2[2]
                temp_bc = seq1[self.pattern_dict["C"][0].start:self.pattern_dict["C"][0].stop]

                for j in BC_dict.keys():
                    if temp_bc in BC_dict[j]["subBC"].keys():
                        out_dict[j]["read_num"][j] += 1
                        out_fq1[j]["subsample"].write(f'@{name1}\n{seq1}\n+\n{qual1}\n'.encode())
                        out_fq2[j]["subsample"].write(f'@{name2}\n{seq2}\n+\n{qual2}\n'.encode())
                        if self.args.split_to_well:
                            wellBC = BC_dict[j]["subBC"][temp_bc]
                            out_dict[j]["read_num"][wellBC] += 1
                            out_fq1[j]["subwell"][wellBC].write(f'@{name1}\n{seq1}\n+\n{qual1}\n'.encode())
                            out_fq2[j]["subwell"][wellBC].write(f'@{name2}\n{seq2}\n+\n{qual2}\n'.encode())
        # close files
        for i in out_fq1.keys():
            out_fq1[i]["subsample"].close()
            out_fq2[i]["subsample"].close()
            if self.args.split_to_well:
                for j in out_fq1[i]["subwell"]:
                    out_fq1[i]["subwell"][j].close()
                    out_fq2[i]["subwell"][j].close()
        
        out_json = raw_sample + ".bulk_rna.well_bc.json"
        utils.write_json(BC_dict,out_json)
        out_json = raw_sample + ".bulk_rna.fastq_inf.json"
        utils.write_json(out_dict,out_json)

        logger.info(out_dict)
        logger.info("Analysis finish!")


    def get_bc_dict(self,args):
        df = self.load_split_file(args.split_inf,args.sample)
        sub_dict = {}
        for _,row in df.iterrows():
            sub_sample = row["sub_sample"]
            sub_dict[sub_sample] =  {"subBC":{},"wellBC":{}}
            split_BC = utils.read_one_col(row['wellBC'])
            sub_dict[sub_sample]["subBC"] = parse_protocol.get_mismatch_dict(split_BC, 1)
            if args.split_to_well:
                for i in split_BC:
                    sub_dict[sub_sample]["wellBC"][i] = parse_protocol.get_mismatch_dict([i], 1)
        return sub_dict

    def load_split_file(self,file,sample):
        logger.info("infromation load!")
        df = pd.read_csv(file,sep=",",header=0)
        # file check
        if '_'.join(df.columns) != "sample_sub_sample_wellBC":
            sys.exit(f'Wrong file,{args.split_inf} header should be: sample,sub_sample,wellBC')
        if sample not in df["sample"].values:
            sys.exit(f"{sample} doesn't in split csv")

        df = df[df["sample"] == sample]

        for i in df["wellBC"]:
            if not os.path.exists(i):
                sys.exit(f"[Error]: {i}  not exists, please check!.\n")
        return df
 
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
    parser.add_argument("--pattern")
    parser.add_argument('--version', action='version', version='1.0')
    args = parser.parse_args()

    logger.info("Analysis Start!")
    runner = Split_Fastq(args)
    runner.run()   
