import json
import logging
from collections import defaultdict

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import table

# Initialise the logger
log = logging.getLogger("multiqc")
ASSAY = "bulk_rna"

class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super().__init__(
            name = ASSAY,
            anchor = ASSAY,
        )
        log.info(f"Running module: {ASSAY}")
        
        stat_data = self.parse_json(ASSAY, "stats")
        well_data = self.parse_json(ASSAY, "well_count")
        if all(len(x) == 0 for x in [stat_data,well_data]):
            raise ModuleNoSamplesFound
        
        sample_list = list(stat_data.keys())
        sample_list.sort()
        # Basic Stats Table
        self.general_stats_table(stat_data)

        # well detail
        helptext = """
        The rowname of the table is the well barcode sequence number. Only output wells that meet the conditions(default: UMI>=500).  
        If no well pass the filter,output wells that UMI,read and gene > 0.
        """
        for sample in sample_list:
            self.add_section(name = f"{sample} - per well", anchor = f'{sample}_well_table',helptext=helptext,plot = self.well_table(well_data[sample],sample))

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        
        self.add_software_version(None)
    
    def parse_json(self, assay, seg):
        data_dict = defaultdict(dict)
        n = 0
        for f in self.find_log_files(f"{assay}/{seg}"):
            log.debug(f"Found file: {f['fn']}")
            n += 1
            parsed_data = json.loads(f["f"])
            if parsed_data is not None:
                x = f["s_name"]
                s_name = x[: x.find(f".{assay}")]
                if s_name in data_dict:
                    log.info(f"Duplicate sample name found! Update: {s_name}")
                self.add_data_source(f, s_name=s_name, section=seg)
                data_dict[s_name].update(parsed_data)

        data_dict = self.ignore_samples(data_dict)

        log.info(f"Found {n} {assay} {seg} reports")
        # Write parsed report data to a file
        self.write_data_file(data_dict, f"multiqc_{assay}_{seg}")
        return data_dict
        
    def general_stats_table(self, summary_data):
        headers = {
            "Protocol": {
                "title": "Protocol",
                "description": "Predefined pattern of barcode and UMI",
                "scale": "purple",
                "hidden": True
            },
            "Raw Reads": {
                "title": "Raw Reads",
                "description": "Number of reads in the input file",
                "scale": "blue",
                "format": "{:,.0f}",
                "hidden": True
            },
            "Valid Reads": {
                "title": "Valid Reads",
                "description": "Percent of reads with valid barcodes",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "green",
                "hidden": False
            },
            "Corrected Barcodes": {
                "title": "Corrected Barcodes",
                "description": "Percent of corrected barcodes",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "green",
                "hidden": True
            },
            "Reads Mapped To Unique Loci": {
                "title": "Unique Reads",
                "description": "Percent of reads mapped uniquely to genome",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "green",
                "hidden": False
            },
            "Reads Mapped To Multiple Loci": {
                "title": "Multi-Loci Reads",
                "description": "Percent of reads mapped to multiple loci",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "green",
                "hidden": True
            },
            "Reads Mapped Uniquely To Transcriptome": {
                "title": "Counted Unique Reads",
                "description": "Percent of reads mapped uniquely to transcriptome; These reads are used for UMI counting",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "green",
                "hidden": False
            },
            "Mapped Reads Assigned To Exonic Regions": {
                "title": "Exonic Reads",
                "description": "Percent of mapped reads assigned to exonic regions",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "green",
                "hidden": True
            },
            "Mapped Reads Assigned To Intronic Regions": {
                "title": "Intronic Reads",
                "description": "Percent of mapped reads assigned to intronic regions",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "green",
                "hidden": True
            },
            "Mapped Reads Assigned To Intergenic Regions": {
                "title": "Intergenic Reads",
                "description": "Percent of mapped reads assigned to intergenic regions",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "green",
                "hidden": True
            },
            "Mapped Reads Assigned Antisense To Gene": {
                "title": "Antisense Reads",
                "description": "Percent of mapped reads assigned antisense to gene",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "green",
                "hidden": True
            },
            "Median Reads per Well": {
                "title": "Median Reads",
                "description": "Median number of reads per well",
                "scale": "blue",
                "format": "{:,.0f}",
                "hidden": False
            },
            "Median UMI per Well": {
                "title": "Median UMI",
                "description": "Median number of umi per well",
                "scale": "blue",
                "format": "{:,.0f}"
            },
            "Median Genes per Well": {
                "title": "Median Genes",
                "description": "Median number of genes per well",
                "scale": "blue",
                "format": "{:,.0f}",
                "hidden": False
            },
            "Mean Reads per Well": {
                "title": "Mean Reads",
                "description": "Mean number of reads per well",
                "scale": "blue",
                "format": "{:,.0f}",
                "hidden": False
            },
            "Mean UMI per Well": {
                "title": "Mean UMI",
                "description": "Mean number of umi per well",
                "scale": "blue",
                "format": "{:,.0f}",
                "hidden": False
            },
            "Mean Genes per Well": {
                "title": "Mean Genes",
                "description": "Mean number of genes per well",
                "scale": "blue",
                "format": "{:,.0f}",
                "hidden": False
            }
        }
        self.general_stats_addcols(summary_data, headers=headers)

    def well_table(self,well_data,sample):
        table_config = {
            "id": f"{sample}_Well_table",
            "namespace": f"{sample}",
            "title": "Detailed information per Well",
            "col1_header": "Well",
            "sort_rows":False
        }
        headers = {
            "UMI": {
                "title": "UMI",
                'description':"UMI counts",
                "format": "{:,.0f}"
            },
            "read": {
                "title": "read",
                'description':"read counts",
                "format": "{:,.0f}"
            },
            "gene": {
                "title": "gene",
                'description':"gene counts",
                "format": "{:,.0f}"
            },
        }
        return table.plot(well_data,pconfig=table_config,headers=headers)