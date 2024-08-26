def multiqc_sgr_config():
    from multiqc import config
    
    """ Set up MultiQC config defaults for this package """
    sgr_search_patterns = {
        "bulk_rna/stats": {"fn": "*bulk_rna.*stats.json"},
        "bulk_rna/well_count": {"fn": "*bulk_rna.counts_report.json"},
    }
    config.update_dict(config.sp, sgr_search_patterns)