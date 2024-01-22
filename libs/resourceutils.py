import gffutils


def generate_intoron_gtf(db: gffutils.FeatureDB, output: str) -> None:
    introns = db.create_introns(exon_featuretype='exon', 
                                new_featuretype='intron', 
                                merge_attributes=True, 
                                numeric_sort=True)
    pybed = gffutils.pybedtools_integration.to_bedtool(introns)
    pybed.saveas(output)
    
    return