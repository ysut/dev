import pandas as pd
from .deco import print_filtering_count


def rename_like_vcf_format(df: pd.DataFrame) -> pd.DataFrame:
    df = df.rename(columns={'chr': 'CHROM', 'pos': 'POS', 
                            'ref': 'REF', 'alt': 'ALT'})

    return df


def generate_variant_id_col(row):
    """Generate variant ID column from CHROM, POS, REF, and ALT

    Args:
        row (pd.DataFrame): Using with pandas.DataFrame.apply

    Returns:
        str: Return variant ID (e.g. 1:23456-78910-A-T)
    """
    variant_id = f'{row["CHROM"]}:{row["POS"]}-{row["REF"]}-{row["ALT"]}'

    return variant_id


def variant_id_ck(df: pd.DataFrame, variant_id: str) -> None:
    if variant_id in df.columns:
        pass
    else:
        df['varinat_id'] = df['CHROM'] + ':' + df['POS'] + '-' + df['REF'] + '-' + df['ALT']
        pass
    
    return df
     
############ Functions for cleansing and adjusting HGMD data ############
def adjust_enst_for_hgmd(df: pd.DataFrame) -> pd.DataFrame:
    result = df.replace(
        {'ENST': {'ENST00000263201': 'ENST00000437685'},
         'ENST_Full': {'ENST00000263201.7_4': 'ENST00000437685.6_1',
                       'ENST00000361547.7_7': 'ENST00000361547.7_8',
                       'ENST00000609375.1_7': 'ENST00000347364.7_5',
                       'ENST00000649912.1_4': 'ENST00000347364.7_5'}})
    return result


@print_filtering_count
def remove_unkown_refalt(df: pd.DataFrame) -> pd.DataFrame:
    result = df.dropna(subset='REF', axis=0)
    return result
