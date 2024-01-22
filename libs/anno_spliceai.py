import re
import numpy as np
import pandas as pd
from . import utils, preprocess


def load_result_splai(vcf: str) -> pd.DataFrame:
    vcf_columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
    splai = pd.read_table(vcf, sep='\t', header=None, names=vcf_columns)
 
    return splai


def _format_splai(splai: pd.DataFrame) -> pd.DataFrame:
    splai = pd.concat(
        [splai[['ID']], splai['INFO'].str.split('=', expand=True)], axis=1)
    sr = splai[1].str.contains(',')
    splai = pd.concat([splai, sr], axis=1)
    splai.columns = ['ID', 'dummy', 'INFO', 'is_Multi']
    splai = pd.concat(
        [splai[['ID', 'is_Multi']], splai['INFO'].str.split(',', expand=True)], 
        axis=1)

    return splai


def _extract_splai_result(row, genecol: str, colsnum: int):
    if row['is_Multi']:
        for i in range(colsnum):
            info = row[i]
            try:
                info: str = re.sub(r'\w+\|', '', info, 1)
            except:
                pass
            else:
                gene: str = re.match(r'[^|]+', info).group()
                if row[genecol] == gene:
                    # print(f'column: {i}, {gene}')
                    return info
                else:
                    pass
    else:
        try:
            info: str = re.sub(r'\w+\|', '', row[0], 1)
        except:
            return 'No SpliceAI predictions'
        else:
            return info


def anno(df: pd.DataFrame, splai: pd.DataFrame, genecol: str) -> pd.DataFrame:    
    if splai.loc[1, 'ID'] == '.':
        splai['ID'] = splai.parallel_apply(preprocess.generate_variant_id_col, 
                                           axis=1)
    else:
        pass

    splai = _format_splai(splai)

    df = pd.merge(df, splai, how='left', 
                  left_on='variant_id', right_on='ID')
    df = df.drop_duplicates()

    droplist = [s for s in splai.columns if isinstance(s, int)]
    df['SpliceAI'] = df.apply(_extract_splai_result, 
                              genecol=genecol, colsnum=len(droplist) ,axis=1)
    # df = df.drop(
    #     [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14 ,15, 16, 17], axis=1)
    df = df.drop(droplist, axis=1)
    df = pd.concat([df, df['SpliceAI'].str.split('|', expand=True)], 
                    axis=1).drop(0, axis=1)
    
    df.rename(columns={1: 'DS_AG', 2: 'DS_AL', 3: 'DS_DG', 4: 'DS_DL',
                       5: 'DP_AG', 6: 'DP_AL', 7: 'DP_DG', 8: 'DP_DL'}, 
                       inplace=True)
    
    return df


def _extract_maxsplai(row):
    maxvalue = max(
        row['DS_AG'], row['DS_AL'], row['DS_DG'], row['DS_DL'])
    
    return maxvalue


def insert_maxsplai(df: pd.DataFrame) -> pd.DataFrame:
    df = df.astype({'DS_AG': float, 'DS_AL': float, 
                    'DS_DG': float, 'DS_DL': float})

    df['maxsplai'] = df.parallel_apply(_extract_maxsplai, axis=1)

    return df


def extract_splai_result_lite(row):
    for i in range(9):  
        info: str = row[i]
        if info:
            info: str = info[2:]
            gene: str = re.match(r'[^|]+', info).group()
            if row['gene'] == gene:
                return info
            else:
                pass
        else:
            pass
    return 'No SpliceAI predictions'
