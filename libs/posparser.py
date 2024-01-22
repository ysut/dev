import re
import numpy as np
import pandas as pd
import gffutils
# import gffutils.pybedtools_integration
import pysam
# import pybedtools 
# import swifter
from Bio.Seq import Seq
from liftover import get_lifter



############ Functions for analysis ############
def classifying_canonical(df: pd.DataFrame, cdot: str) -> pd.DataFrame:
   df['is_Canonical'] = 'False'
   df['Int_loc'] = df[cdot].str.extract('([+-]\d+)')
   df.loc[(df['Int_loc'] == '-2') 
          | (df['Int_loc'] == '-1') 
          | (df['Int_loc'] == '+1') 
          | (df['Int_loc'] == '+2'), 
          'is_Canonical'] = 'True'
   canonicallist = []
   
   for d in ['-2', '-1', '+1', '+2']:
      cano_count = len(df.loc[df['Int_loc'] == d])
      print(f"{d}: {cano_count}")
      canonicallist.append(cano_count)

   print(f'Total variants      : {len(df)}')
   print(f'Canonical variants  : {sum(canonicallist)}')
   print(f'non-Canon variants  : {len(df) - sum(canonicallist)}\n')
   result = df.fillna({'Int_loc': 'Exonic'})

   return result


############ Functions for apply method ############
def calc_exon_loc(row, 
                  tabixfile: pysam.pysam.libctabix.TabixFile, 
                  enstcolname: str):
    """Calculate exon location from start of exon or end of exon.

    Args:
        row (pd.Series)  : _description_
        tabixfile        : 
        enstcolname (str): The name of ENST ID col

    Returns:
        str: "Upstream distance : Downstream distance"
    """
    query_chr: str = f"chr{row['CHROM']}"
    query_pos: int = int(row['POS'])
    query_start: int = int(query_pos) - 1
    query_end: int = int(query_pos)
    query_enst: int = row[enstcolname]

    for r in tabixfile.fetch(query_chr, query_start, query_end, 
                            parser=pysam.asGFF3()):
        
        try:
            enst = re.match(r'ENST\d+', r.transcript_id).group()
        except KeyError:
            pass
        else:
            if enst == query_enst:
                if r.feature == 'exon':
                    if r.strand == '+':
                        upd = r.end - query_pos
                        downd = query_pos - (r.start + 1)
                    elif r.strand == '-':
                        upd = query_pos - (r.start + 1)
                        downd = r.end - query_pos
                    else:
                        return 'unk_strand'
                    return f'{upd}:{downd}'
                else:
                    pass
            else:
                pass


def extract_splicing_region(row):
    if row['ex_up_dist'] is None:
        pass
    else: 
        if int(row['ex_down_dist']) == 0:
            return 'ex_acceptor_site'
        elif int(row['ex_up_dist']) <= 2:
            return 'ex_donor_site'
        else:
            return 'non_SplExon'

def select_donor_acceptor(row):
    """Select Donor site or Acceptor site

    Args:
        row (pd.DataFrame): Required columns are 
                            'Int_loc', 'ex_down_dist', and 'ex_up_dist'.
    Returns:
        str: 'Donor' or 'Acceptor' with 'ex' or 'int'
    """    
    if row['Int_loc'] == 'Exonic':
        try:
            int(row['ex_down_dist'])
        except TypeError:
            return 'Unknown'
        else:               
            d = int(row['ex_down_dist'])
            u = int(row['ex_up_dist'])
            if d < u:
                return 'Acceptor_ex'
            elif d > u:
                return 'Donor_ex'
            else:
                return 'Center_of_Exon'
    elif int(row['Int_loc']) < 0:
        return 'Acceptor_int'
    elif int(row['Int_loc']) > 0:
        return 'Donor_int'
    else:
        return 'Unkown'


def extract_splai_result(row):
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


def extract_splai_result_2(row, genecol: str):
    if row['is_Multi']:
        for i in range(17):
            info = row[i]
            # print(info)
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


def calc_ex_int_num(row, 
                    db: gffutils.interface.FeatureDB,
                    db_intron: gffutils.interface.FeatureDB):
    # print(f'{row["ENST_Full"]}-{row["gene"]}:{row["variant_id"]}:{row["Int_loc"]}')
    if (row['SpliceType'] == 'Donor_int' 
        or row['SpliceType'] == 'Acceptor_int'):

        introns = db_intron.children(row['ENST_Full'], featuretype='intron')
        max_intron = 0
        intron_num = 0
        while 1:
            try:
                intron = next(introns)
            except StopIteration:
                break
            else:
                max_intron += 1
                if (intron.start <= int(row['POS']) and int(row['POS']) <= intron.end):
                    intron_num = intron.attributes['exon_number'][0]
                else:
                    pass
        return f'{intron_num}/{max_intron}'
            
    elif (row['SpliceType'] == 'Donor_ex' 
          or row['SpliceType'] == 'Acceptor_ex'):
        
        exons = db.children(row['ENST_Full'], featuretype='exon')
        max_exon = 0
        exon_num = 0
        while 1:
            try:
                exon = next(exons)
            except StopIteration:
                break
            else:
                max_exon += 1
                if (exon.start <= int(row['POS']) and int(row['POS']) <= exon.end):
                    exon_num = exon.attributes['exon_number'][0]
                else:
                    pass
        return f'{exon_num}/{max_exon}'
    
    else:
        return 'unknown'
            

def select_exon_pos(row):
    if row['ex_up_dist']:
        return min(int(row['ex_up_dist']), int(row['ex_down_dist']))
    else:
        pass


def extract_splicing_region(row):
    if row['ex_up_dist'] is None:
        pass
    else: 
        if int(row['ex_down_dist']) == 0:
            return 'ex_acceptor_site'
        elif int(row['ex_up_dist']) <= 2:
            return 'ex_donor_site'
        else:
            return 'non_SplExon'


def select_donor_acceptor(row):
    """Select Donor site or Acceptor site

    Args:
        row (pd.DataFrame): Required columns are 
                            'Int_loc', 'ex_down_dist', and 'ex_up_dist'.
    Returns:
        str: 'Donor' or 'Acceptor' with 'ex' or 'int'
    """    
    if row['Int_loc'] == 'Exonic':
        try:
            int(row['ex_down_dist'])
        except TypeError:
            return 'Warning: Unknown'
        else:               
            d = int(row['ex_down_dist'])
            u = int(row['ex_up_dist'])
            if d < u:
                return 'Acceptor_ex'
            elif d > u:
                return 'Donor_ex'
            else:
                return 'Center_of_Exon'
    elif int(row['Int_loc']) < 0:
        return 'Acceptor_int'
    elif int(row['Int_loc']) > 0:
        return 'Donor_int'
    else:
        return 'Warning: Unkown'
    

def calc_prc_exon_loc(row):  
    if row['Int_loc'] == 'Exonic':
        try:
            curt_ex_length = int(row['ex_up_dist']) + int(row['ex_down_dist'])
        except TypeError as e:
            return np.nan
        
        try:
            row['ExInt_INFO']['strand']
        except:
            return np.nan


        if (row['ExInt_INFO']['strand'] == '+' 
                and int(row['ex_up_dist']) <= int(row['ex_down_dist'])):
            return int(row['exon_pos']) / curt_ex_length * 100
        
        elif (row['ExInt_INFO']['strand'] == '-' 
              and int(row['ex_up_dist']) >= int(row['ex_down_dist'])):
            return int(row['exon_pos']) / curt_ex_length * 100
        
        elif (row['ExInt_INFO']['strand'] == '+' 
              and int(row['ex_up_dist']) > int(row['ex_down_dist'])):
            return (1 - (int(row['exon_pos']) / curt_ex_length)) * 100
    
        elif (row['ExInt_INFO']['strand'] == '-' 
              and int(row['ex_up_dist']) < int(row['ex_down_dist'])):
            return (1 - (int(row['exon_pos']) / curt_ex_length)) * 100
    
        else:
            return 'Waring: unexpected conditions'
    
    else:
        return 'Intron_variant'

# def anno_strand_old(row, 
#                     db: gffutils.interface.FeatureDB,
#                     db_intron: gffutils.interface.FeatureDB):
#     """Calculate exon location from start of exon or end of exon.

#     Args:
#         row (_type_): _description_

#     Returns:
#         str: "Upstream distance : Downstream distance"
#     """
#     tbx_anno = tabixfile
#     query_chr: str = f"chr{row['CHROM']}"
#     query_pos: int = int(row['POS'])
#     query_start: int = int(query_pos) - 1
#     query_end: int = int(query_pos) + 1
#     query_enst: str = row[queryenst]

#     for r in tbx_anno.fetch(query_chr, query_start, query_end, 
#                             parser=pysam.asGFF3()):
#         try:
#             enst = re.match(r'ENST\d+', r.transcript_id).group()
#         except KeyError:
#             pass
#         else:
#             if enst == query_enst:
#                 return r.strand
#             else:
#                 return 'unk'


# def anno_strand_old(row, tabixfile: pysam.pysam.libctabix.TabixFile, queryenst: str):
#     """Calculate exon location from start of exon or end of exon.

#     Args:
#         row (_type_): _description_

#     Returns:
#         str: "Upstream distance : Downstream distance"
#     """
#     tbx_anno = tabixfile
#     query_chr: str = f"chr{row['CHROM']}"
#     query_pos: int = int(row['POS'])
#     query_start: int = int(query_pos) - 1
#     query_end: int = int(query_pos) + 1
#     query_enst: str = row[queryenst]

#     for r in tbx_anno.fetch(query_chr, query_start, query_end, 
#                             parser=pysam.asGFF3()):
#         try:
#             enst = re.match(r'ENST\d+', r.transcript_id).group()
#         except KeyError:
#             pass
#         else:
#             if enst == query_enst:
#                 return r.strand
#             else:
#                 return 'unk'
            

