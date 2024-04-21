import os
import subprocess
import numpy as np
import pandas as pd
from pybedtools import BedTool
from pandarallel import pandarallel

########   Initialize and setup pandas methods   ########
pandarallel.initialize(nb_workers=os.cpu_count()-1, progress_bar=False, 
                       verbose=0, use_memory_fs=False) 
os.environ['JOBLIB_TEMP_FOLDER'] = '/tmp' 


# Set the path of CCRs
autoccr = '../../Resources/04_CCRs/ccrs.autosomes.v2.20180420.bed.gz'
xccr = '../../Resources/04_CCRs/ccrs.xchrom.v2.20180420.bed.gz'

canonlist = '../../Resources/01_CanonicalTranscripts/CanonicalTranscripts.exoncount.tsv'
canon = pd.read_table(canonlist, sep='\t', header=0)


# Calculate the length of CDS
def calc_cds_len(row, db) -> int:
    cds_length: int = 0
    query_enst = row['ENST_Full']
    for d in db.children(query_enst, featuretype='CDS'):
        cds = np.abs(d.end - d.start) + 1
        cds_length += cds
    
    return cds_length

def calc_cds_len_shorten(row):
    if row['Exon_skipping']:
        skipped = float(row['Size_skipped_exon'])
    elif row['Part_ExDel']:
        deleted = float(row['Size_Part_ExDel'])
    else:
        return False
    
    try:
        skipped
    except NameError:
        skipped = np.nan
    else:
        pass

    try:
        deleted
    except NameError:
        deleted = np.nan
    else:
        pass
    
    if row['CDS_Length'] == 0:
        print(f"Warning: CDS_Length == 0 in {row['variant_id']}")
        return False
    
    shorten_len = skipped + deleted
    shorten_parcent = shorten_len / float(row['CDS_Length'])
    if shorten_parcent > 0.1:
        return True
    else:
        return False


# Determine if the gene is included in eLoFs genes
elofs = pd.read_table('../../Resources/02_EstimatedLoFGenes/lof_genes.txt', 
                      header=None, names=['gene'], sep='\t')
elofs_genes = elofs['gene'].unique().tolist()

def elofs_judge(row):
    if row['gene'] in elofs_genes:
        return True
    else:
        return False


# Determine causing NMD or escape NMD

def nmd_judge(row):
    try:
        curt_int = row['ExInt_INFO']['curt_Int']
    # except KeyError:
    except:
        return 'Exonic(Non-Canonical)'
    else:
        query_enst = row['ENST_Full']
        try:
            max_exon = canon.loc[canon['ENST_Full'] == query_enst, 'MaxExon'].values[0]
        except:
            return 'No_intron_info'
        else:
            max_intron = max_exon - 1
            if curt_int == max_intron:
                return 'Escape_NMD'
            elif curt_int > max_intron:
                return 'Warning: current_int > max_intron'
            else:
                return 'Possibly_NMD'


# Determine inframe or frameshift
def frame_check(x):
    if np.isnan(x):
        return False
    else:   
        if x % 3 == 0:
            return False
        elif x % 3 != 0:
            return True
        else:
            print('Error: frame_check()')
            return False


def anno_ccr_score(df: pd.DataFrame) -> pd.DataFrame:
    def fetch_ccr_score(row, col):
        region = row[col]
        if isinstance(region, str):
            split_region = region.split(' ')
            region_tuple = (split_region[0], split_region[1], split_region[2])
            pass
        else:
            return np.nan
        
        return results_dict.get(region_tuple, None)
    
    df_skip = df[df['skipped_region'].notnull()] 
    df_del = df[df['deleted_region'].notnull()]
    sr_skip = df_skip['skipped_region']
    sr_del = df_del['deleted_region']

    sr = pd.concat([sr_skip, sr_del], axis=0)
    bedstr = '\n'.join(sr)
    all_regions = BedTool(bedstr, from_string=True)
    try:
        ccr_auto = BedTool(autoccr)
        ccr_x  = BedTool(xccr)
    except FileNotFoundError:
        subprocess.run(['bash', '../../Resources/04_CCRs/dlccrs.sh'])
        ccr_auto = BedTool(autoccr)
        ccr_x  = BedTool(xccr)

    intersected_auto = all_regions.intersect(ccr_auto, wa=True, wb=True)
    intersected_x = all_regions.intersect(ccr_x, wa=True, wb=True)

    results_dict = {}
    for feature in intersected_auto:
        query_key = (feature.fields[0], feature.fields[1], feature.fields[2])
        score = float(feature.fields[6])  # fetch prcent CCR

        if query_key not in results_dict or score > results_dict[query_key]:
            results_dict[query_key] = score

    for feature in intersected_x:
        query_key = (feature.fields[0], feature.fields[1], feature.fields[2])
        score = float(feature.fields[6])  # fetch prcent CCR

        if query_key not in results_dict or score > results_dict[query_key]:
            results_dict[query_key] = score

    df['skipped_ccrs'] = df.parallel_apply(
        fetch_ccr_score, col='skipped_region', axis=1)
    df['deleted_ccrs'] = df.parallel_apply(
        fetch_ccr_score, col='deleted_region', axis=1)

    return df


