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

### Logging setup
from logging import getLogger, config
import yaml
parent_directory = os.path.dirname(os.path.dirname('__file__'))
config_path: str = os.path.join(parent_directory, '../../../config/logging.yaml')
with open(config_path, 'r') as f:
    config.dictConfig(yaml.safe_load(f))
logger = getLogger(__name__)


# Set the path of CCRs
autoccr = '../../../Resources/04_CCRs/ccrs.autosomes.v2.20180420.bed.gz'
xccr = '../../../Resources/04_CCRs/ccrs.xchrom.v2.20180420.bed.gz'

# canonlist = '../../Resources/01_CanonicalTranscripts/CanonicalTranscripts.exoncount.tsv'
# canon = pd.read_table(canonlist, sep='\t', header=0)


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
        logger.debug(f"Warning: CDS_Length == 0 in {row['variant_id']}")
        return False
    
    shorten_len = skipped + deleted
    shorten_parcent = shorten_len / float(row['CDS_Length'])
    if shorten_parcent > 0.1:
        return True
    else:
        return False


# Determine if the gene is included in eLoFs genes

import re
elofs = pd.read_table(
	'../../../Resources/02_EstimatedLoFGenes/Final_eLoF_genes_list/Supplementary_Data_Estimated_LoF_genes.tsv', 
	usecols=['HGNC_ID'], sep='\t')

elofs_hgnc_ids_with_prefix = elofs['HGNC_ID'].unique().tolist()
elofs_hgnc_ids = [re.sub('HGNC:', '', hgnc) for hgnc in elofs_hgnc_ids_with_prefix]

def elofs_judge(row):
    if row['HGNC_ID'] in elofs_hgnc_ids:
        return True
    else:
        return False


# Determine causing NMD or escape NMD

def nmd_judge(row):
    if row['EXON']:
        max_exon: int = int(row['EXON'].split('/')[1])
        max_intron: int = max_exon - 1
    elif row['INTRON']:
        max_intron: int = int(row['INTRON'].split('/')[1])
    else:
        max_intron: int = -1
	
    try:
        curt_int = row['ExInt_INFO']['curt_Int']
    except:
        return 'Exonic (Non-Canonical)'
    else:
        # query_enst = row['ENST_Full']

        if max_intron == -1:
            return '[Warning] No intron info'
        # try:
        #     max_exon = canon.loc[canon['ENST_Full'] == query_enst, 'MaxExon'].values[0]
        # except:
        else:
            if curt_int == max_intron:
                return 'Escape_NMD'
            elif curt_int > max_intron:
                return '[Warning] current_int > max_intron'
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
    
    # Create skipped or deleted region bed file
    df_skip = df[df['skipped_region'].notnull()].copy()
    df_del = df[df['deleted_region'].notnull()].copy()
    
    df_skip['skipped_region'] = df['skipped_region'].replace(
        'Cannot predict splicing event', np.nan).copy()
    df_del['deleted_region'] = df['deleted_region'].replace(
        'Cannot predict splicing event', np.nan).copy()
    sr_skip = df_skip['skipped_region']
    sr_del = df_del['deleted_region']
    
    sr = pd.concat([sr_skip, sr_del], axis=0)
    sr.dropna(inplace=True)
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


