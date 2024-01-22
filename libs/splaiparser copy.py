import numpy as np
import pandas as pd
import gffutils


def calc_exint_info(row, db, db_intron):
    ## Set up variables
    query_enst = row['ENST_Full'] 
    chrom, pos = f'chr{row["CHROM"]}', int(row['POS'])

    ## Strand check
    strand = next(db.children(query_enst, featuretype='transcript')).strand
    if strand == '+':
        region: tuple = (chrom, pos-1, pos)
    elif strand == '-':
        region: tuple = (chrom, pos, pos+1)
    else:
        print(f'Warning: unkown strand -> {row["ENST_Full"]}')
        region: tuple = (chrom, pos-1, pos)

    ## Fetch exon or intron information from each GENCODE DBs
    try:
        fetched_data = db.children(query_enst, 
                                   limit=region ,featuretype='exon')
        d = next(fetched_data)
    except StopIteration: 
        fetched_data = db_intron.children(query_enst, 
                                          limit=region ,featuretype='intron')
        d = next(fetched_data)
    else:
        pass
    
    ## Set attributes and current featuretype
    d_attr: list = d.attributes
    curtFeature = d.featuretype

    ## This step is divided into two parts (Exon or Intron)
    if curtFeature == 'exon':
        #1. Set current exon coordinates
        curtExNum = int(d_attr['exon_number'][0])
        curtExStart, curtExEnd = d.start, d.end
        
        #2. Set previous exon coordinates
        if curtExNum > 1:
            exons = db.children(query_enst, featuretype='exon')
            for e in exons:
                if int(e.attributes['exon_number'][0]) == curtExNum - 1:
                    prevExStart, prevExEnd = e.start, e.end
                    break
                else:
                    pass
        else:
            prevExStart, prevExEnd = '1st_Exon', '1st_Exon'

        #3. Set next exon coordinates
        exons = db.children(query_enst, featuretype='exon')
        nextExStart, nextExEnd = 'Last_Exon', 'Last_Exon'
        for e in exons:
            if int(e.attributes['exon_number'][0]) == curtExNum + 1:
                nextExStart, nextExEnd = e.start, e.end
                break
            else:
                pass

        #4. Set eStart & eEnd
        eStart = curtExStart
        eEnd = curtExEnd

    elif curtFeature == 'intron':
        #1. Set current intron coordinates
        curtIntNum = int(d_attr['exon_number'][0])
        curtIntStart, curtIntEnd = d.start, d.end

        #2. Set previous & next exon coordinates
        exons = db.children(query_enst, featuretype='exon')
        for e in exons:
            if int(e.attributes['exon_number'][0]) == curtIntNum:
                prevExStart, prevExEnd = e.start, e.end
            elif int(e.attributes['exon_number'][0]) == curtIntNum + 1:
                nextExStart, nextExEnd = e.start, e.end
            else:
                pass 
        
        #3. Check Acceptor site or Donor site and return close Exon info
        up = pos - d.start + 1
        down = d.end - pos + 1
        if (((up > down) & (strand == '+')) 
            |((up < down) & (strand == '-'))): 
                eStart = nextExStart
                eEnd = nextExEnd
        elif (((up < down) & (strand == '+')) 
            | ((up > down) & (strand == '-'))):
                eStart = prevExStart
                eEnd = prevExEnd
        else: # center of intron
                eStart = 'unk'
                eEnd = 'unk'
            
    else:
        print(f'Warning: Not exon or intron -> {row["ID"]}')
        pass

    ## Return results as dict type
    if d.featuretype == 'exon':
        results = {'strand': strand, 
                   'eStart': eStart, 
                   'eEnd': eEnd,
                   'curt_Ex': curtExNum, 
                   'curt_ExStart': curtExStart,
                   'curt_ExEnd': curtExEnd,
                   'prev_Ex': int(curtExNum - 1),
                   'prev_ExStart': prevExStart,
                   'prev_ExEnd': prevExEnd,
                   'next_Ex': int(curtExNum + 1), 
                   'next_ExStart': nextExStart,
                   'next_ExEnd': nextExEnd}

    else:
        results = {'strand': strand,
                   'eStart': eStart,
                   'eEnd': eEnd,
                   'curt_Int': curtIntNum, 
                   'curt_IntStartEnd': curtIntStart, 
                   'curtIntEnd': curtIntEnd,
                   'prev_Ex': curtIntNum, 
                   'prev_ExStart': prevExStart,
                   'prev_ExEnd': prevExEnd,
                   'next_Ex': int(curtIntNum + 1), 
                   'next_ExStart': nextExStart,
                   'next_ExEnd': nextExEnd}
    
    return results


#1.   Calculate gained exon size for pusedoexon activation
#1-1. Filters
def _filtering_DS_Loss_threshold(thresholds: dict, **kwargs):
    min_sALDL = min(kwargs['DS_AL'], kwargs['DS_DL'])
    max_sALDL = max(kwargs['DS_AL'], kwargs['DS_DL'])

    if ((min_sALDL >= thresholds['TH_min_sALDL']) 
        & (max_sALDL >= thresholds['TH_max_sALDL'])):
        return 'PASS'
    else:
        return 'FAIL'

def _filtering_DS_Gain_threshold(thresholds: str, **kwargs):
    min_sAGDG = min(kwargs['DS_AG'], kwargs['DS_DG'])
    max_sAGDG = max(kwargs['DS_AG'], kwargs['DS_DG'])

    if ((min_sAGDG >= thresholds['TH_min_sAGDG'])
        & (max_sAGDG >= thresholds['TH_max_sAGDG'])):
        return 'PASS'
    else:
        return 'FAIL'

def _is_partial_effect(**kwargs):
    strand = kwargs['ExInt_INFO']['strand']
    pAG, pDG = kwargs['DP_AG'], kwargs['DP_DG']

    if ((strand == '+') & (pAG < pDG)) | ((strand == '-') & (pAG > pDG)):
        return True
    else:
        return False

def _calc_gained_exon_size(thresholds: dict, **kwargs):
    if ((_filtering_DS_Gain_threshold(thresholds, **kwargs) == 'PASS')
        & (_is_partial_effect(**kwargs))):
        return np.abs(int(kwargs['DP_AG']) - int(kwargs['DP_DG'])) + 1
    else:
        return None

#.1-2 Verify the pseudoexon location
def _verify_pseudoexon_location(db_intron, **kwargs):
    """Verify the pseudoexon location
    - Both Acceptor gain site and Donor gain site 
      are located in the same intron.
    - AG site is located at >50 bp from the start of intron.
    - DG site is located at >50 bp from the end of intron.
    """

    posAG, posDG = int(kwargs['DP_AG']), int(kwargs['DP_DG'])
    introns = db_intron.children(kwargs['ENST_Full'], featuretype='intron')

    for i in introns:
        if i.start < posAG < i.end:
            if (i.strand == '+'
                and i.start + 50 < posAG
                and i.end - 50 > posDG):
                return True

            elif (i.strand == '-'
                and i.start + 50 < posDG
                and i.end - 50 > posAG):
                return True

            else:
                return False
        else:
            return False


##. Validate cryptic splice site activation for partial deletion or retention
def _is_cryptic_Acp_activation(thresholds: dict, **kwargs):
    sAG, sDG = kwargs['DS_AG'], kwargs['DS_DG']

    if ((sAG >= thresholds['TH_sAG']) & (sAG > sDG)):
        return True
    else:
        return False
    
def _is_cryptic_Dnr_activation(thresholds: dict, **kwargs):
    sAG, sDG = kwargs['DS_AG'], kwargs['DS_DG']

    if ((sDG >= thresholds['TH_sDG']) & (sAG > sDG)):
        return True
    else:
        return False
    

##. Orientation filters
def _filtering_Acp_orientation(**kwargs): 
    # 1-based
    posAG: int = kwargs['POS'] + int(kwargs['DP_AG'])
    info: dict = kwargs['ExInt_INFO']
    strand = info['strand']
    prevExStart = info['prev_ExStart']
    prevExEnd = info['prev_ExEnd']

    if prevExStart == '1st_Exon':
        return '1st_Exon'
    else:
        pass

    if (((strand == '+') & (int(prevExEnd) < posAG))
        |((strand == '-') & (posAG < int(prevExStart)))):
        return 'PASS'
    else:
        return 'FAIL'

def _filtering_Dnr_orientation(**kwargs):
    # 1-based
    posDG: int = kwargs['POS'] + int(kwargs['DP_DG'])
    info: dict = kwargs['ExInt_INFO']
    strand = info['strand']
    nextExStart: int = info['next_ExStart']
    nextExEnd: int = info['next_ExEnd']

    if nextExStart == 'Last_Exon':
        return 'Last_Exon'
    else:
        pass
    
    if (((strand == '+') & (posDG < int(nextExStart)))
        |((strand == '-') & (int(nextExEnd) < posDG))):
        return 'PASS'
    else:
        return 'FAIL'


##. Predicted changed exon size in 5-prime side and 3-prime side
def _bp_5prime(thresholds: str, **kwargs) -> int:
    posAG: int = kwargs['POS'] + int(kwargs['DP_AG'])
    info: dict = kwargs['ExInt_INFO']
    strand: str = info['strand']
    eStart: int = info['eStart']
    eEnd: int = info['eEnd']

    if ((strand == '+') 
        & (_is_cryptic_Acp_activation(thresholds, **kwargs))
        & (_filtering_Acp_orientation(**kwargs) == 'PASS')):
        return posAG - eStart # 1-based
    elif ((strand == '-') 
        & (_is_cryptic_Acp_activation(thresholds, **kwargs))
        & (_filtering_Acp_orientation(**kwargs) == 'PASS')):
        return eEnd - posAG # 1-based
    else:
        return 0

def _bp_3prime(thresholds: str, **kwargs) -> int:
    posDG: int = kwargs['POS'] + int(kwargs['DP_DG'])
    info: dict = kwargs['ExInt_INFO']
    strand: int = kwargs['ExInt_INFO']['strand']
    eStart: int = info['eStart']
    eEnd: int = info['eEnd']

    if ((strand == '+') 
        & (_is_cryptic_Dnr_activation(thresholds, **kwargs))
        & (_filtering_Dnr_orientation(**kwargs) == 'PASS')):
        return posDG - eEnd # 1-based
    elif ((strand == '-') 
        & (_is_cryptic_Dnr_activation(thresholds, **kwargs))
        & (_filtering_Dnr_orientation(**kwargs) == 'PASS')):
        return eStart - posDG # 1-based
    else:
        return 0

        

##. Evaluate orientation and classify Lost exon or Reteined intron
def _classify_LEX_RIT(**kwargs):
    strand = kwargs['ExInt_INFO']['strand']
    pAL, pDL = kwargs['DP_AL'], kwargs['DP_DL']

    if ((strand == '+') & (pAL < pDL)) | ((strand == '-') & (pAL > pDL)):
        return 'LEX'
    else:
        return 'RIT'


##. Varidate variant position from close exon boundary (50 bp or 250 bp) 
def _calc_dist_from_exon(**kwargs):
    pos = kwargs['POS']
    eStart, eEnd = kwargs['ExInt_INFO']['eStart'], kwargs['ExInt_INFO']['eEnd']
    dist_exon_start: int = pos - eStart
    dist_exon_end: int = pos - eEnd
    if ((dist_exon_start <= 0) & (dist_exon_end < 0)):
        return np.abs(dist_exon_start)
    elif ((dist_exon_start > 0) & (dist_exon_end >= 0)):
        return np.abs(dist_exon_end)
    else:
        return 0

def _varidate_var_pos_250bp(**kwargs):
    if _calc_dist_from_exon(**kwargs) > 250:
        return 'outside_250bp'
    else:
        return 'within_250bp'

def _varidate_var_pos_50bp(**kwargs):
    if _calc_dist_from_exon(**kwargs) > 50:
        return 'outside_50bp'
    else:
        return 'within_50bp'


##. Predictions
def predict_gained_exon(thresholds: dict, **kwargs):
    gained_exon_size = _calc_gained_exon_size(thresholds, **kwargs)
    if gained_exon_size:
        if ((gained_exon_size >= thresholds['TH_min_GExon']) 
            & (gained_exon_size <= thresholds['TH_max_GExon'])):
            return True
        else:
            return False
    else:
        return False

def predict_lost_exon(thresholds: dict, **kwargs):
    if ((_filtering_DS_Loss_threshold(thresholds, **kwargs) == 'PASS') 
        & (_classify_LEX_RIT(**kwargs) == 'LEX')):
        return np.abs(int(kwargs['DP_AL'])- int(kwargs['DP_DL'])) + 1
    else:
        return None

def predict_retein_intron(**kwargs):
    if ((_filtering_DS_Loss_threshold(**kwargs) == 'PASS') 
        & (_classify_LEX_RIT(**kwargs) == 'RIT')):
        return np.abs(int(kwargs['DP_DL']) - int(kwargs['DP_AL'])) - 1
    else:
        return None



################################################################################
##                          Summrize splicing events                          ##
################################################################################

def pseudoexon_activation(row, thresholds, db_intron):
    if (_varidate_var_pos_50bp(**row) == 'outside_50bp'
        and predict_gained_exon(thresholds=thresholds, **row)
        and _calc_gained_exon_size(thresholds=thresholds, **row)
        and _verify_pseudoexon_location(db_intron=db_intron, **row)):   
        return True
    else:
        return False
    

def partial_intron_retention(row, thresholds):
    if ((_varidate_var_pos_250bp(**row) == 'within_250bp')
        and (_is_cryptic_Acp_activation(thresholds=thresholds, **row)) 
             or (_is_cryptic_Dnr_activation(thresholds=thresholds, **row))
        and (-251 < _bp_5prime(thresholds, **row) < 0) 
             or (0 < _bp_3prime(thresholds, **row) < 251)):
        return True
    else:
        return False


def partial_exon_deletion(row, thresholds):
    if ((_varidate_var_pos_250bp(**row) == 'within_250bp')
        and ((_bp_5prime(thresholds, **row) > 0) 
             or (_bp_3prime(thresholds, **row) < 0))):
        return True
    else:
        return False


def exon_skipping(row, thresholds):
    """Varidate exon skipping
    When the variant is located outside >50 bp from 
    closest exon-intron boundary, exon skipping may not occur.
    """
    lost_exon_size = predict_lost_exon(thresholds=thresholds, **row)
    if ((_varidate_var_pos_50bp(**row) == 'outside_50bp')
        or (lost_exon_size is None)):
        return False
    elif ((_varidate_var_pos_50bp(**row) == 'within_50bp')
          and (lost_exon_size)):            
        info = row['ExInt_INFO']
        native_exon_length = int(info['eEnd']) - int(info['eStart']) + 1
        if lost_exon_size == native_exon_length:
            return native_exon_length
        else:
            return False
    else:
        return False
    

def intron_retention(row, thresholds):
    """Varidate intron retention
    When the variant is located outside >50 bp from 
    close exon-intron boundary, intron retention may not occur.
    """
    if ((_varidate_var_pos_50bp(**row) == 'outside_50bp' 
        or predict_retein_intron(thresholds=thresholds, **row) is None)):
        return False
    elif ((_varidate_var_pos_50bp(**row) == 'within_50bp' 
        or predict_retein_intron(thresholds=thresholds, **row))):
        return True
    else:
        return False


################################################################################
##                   Calculate Aberrant splicing event size                   ##
################################################################################

def anno_gained_exon_size(row, thresholds):
    if row['Pseudoexon']:
        pseudoexon_size = _calc_gained_exon_size(thresholds=thresholds, **row)
        return pseudoexon_size
    else:
        return np.nan
    
def anno_skipped_exon_size(row, thresholds):
    if row['Exon_skipping']:
        info = row['ExInt_INFO']
        lost_exon_size = predict_lost_exon(thresholds=thresholds, **row)
        nativeExonLength = int(info['eEnd']) - int(info['eStart']) + 1

        if lost_exon_size == nativeExonLength:
            return nativeExonLength
        else:
            return np.nan
    else:
        return np.nan
        
def anno_skipped_exon_size(row, thresholds):
    if row['Exon_skipping']:
        info = row['ExInt_INFO']
        lost_exon_size = predict_lost_exon(thresholds=thresholds, **row)
        nativeExonLength = int(info['eEnd']) - int(info['eStart']) + 1

        if lost_exon_size == nativeExonLength:
            return nativeExonLength
        else:
            return np.nan
        

def anno_partial_deleted_exon_size(row, thresholds):
    if row['partial_exon_deletion']:
        return _bp_5prime(thresholds=thresholds, **row)


def _bp_5prime(thresholds: str, **kwargs) -> int:
    posAG: int = kwargs['POS'] + int(kwargs['DP_AG'])
    info: dict = kwargs['ExInt_INFO']
    strand: str = info['strand']
    eStart: int = info['eStart']
    eEnd: int = info['eEnd']

    if ((strand == '+') 
        & (_is_cryptic_Acp_activation(thresholds, **kwargs))
        & (_filtering_Acp_orientation(**kwargs) == 'PASS')):
        return posAG - eStart # 1-based
    elif ((strand == '-') 
        & (_is_cryptic_Acp_activation(thresholds, **kwargs))
        & (_filtering_Acp_orientation(**kwargs) == 'PASS')):
        return eEnd - posAG # 1-based
    else:
        return 0






