import numpy as np
import pandas as pd
import gffutils

def __exits_spliceai_scores(row):
    if row['MaxSpliceAI'] == 'NA' or row['MaxSpliceAI'] == '.':
        return False
    else:
        return True


def calc_exint_info(row, db, db_intron):
    ## Set up variables
    query_enst = row['ENST_Full'] 
    chrom, pos = f'chr{row["CHROM"]}', int(row['POS'])
    strand = row['Strand']

    ## Strand check
    # try:
    #     strand = next(db.children(query_enst, featuretype='transcript')).strand
    # except StopIteration:
    #     # print(f'Warning: No strand information -> {row["ENST_Full"]}')
    #     return 'Warning'
    # else:
    #     pass
    
    if strand == '+':
        region: tuple = (chrom, pos-1, pos)
        # region: tuple = (chrom, pos, pos+1)
    elif strand == '-':
        region: tuple = (chrom, pos, pos+1)
    else:
        print(f'Warning: unkown strand -> {row["ENST_Full"]}')
        region: tuple = (chrom, pos-1, pos)

    ## Fetch exon or intron information from each GENCODE DBs

    # Fixed
    # if ((row['SpliceType'] == 'Donor_ex') | (row['SpliceType'] == 'Acceptor_ex')):
    #     fetched_data = db.children(
    #         query_enst, limit=region, featuretype='exon')
    #     d = next(fetched_data)
    # else:
    #     fetched_data = db_intron.children(
    #         query_enst, limit=region, featuretype='intron')
    #     d = next(fetched_data)


    # Not Used
    try:
        fetched_data = db.children(
            query_enst, limit=region, featuretype='exon')
        d = next(fetched_data)
    except StopIteration: 
        try:
            fetched_data = db_intron.children(
                query_enst, limit=region, featuretype='intron')
            d = next(fetched_data)
        except StopIteration:
            # print(f'Warning: No exon or intron found -> {row["variant_id"]}: {row["ENST_Full"]}')
            return 'Warning'
        else:
            pass
    else:
        pass
    
    if row['is_Canonical'] == 'True':
        fetched_data = db_intron.children(
            query_enst, limit=region, featuretype='intron')
        d = next(fetched_data)

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
        # print(f'Warning: Not exon or intron -> {row["ID"]}:{row["ENST_Full"]}')
        return 'Warning'

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
                   'curt_IntStart': curtIntStart, 
                   'curt_IntEnd': curtIntEnd,
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
    min_sALDL = min(float(kwargs['DS_AL']), float(kwargs['DS_DL']))
    max_sALDL = max(float(kwargs['DS_AL']), float(kwargs['DS_DL']))

    if ((min_sALDL >= float(thresholds['TH_min_sALDL'])) 
        & (max_sALDL >= float(thresholds['TH_max_sALDL']))):
        # print('PASS')
        return 'PASS'
    else:
        return 'FAIL'

def _filtering_DS_Gain_threshold(thresholds: str, **kwargs):
    min_sAGDG = min(float(kwargs['DS_AG']), float(kwargs['DS_DG']))
    max_sAGDG = max(float(kwargs['DS_AG']), float(kwargs['DS_DG']))

    if ((min_sAGDG >= float(thresholds['TH_min_sAGDG']))
        & (max_sAGDG >= float(thresholds['TH_max_sAGDG']))):
        # print('PASS')
        return 'PASS'
    else:
        # print('FAIL')
        return 'FAIL'

def _is_partial_effect(**kwargs):
    strand = kwargs['ExInt_INFO']['strand']
    pAG, pDG = int(kwargs['DP_AG']), int(kwargs['DP_DG'])

    if ((strand == '+') & (pAG < pDG)) | ((strand == '-') & (pAG > pDG)):
        # print('partial effect')
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

    pAG, pDG = int(kwargs['DP_AG']), int(kwargs['DP_DG'])
    introns = db_intron.children(kwargs['ENST_Full'], featuretype='intron')
    posAG: int = int(kwargs['POS']) + pAG
    posDG: int = int(kwargs['POS']) + pDG

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
            pass

##. Validate cryptic splice site activation for partial deletion or retention
def _is_cryptic_Acp_activation(thresholds: dict, **kwargs):
    sAG, sDG = float(kwargs['DS_AG']), float(kwargs['DS_DG'])

    if ((sAG >= float(thresholds['TH_sAG'])) & (sAG > sDG)):
        return True
    else:
        return False
    
def _is_cryptic_Dnr_activation(thresholds: dict, **kwargs):
    sAG, sDG = float(kwargs['DS_AG']), float(kwargs['DS_DG'])

    if ((sDG >= float(thresholds['TH_sDG'])) & (sAG > sDG)):
        return True
    else:
        return False
    

##. Orientation filters
def _filtering_Acp_orientation(**kwargs): 
    # 1-based
    posAG: int = int(kwargs['POS']) + int(kwargs['DP_AG'])
    info: dict = kwargs['ExInt_INFO']

    if info == 'Warning':
        return 0

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
    posDG: int = int(kwargs['POS']) + int(kwargs['DP_DG'])
    info: dict = kwargs['ExInt_INFO']
    if info == 'Warning':
        return 0    

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
    posAG: int = int(kwargs['POS']) + int(kwargs['DP_AG'])
    info: dict = kwargs['ExInt_INFO']
    if info == 'Warning':
        return 0
    if ((kwargs['ExInt_INFO']['eStart'] == 'unk') 
        | (kwargs['ExInt_INFO']['eEnd'] == 'unk')):
        return 0
    strand: str = info['strand']
    eStart: int = int(info['eStart'])
    eEnd: int = int(info['eEnd'])
    

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
    try:
        strand: int = kwargs['ExInt_INFO']['strand']
    except:
        return 0
    
    if ((kwargs['ExInt_INFO']['eStart'] == 'unk') 
        | (kwargs['ExInt_INFO']['eEnd'] == 'unk')):
        return 0

    posDG: int = int(kwargs['POS']) + int(kwargs['DP_DG'])
    info: dict = kwargs['ExInt_INFO']
    strand: int = kwargs['ExInt_INFO']['strand']
    eStart: int = int(info['eStart'])
    eEnd: int = int(info['eEnd'])

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
    try:
        strand: int = kwargs['ExInt_INFO']['strand']
    except:
        return 0
    
    strand = kwargs['ExInt_INFO']['strand']
    pAL, pDL = int(kwargs['DP_AL']), int(kwargs['DP_DL'])

    if ((strand == '+') & (pAL < pDL)) | ((strand == '-') & (pAL > pDL)):
        return 'LEX'
    else:
        return 'RIT'


##. Varidate variant position from close exon boundary (50 bp or 250 bp) 
def _calc_dist_from_exon(**kwargs):
    try:
        kwargs['ExInt_INFO']['eStart']
    except:
        return 0
    try:
        kwargs['ExInt_INFO']['eEnd']
    except:
        return 0
    
    if ((kwargs['ExInt_INFO']['eStart'] == 'unk') 
        | (kwargs['ExInt_INFO']['eEnd'] == 'unk')):
        return 0
    
    pos = int(kwargs['POS'])
    eStart, eEnd = int(kwargs['ExInt_INFO']['eStart']), int(kwargs['ExInt_INFO']['eEnd'])
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
            # print('Gained exon')
            return True
        else:
            # print('No gained exon')
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
    if __exits_spliceai_scores(row):
        pass
    else:
        return "Cannot predict splicing event"

    if (_varidate_var_pos_50bp(**row) == 'outside_50bp'
        and predict_gained_exon(thresholds=thresholds, **row)
        and _calc_gained_exon_size(thresholds=thresholds, **row)
        and _verify_pseudoexon_location(db_intron=db_intron, **row)):
        # print('Pseudoexon activation')
        return True
    else:
        # print('No pseudoexon activation')
        return False
    

def partial_intron_retention(row, thresholds):
    if __exits_spliceai_scores(row):
        pass
    else:
        return "Cannot predict splicing event"
    
    if ((_varidate_var_pos_250bp(**row) == 'within_250bp')
        and (_is_cryptic_Acp_activation(thresholds=thresholds, **row)) 
             or (_is_cryptic_Dnr_activation(thresholds=thresholds, **row))
        and (-251 < _bp_5prime(thresholds, **row) < 0) 
             or (0 < _bp_3prime(thresholds, **row) < 251)):
        return True
    else:
        return False


def partial_exon_deletion(row, thresholds):
    if __exits_spliceai_scores(row):
        pass
    else:
        return "Cannot predict splicing event"
    
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
    if __exits_spliceai_scores(row):
        pass
    else:
        return "Cannot predict splicing event"
    
    lost_exon_size = predict_lost_exon(thresholds=thresholds, **row)
    if ((_varidate_var_pos_50bp(**row) == 'outside_50bp')
        or (lost_exon_size is None)):
        return False
    elif ((_varidate_var_pos_50bp(**row) == 'within_50bp')
          and (lost_exon_size)):            
        info = row['ExInt_INFO']
        native_exon_length = int(info['eEnd']) - int(info['eStart']) + 1
        if lost_exon_size == native_exon_length:
            return True
        else:
            return False
    else:
        return False
    

def intron_retention(row, thresholds):
    """Varidate intron retention
    When the variant is located outside >50 bp from 
    close exon-intron boundary, intron retention may not occur.
    """
    if __exits_spliceai_scores(row):
        pass
    else:
        return "Cannot predict splicing event"
    
    if ((_varidate_var_pos_50bp(**row) == 'outside_50bp' 
        or predict_retein_intron(thresholds=thresholds, **row) is None)):
        return False
    elif ((_varidate_var_pos_50bp(**row) == 'within_50bp' 
        or predict_retein_intron(thresholds=thresholds, **row))):
        return True
    else:
        return False

## Multi-exon skipping
def multi_exon_skipping(row, thresholds: dict):
    if __exits_spliceai_scores(row):
        pass
    else:
        return "Cannot predict splicing event"
    
    if row['Exon_skipping']:
        info = row['ExInt_INFO']
        lost_exon_size = predict_lost_exon(thresholds=thresholds, **row)
        native_exon_size = np.abs(int(info['eEnd']) - int(info['eStart']) + 1)

        if lost_exon_size == native_exon_size:
            return 'One exon skipping'
        elif lost_exon_size > native_exon_size:
            print('Assumed multiple exon skipping')
            if info['strand'] == '+':
                if row['SpliceType'] == 'Donor_int':
                    try:
                        two_exons = int(info['eEnd']) - int(info['prev_ExStart']) + 1
                    except:
                        two_exons = np.nan
                elif row['SpliceType'] == 'Acceptor_int':
                    try:
                        two_exons = int(info['next_ExEnd']) - int(info['eStart']) + 1
                    except:
                        two_exons = np.nan
            else:
                if row['SpliceType'] == 'Donor_int':
                    try:
                        two_exons = int(info['prev_ExEnd']) - int(info['eStart']) + 1
                    except:
                        two_exons = np.nan
                elif row['SpliceType'] == 'Acceptor_int':
                    try:
                        two_exons = int(info['eEnd']) - int(info['next_ExStart']) + 1
                    except:
                        two_exons = np.nan
            
            if lost_exon_size == two_exons:
                return 'Double exon skipping'
            else:
                return 'unk'

################################################################################
##                   Calculate Aberrant splicing event size                   ##
################################################################################

def anno_intron_retention_size(row, thresholds):
    if __exits_spliceai_scores(row):
        pass
    else:
        return "Cannot predict splicing event"
    
    if row['Int_Retention']:
        return predict_retein_intron(thresholds=thresholds, **row)
    else:
        return np.nan

def anno_partial_intron_retention_size(row, thresholds):
    if __exits_spliceai_scores(row):
        pass
    else:
        return "Cannot predict splicing event"
    
    if row['Part_IntRet']:
        if -251 < _bp_5prime(thresholds, **row) < 0:
             return np.abs(_bp_5prime(thresholds, **row))
        elif 0 < _bp_3prime(thresholds, **row) < 251: 
            return _bp_3prime(thresholds, **row)
    else:
        return np.nan

def anno_gained_exon_size(row, thresholds):
    if __exits_spliceai_scores(row):
        pass
    else:
        return "Cannot predict splicing event"
    
    if row['Pseudoexon']:
        return _calc_gained_exon_size(thresholds=thresholds, **row)
    else:
        return np.nan

def anno_partial_exon_del_size(row, thresholds):
    if __exits_spliceai_scores(row):
        pass
    else:
        return "Cannot predict splicing event"
    
    if row['Part_ExDel']:
        if _bp_5prime(thresholds, **row) > 0:
            return _bp_5prime(thresholds, **row)
        elif _bp_3prime(thresholds, **row) < 0:
            return np.abs(_bp_3prime(thresholds, **row))
    else:
        return np.nan


def anno_skipped_exon_size(row, thresholds: dict):
    if row['Exon_skipping']:
        if row['multiexs'] == 'One exon skipping':
            return predict_lost_exon(thresholds=thresholds, **row)
        elif row['multiexs'] == 'Two exons skipping':
            info = row['ExInt_INFO']
            current_exon = int(info['eEnd']) - int(info['eStart']) + 1
            if info['strand'] == '+':
                if row['SpliceType'] == 'Donor_int':
                    try:
                        prev_exon = int(info['prev_ExEnd']) - int(info['prev_ExStart']) + 1
                    except:
                        return np.nan
                    else:
                        return current_exon + prev_exon

                elif row['SpliceType'] == 'Acceptor_int':
                    try:
                        next_exon = int(info['next_ExEnd']) - int(info['eStart']) + 1
                    except:
                        return np.nan
                    else:
                        return current_exon + next_exon
            else:
                if row['SpliceType'] == 'Donor_int':
                    try:
                        prev_exon = int(info['prev_ExEnd']) - int(info['prev_ExStart']) + 1
                    except:
                        return np.nan
                    else:
                        return current_exon + prev_exon 
                elif row['SpliceType'] == 'Acceptor_int':
                    try:
                        next_exon = int(info['next_ExEnd']) - int(info['eStart']) + 1
                    except:
                        return np.nan
                    else:
                        return current_exon + next_exon
        else:
            return np.nan
    else:
        return np.nan


# Truncated regions
def anno_skipped_regions(row):
    if __exits_spliceai_scores(row):
        pass
    else:
        return "Cannot predict splicing event"
    
    # info: dict = row['ExInt_INFO']
    try:
        strand: str = row['Strand']
    except:
        return np.nan
    posVar: int = int(row['POS'])

    if row['Exon_skipping']:
        if strand == '+':
            start = posVar + int(row['DP_AL'])
            end = posVar + int(row['DP_DL'])
        elif strand == '-':
            start = posVar + int(row['DP_DL'])
            end = posVar + int(row['DP_AL'])
        else:
            print('Warning: unkown strand')
            return np.nan
    else:
        return np.nan
    # print(f"{row['CHROM']} {str(start)} {str(end)}")
    return f"{row['CHROM']} {str(start)} {str(end)}"

def anno_deleted_regions(row, thresholds: dict):
    if __exits_spliceai_scores(row):
        pass
    else:
        return "Cannot predict splicing event"
    
    info: dict = row['ExInt_INFO']
    try:
        strand: str = row['Strand']
        eStart: int = int(info['eStart'])
        eEnd: int = int(info['eEnd'])
    except:
        return np.nan 
    
    eStart: int = int(info['eStart'])
    eEnd: int = int(info['eEnd'])
    posVar: int = int(row['POS'])

    if row['Part_ExDel']:
        if _bp_5prime(thresholds, **row) > 0:
            if strand == '+':
                start = eStart
                end = posVar + int(row['DP_AG'])
            elif strand == '-':
                start = posVar + int(row['DP_AG'])
                end = eEnd
            else:
                print('Warning: unkown strand')
                return np.nan
        elif _bp_3prime(thresholds, **row) < 0:
            if strand == '+':
                start = posVar + int(row['DP_DG'])
                end = eEnd
            elif strand == '-':
                start = eStart
                end = posVar + int(row['DP_DG'])
            else:
                print('Warning: unkown strand')
                return np.nan
        else:
            print('Warning: unkown deletion conditions')
            return np.nan
    else:
        return np.nan
    # print(f"{row['CHROM']} {str(start)} {str(end)}")
    return f"{row['CHROM']} {str(start)} {str(end)}"
    
