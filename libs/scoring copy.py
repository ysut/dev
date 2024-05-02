import pandas as pd

bp7_csq: set = {
    'intron_variant', 'synonymous_variant', 
    'splice_region_variant&intron_variant', '5_prime_UTR_variant',
    '3_prime_UTR_variant', 'splice_acceptor_variant&intron_variant',
    'splice_region_variant&synonymous_variant'
    }

def _calc_canon_prescore(row) -> int:
    if row['maxsplai'] < 0.1:
        return -2
    elif row['maxsplai'] < 0.2:
        return -1
    else:
        return 0

def insilico_screening(row) -> int:
    #1. Non-canonical
    if row['is_Canonical'] == 'False':
        if row['maxsplai'] >= 0.2:
            return 3
        elif row['maxsplai'] <= 0.1:
            if ((row['SpliceType'] == 'Acceptor_int') | (row['SpliceType'] == 'Donor_int')):
                if ((int(row['Int_loc']) <= -21) | (int(row['Int_loc']) >= 7)):
                    return 0
                else:
                    return 1
            elif ((row['SpliceType'] == 'Acceptor_ex') | (row['SpliceType'] == 'Donor_ex')):
                if row['csq'] in bp7_csq:
                    if ((int(row['ex_up_dist']) >= 1) & (int(row['ex_down_dist']) >= 3)):
                        return 0
                    else:
                        return 1
                else:
                    return 1
            else:
                return 1
        else:
            return 2
  
    #2. Canonical
    else:
        pre_score = _calc_canon_prescore(row)
        # Frameshift variants
        if row['is_Frameshift']:
            if row['is_NMD_at_Canon'] == 'Possibly_NMD':
                if row['is_eLoF']:
                    return pre_score + 7
                else:
                    return pre_score + 4
            else:
                if ((float(row['skipped_ccrs']) >= 0.95) | (float(row['deleted_ccrs']) >= 0.95)):
                    return pre_score + 6
                else:
                    if row['is_10%_truncation']:
                        return pre_score + 6
                    else:
                        return pre_score + 5
        # In-frame
        else:
            if ((float(row['skipped_ccrs']) >= 0.95) | (float(row['deleted_ccrs']) >= 0.95)):
                return pre_score + 6
            else:
                if row['is_10%_truncation']:
                    return pre_score + 6
                else:
                    return pre_score + 5


def clinvar_screening(row) -> int:
    if row['clinvar_same_motif'] == 'unk_SpliceType':
        return 0
    else:        
        if row['clinvar_same_pos']:
            return 2
        else:
            if row['clinvar_same_motif']:
                return 1
            else:
                return 0


def calc_priority_score(df: pd.DataFrame) -> pd.DataFrame:
    df['PriorityScore'] = df['insilico_screening'] + df['clinvar_screening']
    return df
