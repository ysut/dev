

def apply_PP3(row):
    if row['is_Canonical'] == 'False':
        if row['maxsplai'] >= 0.2:
            return 'PP3'

def apply_BP4(row):
    if row['is_Canonical'] == 'False':
        if row['maxsplai'] <= 0.1:
            return 'BP4'

bp7_csq: set = {
    'intron_variant', 'synonymous_variant', 
    'splice_region_variant&intron_variant', '5_prime_UTR_variant',
    '3_prime_UTR_variant', 'splice_acceptor_variant&intron_variant',
    'splice_region_variant&synonymous_variant'
    }

def apply_BP7(row):
    if row['BP4'] == 'BP4':
        if ((row['SpliceType'] == 'Acceptor_int') 
            | (row['SpliceType'] == 'Donor_int')):
            if ((int(row['Int_loc']) <= -21) | (int(row['Int_loc']) >= 7)):
                return 'BP7'
            
        elif ((row['SpliceType'] == 'Acceptor_ex') 
              | (row['SpliceType'] == 'Donor_ex')):
            if row['csq'] in bp7_csq:
                if ((int(row['ex_up_dist']) >= 1) & (int(row['ex_down_dist']) >= 3 )):
                    return 'BP7'
            else:
                return None
    else:
        return None

def apply_PVS1s(row):
    def _elof_judge(row) -> str:
        if row['is_eLoF']:
            return 'PVS1'
        else:
            return None

    def _truncated_len_judge(row) -> str:
        # print(row['is_10%_truncation'])
        if row['is_10%_truncation']:
            return 'PVS1_Strong'
        else:
            return 'PVS1_Moderate'

    def _ccrs_judge(row) -> str:
        if ((float(row['skipped_ccrs']) >= 0.95) 
            |(float(row['deleted_ccrs']) >= 0.95)):
            return 'PVS1_Strong'
        else:
            return _truncated_len_judge(row)

    if row['is_Canonical'] == 'True':
        if row['is_Frameshift']:
            if row['is_NMD_at_Canon'] == 'Possibly_NMD':
                return _elof_judge(row)
            else:
                return _ccrs_judge(row)
        else:
            return _ccrs_judge(row)

def apply_PS1s(row):
    if row['clinvar_same_pos']:
        return 'PS1'
    else:
        if row['clinvar_same_motif']:
            return 'PS1_Moderate'
        else:
            return None


def final_evaluation(row):
    if row['PVS1s'] == 'PVS1':
        if ((row['PS1s'] == 'PS1') | (row['PS1s'] == 'PS1_Moderate')):
            return 'Very high priority'
        else:
            return 'High priority'
    elif row['PVS1s'] == 'PVS1_Strong':
        if ((row['PS1s'] == 'PS1') | (row['PS1s'] == 'PS1_Moderate')):
            return 'High priority'
        else:
            return 'Moderate priority'
    elif row['PVS1s'] == 'PVS1_Moderate':
        if ((row['PS1s'] == 'PS1') | (row['PS1s'] == 'PS1_Moderate')):
            return 'Moderate priority'
        else:
            return 'Low priority'
    elif row['PP3'] == 'PP3':
        if ((row['PS1s'] == 'PS1') | (row['PS1s'] == 'PS1_Moderate')):
            return 'High priority'
        else:
            return 'Low priority'
    elif row['BP4'] == 'BP4':
        if row['BP7'] != 'BP7':
            return 'Low priority'
        elif row['BP7'] == 'BP7':
            if row['PS1s'] == 'PS1_Moderate':
                return 'Low priority'
            else:
                return 'Very low priority'
    else:
        return 'Low priority'
        
    
    

