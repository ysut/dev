import pandas as pd

bp7_csq: set = {
    'intron_variant', 'synonymous_variant', 
    'splice_region_variant&intron_variant', '5_prime_UTR_variant',
    '3_prime_UTR_variant', 'splice_acceptor_variant&intron_variant',
    'splice_region_variant&synonymous_variant'
    }

# scores = {'clinvar_same_pos': 2,                    #1
#           'clinvar_same_motif': 1,                  #2
#           'clinvar_else': 0,                        #3 
#           'non_canon_splai_gte_0.2': 3,             #4
#           'non_canon_splai_bet_0.1_0.2': 0,         #5
#           'non_canon_splai_lte_0.1_outside': -2,    #6
#           'non_canon_splai_lte_0.1_other': -1,      #7
#           'frameshift_nmd_eloF': 7,                 #8
#           'frameshift_nmd_not_eloF': 3,             #9
#           'canon_strong': 6,                        #10
#           'canon_moderate': 5,                      #11
#           'canon_splai_lte_0.1': -2,                #12
#           'canon_splai_bet_0.1_0.2': -1,            #13
#           'canon_splai_gte_0.2': 1,                 #14
#           }

class Scoring:
    def __init__(self, ths: dict) -> None: 
        self.scores: dict = ths

    def _calc_canon_prescore(self, row) -> int:
        if row['maxsplai'] < 0.1:
            return self.scores['canon_splai_lte_0.1']
        elif row['maxsplai'] < 0.2:
            return self.scores['canon_splai_bet_0.1_0.2']
        else:
            return self.scores['canon_splai_gte_0.2']

    def insilico_screening(self, row) -> int:
        #1. Non-canonical
        if row['is_Canonical'] == 'False':
            if row['maxsplai'] >= 0.2:
                return self.scores['non_canon_splai_gte_0.2']
            elif row['maxsplai'] <= 0.1:
                if ((row['SpliceType'] == 'Acceptor_int') | (row['SpliceType'] == 'Donor_int')):
                    if ((int(row['Int_loc']) <= -21) | (int(row['Int_loc']) >= 7)):
                        raw_score = self.scores['non_canon_splai_lte_0.1_outside']
                    else:
                        raw_score = self.scores['non_canon_splai_lte_0.1_other']
                elif ((row['SpliceType'] == 'Acceptor_ex') | (row['SpliceType'] == 'Donor_ex')):
                    if row['csq'] in bp7_csq:
                        if ((int(row['ex_up_dist']) >= 1) & (int(row['ex_down_dist']) >= 3)):
                            raw_score = self.scores['non_canon_splai_lte_0.1_outside']
                        else:
                            raw_score = self.scores['non_canon_splai_lte_0.1_other']
                    else:
                        raw_score = self.scores['non_canon_splai_lte_0.1_other']
                else:
                    raw_score = self.scores['non_canon_splai_lte_0.1_other']
            else:
                raw_score = self.scores['non_canon_splai_bet_0.1_0.2']
    
        #2. Canonical
        else:
            pre_score = self._calc_canon_prescore(row)
            # Frameshift variants
            if row['is_Frameshift']:
                if row['is_NMD_at_Canon'] == 'Possibly_NMD':
                    if row['is_eLoF']:
                        raw_score = pre_score + self.scores['frameshift_nmd_eloF']
                    else:
                        raw_score = pre_score + self.scores['frameshift_nmd_not_eloF']
                else:
                    if ((float(row['skipped_ccrs']) >= 95) | (float(row['deleted_ccrs']) >=　95)):
                        raw_score = pre_score + self.scores['canon_strong']
                    else:
                        if row['is_10%_truncation']:
                            raw_score = pre_score + self.scores['canon_strong']
                        else:
                            raw_score = pre_score + self.scores['canon_moderate']
            # In-frame
            else:
                if ((float(row['skipped_ccrs']) >= 95) | (float(row['deleted_ccrs']) >= 95)):
                    raw_score = pre_score + self.scores['canon_strong']
                else:
                    if row['is_10%_truncation']:
                        raw_score = pre_score + self.scores['canon_strong']
                    else:
                        raw_score = pre_score + self.scores['canon_moderate']

            # Calibrate minus scores to 0
        if raw_score < 0:
            raw_score = 0

        return raw_score


    def clinvar_screening(self, row) -> int:
        if row['insilico_screening'] >= 0:
            if row['clinvar_same_pos']:
                return self.scores['clinvar_same_pos']
            else:
                if row['clinvar_same_motif']:
                    return self.scores['clinvar_same_motif']
                else:
                    return self.scores['clinvar_else']
        else:
            return self.scores['clinvar_else']


    def calc_priority_score(self, df: pd.DataFrame) -> pd.DataFrame:
        df['PriorityScore'] = df['insilico_screening'] + df['clinvar_screening']
        return df
