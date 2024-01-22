import pnadas as pd

def apply_bp4(row):
    if row['maxsplai'] <= 0.1:
        return 'BP4'
    else:
        return '.'

def apply_bp7(row):
    if row['BP4'] == 'BP4':
        if (row['exon_splice_site'] == 'exon_splice_site' 
            or int(row['Int_loc']) >= 7 
            or int(row['Int_loc']) <= -21):
            return 'BP7'
        else:
            return '.'
    return '.'

def apply_pp3(row):
    if row['maxsplai'] >= 0.2:
        return 'PP4'
    else:
        return '.'

def apply_ps1(row):
    if row['clinvar_same_pos']:
        return 'PS1'
    else:
        return '.'

def apply_ps1_moderate(row):
    if row['clinvar_same_motif']:
        return 'PS1_Moderate'
    else:
        return '.'


def apply_pvs1(row):
    return