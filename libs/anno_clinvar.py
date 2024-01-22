import pysam

def anno_same_pos_vars(row, tabixfile: pysam.pysam.libctabix.TabixFile):
    samepos = []
    query_chr: str = f"{row['CHROM']}"
    query_pos: int = row['POS']
    query_ref: int = row['REF']
    query_alt: int = row['ALT']
    query_variant: str = f"{query_chr}-{query_pos}-{query_ref}-{query_alt}"
    query_start: int = int(query_pos) - 1
    query_end: int = int(query_pos) + 1
    clinvars = tabixfile.fetch(
        query_chr, query_start, query_end, parser=pysam.asBed())

    # main loop 
    while 1:
        try:
            clinvar = next(clinvars)
        except StopIteration:
            break
        else:
            clinvar_variant = f"{clinvar.contig}-{clinvar.start}-{clinvar[3]}-{clinvar[4]}"
            if query_variant == clinvar_variant:
                samepos.append(f'{clinvar[7]}_{clinvar[0]}:{clinvar[1]}')
            else:
                pass

    return samepos

def anno_same_motif_vars(row, tabixfile: pysam.pysam.libctabix.TabixFile):
    """ Fetch the variants within the same splicing motif from ClinVar data

    Args:
        row (_type_): _description_

    Returns:
        list: Return a list of variants registered in ClinVar 
                within the same splicing motif (Prediction_CHROM:POS)
    """
    
    # Initialize result list & set common variables
    samemotifs = []
    query_chr: str = f'{row["CHROM"]}'
    query_pos: int = int(row['POS'])
    
    # Generate query positions 
    if (row['SpliceType'] == 'Donor_int' 
        or row['SpliceType'] == 'Donor_ex'):

        if row['SpliceType'] == 'Donor_int':
            exonside: int = - int(row['Int_loc']) - 3 
            intronside: int = - int(row['Int_loc']) + 7
        else: # row['SpliceType'] == 'Donor_ex' 
            exonside: int = - int(row['exon_pos']) - 3
            intronside: int = - int(row['exon_pos']) + 7
        
        if row['Strand'] == '+':
            query_start: int = query_pos + exonside
            query_end: int = query_pos + intronside
        elif row['Strand'] == '-':
            query_start: int = query_pos - intronside
            query_end: int = query_pos - exonside
        else:
            return 'unk_Strand'
    
    elif (row['SpliceType'] == 'Acceptor_int' 
          or row['SpliceType'] == 'Acceptor_ex'):
        
        if row['SpliceType'] == 'Acceptor_int':
            exonside: int = int(row['Int_loc']) - 1
            intronside: int = int(row['Int_loc']) + 21
        else: # row['SpliceType'] == 'Acceptor_ex'
            exonside: int = int(row['exon_pos']) - 1 
            intronside: int = int(row['exon_pos']) + 21 
            
        if row['Strand'] == '+':
            query_start: int = query_pos - intronside
            query_end: int = query_pos - exonside
        elif row['Strand'] == '-':
            query_start: int = query_pos + exonside
            query_end: int = query_pos + intronside
        else:
            return 'unk_Strand'
        
    else:
        return 'unk_SpliceType'

    clinvars = tabixfile.fetch(
        query_chr, query_start, query_end, parser=pysam.asBed())
        
    # Main loop
    while 1:
        try:
            clinvar = next(clinvars)
        except StopIteration:
            break
        else:
            registered_var = f'{clinvar[7]}_{clinvar[0]}:{clinvar[1]}'
            samemotifs.append(registered_var)

    return samemotifs