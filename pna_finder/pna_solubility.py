def checkWarnings(seq):
    warnings = []
    if len(seq) > 30:
        warnings.append('PNA sequence too long')
    if partialComp(seq, 7):
        warnings.append('Significant self-complementary sequence (>6 bp)')

    purine_check_1 = purine_check_2 = g_check = False

    if len(seq) > 10:
        for ii in range(0, len(seq) - 9):
            sub_seq = seq[ii:ii + 10]
            if 'PPPPP' in purandpyr(sub_seq):
                purine_check_1 = True
            if purandpyr(sub_seq).count('P') > 6:
                purine_check_2 = True
            if 'GGGG' in purandpyr(sub_seq):
                g_check = True
    else:
        if 'PPPPP' in purandpyr(seq):
            purine_check_1 = True
        if purandpyr(seq).count('P') > 6:
            purine_check_2 = True
        if 'GGGG' in purandpyr(seq):
            g_check = True

    if purine_check_1:
        warnings.append('Purine stretch of 5 or greater')
    if purine_check_2:
        warnings.append('High purine content')
    if g_check:
        warnings.append('Guanine stretch of 4 or greater')

    return warnings


def rev_comp(seq):
    if type(seq) is not str:
        raise TypeError('rev_comp takes a sequence string')

    seq = seq.upper()
    rc_seq = ''

    if 'T' in seq and 'U' in seq:
        raise ValueError('Mixed RNA and DNA sequences')
    elif 'U' in seq:
        seq_type = 'RNA'
    else:
        seq_type = 'DNA'  # Set default seq type as DNA

    for letter in seq[::-1]:
        if letter == 'A':
            if seq_type == 'DNA':
                rc_seq += 'T'
            elif seq_type == 'RNA':
                rc_seq += 'U'
        elif letter == 'T' or letter == 'U':
            rc_seq += 'A'
        elif letter == 'C':
            rc_seq += 'G'
        elif letter == 'G':
            rc_seq += 'C'
        else:
            raise ValueError('Non-RNA or non-DNA base included')

    return rc_seq


def purandpyr(seq):
    if type(seq) is not str:
        raise TypeError('rev_comp takes a sequence string')

    seq = seq.upper()
    type_seq = ''

    for letter in seq:
        if letter == 'A' or letter == 'G':
            type_seq += 'P'
        elif letter == 'T' or letter == 'U' or letter == 'C':
            type_seq += 'Y'
        else:
            raise ValueError('Non-RNA or non-DNA base included')

    return type_seq


def partialComp(seq, comp_len):
    if type(seq) is not str:
        raise TypeError('rev_comp takes a sequence string')

    seq = seq.upper()
    slide = len(seq) - comp_len + 1

    for ii in range(slide):
        subseq = seq[ii:ii + comp_len]
        rcompseq = rev_comp(subseq)
        compseq = rcompseq[::-1]

        if rcompseq in seq or compseq in seq:
            return True

    return False