from . import pna_tools
import numpy as np
import warnings


def checkWarnings(seq):
    """

    :param seq:
    :return:
    """
    sequence_warnings = []
    if len(seq) > 30:
        sequence_warnings.append('PNA sequence too long')
    if partialComp(seq, 7):
        sequence_warnings.append('Significant self-complementary sequence (>6 bp)')

    purine_check_1 = purine_check_2 = g_check = False

    if len(seq) > 10:
        for ii in range(0, len(seq) - 9):
            sub_seq = seq[ii:ii + 10]
            if 'PPPPP' in purAndPyr(sub_seq):
                purine_check_1 = True
            if purAndPyr(sub_seq).count('P') > 6:
                purine_check_2 = True
            if 'GGGG' in purAndPyr(sub_seq):
                g_check = True
    else:
        if 'PPPPP' in purAndPyr(seq):
            purine_check_1 = True
        if purAndPyr(seq).count('P') > 6:
            purine_check_2 = True
        if 'GGGG' in purAndPyr(seq):
            g_check = True

    if purine_check_1:
        sequence_warnings.append('Purine stretch of 5 or greater')
    if purine_check_2:
        sequence_warnings.append('High purine content')
    if g_check:
        sequence_warnings.append('Guanine stretch of 4 or greater')

    return sequence_warnings


def purAndPyr(sequence):
    """

    :param sequence:
    :return:
    """
    if type(sequence) is not str:
        raise TypeError('purAndPyr takes a sequence string')

    sequence = sequence.upper()
    type_seq = ''

    for letter in sequence:
        if letter == 'A' or letter == 'G':
            type_seq += 'P'
        elif letter == 'T' or letter == 'U' or letter == 'C':
            type_seq += 'Y'
        else:
            raise ValueError('Non-RNA or non-DNA base included')

    return type_seq


def partialComp(seq, comp_len):
    """

    :param seq:
    :param comp_len:
    :return:
    """
    if type(seq) is not str:
        raise TypeError('partialComp takes a sequence string')

    seq = seq.upper()
    slide = len(seq) - comp_len + 1

    for ii in range(slide):
        subseq = seq[ii:ii + comp_len]
        rcompseq = pna_tools.revComp(subseq)
        compseq = rcompseq[::-1]

        if rcompseq in seq or compseq in seq:
            return True

    return False


def calculateTmPNA(sequence, ct=1e-4, R=1.987e-3, nacl=1, dH_dict=None, dS_dict=None, dG_dict=None):
    """

    :param sequence: N to C PNA sequence
    :param ct:
    :return:
    """
    try:
        sequence = sequence.upper()
        for letter in sequence:
            if letter not in 'ACGT':
                if letter == 'U':
                    raise ValueError('RNA base included. calculateTmPNA takes a DNA sequence string')
                else:
                    raise ValueError('calculateTmPNA takes a DNA sequence string')
    except AttributeError:
        raise TypeError('calculateTmPNA takes a string')

    if R != 1.987e-3:
        warnings.warn('Gas constant value altered, be sure to use appropriate nearest neighbor parameters')

    if not dH_dict:
        dH_dict = {'AA': -8.4, 'AT': -6.5, 'TA': -6.3, 'CA': -7.4, 'GT': -8.6, 'CT': -6.1, 'GA': -7.7, 'CG': -10.1,
                   'GC': -11.1, 'GG': -6.7}

    if not dS_dict:
        dS_dict = {'AA': -23.6, 'AT': -18.8, 'TA': -18.5, 'CA': -19.3, 'GT': -23, 'CT': -16.1, 'GA': -20.3, 'CG': -25.5,
                   'GC': -28.4, 'GG': -15.6}

    if not dG_dict:
        dG_dict = {'AA': -1.02, 'AT': -0.73, 'TA': -0.6, 'CA': -1.38, 'GT': -1.43, 'CT': -1.16, 'GA': -1.46,
                   'CG': -2.09, 'GC': -2.28, 'GG': -1.77}

    dH = dS = dG = 0

    if 'G' in sequence or 'C' in sequence:
        dS += -5.9
        dG += 1.82
    else:
        dS += -9.0
        dG += 2.8

    if sequence == pna_tools.revComp(sequence):
        dS += -1.4
        dG += 0.4
    else:
        ct *= 0.25

    if sequence[-1] == 'T':
        dH += 0.4
        dG += 0.4

    for ii in range(len(sequence) - 1):
        pair = sequence[ii:ii + 2]

        try:
            dH += dH_dict[pair]
        except KeyError:
            dH += dH_dict[pna_tools.revComp(pair)]

        try:
            dS += dS_dict[pair]
        except KeyError:
            dS += dS_dict[pna_tools.revComp(pair)]

        try:
            dG += dG_dict[pair]
        except KeyError:
            dG += dG_dict[pna_tools.revComp(pair)]

    dS *= 0.001
    temp_dna = dH / (dS + R * np.log(ct)) - 273.15
    tm_dna_salt = temp_dna + 12.5 * np.log(nacl)
    c = (20.79, 0.83, -26.13, 0.44)
    temp_pna = c[0] + (c[1] * temp_dna) + (c[2] * purAndPyr(sequence).count('Y') / len(sequence)) + (c[3] * len(sequence))

    return dH, dS, dG, temp_dna, temp_pna
