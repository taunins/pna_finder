# Applies the arguments submitted to the main PNA Finder dialog to the internal PNA Finder functions and external
# alignment programs

from . import string_genes as sg
from .align_tools import *
from . import pna_tools
from . import pna_solubility as sol
import gffutils
import os
import subprocess
import shutil


def get_sequences(submit_dict):
    """
    Executes the PNA Finder Get Sequences function using arguments provided by the dictionary submit_dict:
    :param submit_dict:
        get_sequences submit_dict keys
            startup keys
                bash        -   str
                shell_type  -   str
            job parameter keys
                start       -   int
                end         -   int
                pna_length  -   int
                feature_type -  str
            file keys
                gff_option  -   int
                id_list     -   str
                assembly    -   str
                gff         -   str
                output      -   str
            option keys
                STRING      -   int (optional)
                string_id   -   int (optional)
                warnings    -   int (optional)
    (see manual.pdf in package root directory for more info)
    :return:
    """

    bash_dir = submit_dict['bash']
    shell_type = submit_dict['shell_type']

    # Retrieve job parameters
    window = (submit_dict['start'], submit_dict['end'])
    length = submit_dict['length']

    feature_types = submit_dict['feature_type'].split(',')
    for ii in range(len(feature_types)):
        feature_types[ii] = feature_types[ii].strip()

    gff_option = submit_dict['gff_option']
    warnings_option = submit_dict['warnings']
    temperature_option = submit_dict['temp_option']
    string_option = submit_dict['STRING']

    species = None
    if submit_dict['STRING'] == 1:
        try:
            species = submit_dict['string_id']
        except KeyError:
            raise KeyError('Species Taxonomy ID must be specified')

    # Retrieve job file paths
    id_list = submit_dict['id_list']
    assembly = submit_dict['assembly']
    annotation = submit_dict['gff']

    # Initialize folder and files
    out_dir = submit_dict['output']
    out_base = out_dir
    out_num = 0
    folder_check = True
    name_change = False

    while folder_check:
        try:
            os.mkdir(out_dir)
            folder_check = False
            os.mkdir('%s/error_files' % out_dir)
        except WindowsError:
            name_change = True
            out_dir = out_base + ' (%s)' % str(out_num + 1)
            out_num += 1

    if name_change:
        print('Directory %s already exists, job results written to directory %s' % (out_base, out_dir))

    filebase = submit_dict['id_list'].split('/')[-1].split('.')[0]

    #   Parse and execute gff option
    if gff_option == 0:
        dbfile = annotation.rsplit('.', 1)[0] + '.db'
        gffutils.create_db(annotation, dbfn=dbfile,
                           force=True, keep_order=True,
                           merge_strategy='merge', sort_attribute_values=True)
    elif gff_option == 1:
        dbfile = annotation
    else:
        raise ValueError('Something has gone horribly wrong!')

    bedfile = out_dir + '/' + filebase + '.bed'
    bedfile_edited = out_dir + '/' + filebase + '.edit.bed'
    fasta_outfile = out_dir + '/' + filebase + '.fa'
    sequence_outfile = out_dir + '/' + filebase + '.out'
    string_outfile = None
    if string_option:
        string_outfile = out_dir + '/' + 'string.out'

    # Check version for BEDTools
    bash_check = ['%s/bash' % bash_dir, '--login', '-c', 'bedtools --version']
    version_call = subprocess.check_output(bash_check).decode('utf-8')
    bedtools_version = version_call.split(' ')[1].split('\n')[0]

    # Run pna_tools programs
    pna_tools.findID(dbfile, id_list, bedfile, feature_types=feature_types, id_type='id')
    pna_tools.editBed(bedfile, bedfile_edited, window=window, sequence_length=length)

    # Run BEDTools getfasta
    flags = ['-name', '-s']
    bedtools('getfasta',
             [assembly, bedfile_edited],
             outfile=fasta_outfile,
             flags=flags,
             bash_sub=shell_type, bash_dir=bash_dir,
             err_file='%s/error_files/bedtools.stderr' % out_dir,
             version=bedtools_version)

    # Run sequence solubility checks and STRING interaction finder
    sequence_dict = pna_tools.fastaToDict(fasta_outfile)
    sequence_handle = open(sequence_outfile, 'w')
    sequence_handle.write('Name\tSequence')

    if warnings_option == 1:
        sequence_handle.write('\tWarnings')
    if temperature_option == 1:
        sequence_handle.write('\tTm (DNA/ssPNA)')
    sequence_handle.write('\n')

    string_handle = None
    if submit_dict['STRING'] == 1:
        string_handle = open(string_outfile, 'w')
        string_handle.write('Gene\tNetwork Nodes\tNetwork Edges\tNetwork Genes\n')
    string_names = []

    if not warnings_option and not temperature_option and not string_option:
        pass
    else:
        for key in list(sequence_dict.keys()):
            name = key.rsplit('_', 1)[0]
            sequence = pna_tools.revComp(sequence_dict[key])
            sequence_handle.write('%s\t%s\t' % (key, sequence))

            if warnings_option:
                warnings = sol.checkWarnings(sequence)
                if warnings:
                    first = True
                    for warning in warnings:
                        if first:
                            sequence_handle.write('%s' % warning)
                            first = False
                        else:
                            sequence_handle.write(', %s' % warning)

            if temperature_option:
                dH, dS, dG, temp_dna, temp_pna = sol.calculateTmPNA(sequence)
                sequence_handle.write('\t%s' % round(temp_pna, 1))

            sequence_handle.write('\n')

            if string_option and name not in string_names:
                net_genes, nodes, edges = sg.STRING_net(name, species, exclude_tm=True)
                string_handle.write('%s\t%s\t%s\t' % (name, nodes, edges))
                first = True
                for gene in net_genes:
                    if first:
                        string_handle.write('%s' % gene)
                        first = False
                    string_handle.write(',%s' % gene)
                string_handle.write('\n')
                string_names.append(name)

    sequence_handle.close()
    try:
        string_handle.close()
    except AttributeError:
        pass

    return


def find_off_targets(submit_dict):
    """
    Executes the PNA Finder Find Off-Targets function using arguments provided by the dictionary submit_dict:
    :param submit_dict:
        find_off_targets submit_dict keys
            startup keys
                bash        -   str
                shell_type  -   str
            job parameter keys
                start       -   int
                end         -   int
                mismatch    -   int
                length      -   int
                feature_type -  str
            file keys
                index_option -  int
                targets     -   str
                assembly    -   str
                gff         -   str
                output      -   str
            option keys
                log_mismatch -  int (optional)
                count       -   int (optional)
                homology    -   int (optional)
                remove_files -  int (optional)
    (see manual.pdf in package root directory for more info)
    :return:
    """

    bash_dir = submit_dict['bash']
    shell_type = submit_dict['shell_type']

    # Retrieve job parameters
    window = (submit_dict['start'], submit_dict['end'])
    mismatches = submit_dict['mismatch']
    index_option = submit_dict['index_option']
    count = submit_dict['count']
    check_homology = submit_dict['homology']
    remove_files = submit_dict['remove_files']
    small_genome = False  # TODO

    feature_types = submit_dict['feature_type'].split(',')
    for ii in range(len(feature_types)):
        feature_types[ii] = feature_types[ii].strip()

    # Retrieve job file paths
    targets = submit_dict['targets']
    assembly = submit_dict['assembly']
    annotation = submit_dict['gff']

    # Initialize folder and files
    out_dir = submit_dict['output']
    out_base = out_dir
    out_num = 0
    folder_check = True
    name_change = False

    while folder_check:
        try:
            os.mkdir(out_dir)
            folder_check = False
            os.mkdir('%s/error_files' % out_dir)
        except WindowsError:
            name_change = True
            out_dir = out_base + ' (%s)' % str(out_num + 1)
            out_num += 1

    if name_change:
        print('Directory %s already exists, job results written to directory %s' % (out_base, out_dir))

    filebase = submit_dict['targets'].split('/')[-1].split('.')[0]
    os.mkdir(out_dir + '/sambam/')

    # Parse and execute index option
    if index_option == 0:
        # Check version
        bash_check = ['%s/bash' % bash_dir, '--login', '-c', 'bowtie-build --version']
        version_call = subprocess.check_output(bash_check).decode('utf-8')
        bt_build_version = version_call.split('\n')[0].split()[-1]

        index_base = assembly.rsplit('.', 1)[0]
        bowtie_build(assembly, index_base,
                     bash_sub=shell_type, bash_dir=bash_dir,
                     version=bt_build_version,
                     err_file='%s/error_files/bowtie-build.stderr' % out_dir)
    elif index_option == 1:
        index_base = assembly.rsplit('.', 2)[0]
        if index_base.split('.')[-1] == 'rev':
            index_base = index_base.rsplit('.', 1)[0]
    else:
        raise ValueError('Something has gone horribly wrong!')

    # Check version and run Bowtie
    bash_check = ['%s/bash' % bash_dir, '--login', '-c', 'bowtie --version']
    version_call = subprocess.check_output(bash_check).decode('utf-8')
    bt_version = version_call.split('\n')[0].split()[-1]

    # Find PNA length
    targets_dict = pna_tools.fastaToDict(targets)
    length_dict = {}
    sam_dict = {}

    for key in targets_dict:
        length = str(len(targets_dict[key]))
        try:
            length_dict[length][key] = targets_dict[key]
            if '%s.rev' % key not in length_dict:  # Add reverse sequence for N to 5' parallel alignment
                length_dict[length]['%s.rev' % key] = targets_dict[key][::-1]
        except KeyError:
            length_dict[length] = {}
            length_dict[length][key] = targets_dict[key]
            if '%s.rev' % key not in length_dict:  # Add reverse sequence for N to 5' parallel alignment
                length_dict[length]['%s.rev' % key] = targets_dict[key][::-1]

    # Run Bowtie for each length of PNA, each mismatch tolerance
    for key in length_dict:
        if len(length_dict) > 1:
            print('Separating FASTA file for PNA length %s...' % key)

        temp_targets = out_dir + '/' + filebase + '.%smer.fa' % key
        pna_tools.dictToFasta(length_dict[key], temp_targets)
        temp_len = key

        flags = ['-f', '-l', str(key), '-a']
        if small_genome:
            flags += ['-y']

        for mismatch in mismatches:
            flags += ['-n', str(mismatch)]
            sambase = out_dir + '/sambam/%sMM/' % mismatch + filebase
            os.mkdir('%s/sambam/%sMM/' % (out_dir, mismatch))
            os.mkdir('%s/error_files/%sMM/' % (out_dir, mismatch))

            sam_dict[mismatch] = []
            temp_sam = sambase + '.%smer.sam' % key
            temp_bam = sambase + '.%smer.bam' % key
            temp_sbam = sambase + '.%smer.sorted.bam' % key
            temp_bed = sambase + '.%smer.bed' % key
            out_bed = sambase + '.bed'
            sam_dict[mismatch] += [[temp_len, temp_sam, temp_bam, temp_sbam, temp_bed, out_bed]]

            bowtie(temp_targets,
                   index_base,
                   temp_sam,
                   flags,
                   bash_sub=shell_type, bash_dir=bash_dir,
                   err_file='%s/error_files/%sMM/bowtie.%smer.stderr' % (out_dir, mismatch, key),
                   version=bt_version)

    # Check SAMTools version
    bash_check = ['%s/bash' % bash_dir, '--login', '-c', 'samtools --version']
    version_call = subprocess.check_output(bash_check).decode('utf-8')
    samtools_version = version_call.split(' ')[1].split('\n')[0]

    bed_handle = None
    bed_dict = {}

    for mismatch in sam_dict:
        first = True
        for sam in sam_dict[mismatch]:
            length = sam[0]
            samfile = sam[1]
            bamfile = sam[2]
            sbamfile = sam[3]
            temp_bed = sam[4]
            out_bed = sam[5]

            samtools('view', samfile, outfile=bamfile, flags=['-b'], bash_sub=shell_type, bash_dir=bash_dir,
                     err_file='%s/error_files/%sMM/samtools_view.%smer.stderr' % (out_dir, mismatch, length),
                     version=samtools_version)
            samtools('sort', bamfile, outfile=sbamfile, bash_sub=shell_type, bash_dir=bash_dir,
                     err_file='%s/error_files/%sMM/samtools_sort.%smer.stderr' % (out_dir, mismatch, length),
                     version=samtools_version)
            samtools('index', sbamfile, bash_sub=shell_type, bash_dir=bash_dir,
                     err_file='%s/error_files/%sMM/samtools_index.%smer.stderr' % (out_dir, mismatch, length),
                     version=samtools_version)

            # Check version and run BEDTools
            bash_check = ['%s/bash' % bash_dir, '--login', '-c', 'bedtools --version']
            version_call = subprocess.check_output(bash_check).decode('utf-8')
            bedtools_version = version_call.split(' ')[1].split('\n')[0]

            bedtools_function = 'window'
            infiles = [sbamfile, annotation]
            flags = ['-l', str(-1 * window[0]), '-r', str(-1 * window[0]), '-sw', '-sm', '-bed']

            bedtools(bedtools_function,
                     infiles,
                     outfile=temp_bed,
                     flags=flags,
                     bash_sub=shell_type, bash_dir=bash_dir,
                     err_file='%s/error_files/bedtools.%smer.stderr' % (out_dir, length),
                     version=bedtools_version)

            if len(sam_dict[mismatch]) > 1:
                if first:
                    bed_handle = open(out_bed, 'w')
                    first = False
                    bed_dict[mismatch] = out_bed
                with open(temp_bed) as temp_bed_handle:
                    for line in temp_bed_handle:
                        bed_handle.write(line)
                    temp_bed_handle.close()
                os.remove(temp_bed)
            else:
                os.rename(temp_bed, out_bed)
                bed_dict[mismatch] = out_bed

        try:
            bed_handle.close()
        except AttributeError:
            pass

    # Filter alignments
    first = True
    for mismatch in bed_dict:
        bedfile = bed_dict[mismatch]
        outfile = out_dir + '/' + filebase + '.%smm.out' % mismatch
        if count:
            countfile = out_dir + '/' + filebase + '.%smm.count' % mismatch
        else:
            countfile = None

        if check_homology and first:
            homology_outfile = out_dir + '/' + 'homology.out'
            first = False
        else:
            check_homology = False
            homology_outfile = None

        pna_tools.processBedWindow(bedfile,
                                   outfile,
                                   dist_filter=window[1],
                                   ot_count_file=countfile,
                                   homology_outfile=homology_outfile,
                                   check_homology=check_homology,
                                   feature_types=feature_types)

    if remove_files:
        shutil.rmtree('%s/error_files' % out_dir)
        shutil.rmtree('%s/sambam' % out_dir)

    return


def sequence_warnings(submit_dict):
    """
    Executes the PNA Finder Check Sequence Warnings function using arguments provided by the dictionary submit_dict:
    :param submit_dict:
       check_warnings submit_dict keys
            file keys
                pna_option  -   int
                pna_list    -   str
            option keys
                ss_tm       -  int (optional)
    (see manual.pdf in package root directory for more info)
    :return:
    """

    if submit_dict['pna_option'] == 0:
        pna_targets = submit_dict['pna_list']
        sequence_dict = pna_tools.fastaToDict(pna_targets)
    elif submit_dict['pna_option'] == 1:
        pna_list = submit_dict['pna_list']
        sequence_dict = dict()
        with open(pna_list) as pna_handle:
            for line in pna_handle:
                [key, value] = line.rstrip().split('\t')
                sequence_dict[key] = pna_tools.revComp(value)
    else:
        raise ValueError('Something has gone horribly wrong!')

    temperature_option = submit_dict['temp_option']

    # Initialize folder and files
    out_dir = submit_dict['output']
    out_base = out_dir
    out_num = 0
    folder_check = True
    name_change = False

    while folder_check:
        try:
            os.mkdir(out_dir)
            folder_check = False
            os.mkdir('%s/error_files' % out_dir)
        except WindowsError:
            name_change = True
            out_dir = out_base + ' (%s)' % str(out_num + 1)
            out_num += 1

    if name_change:
        print('Directory %s already exists, job results written to directory %s' % (out_base, out_dir))

    filebase = submit_dict['pna_list'].split('/')[-1].split('.')[0]

    sequence_outfile = out_dir + '/' + filebase + '.out'
    sequence_handle = open(sequence_outfile, 'w')
    sequence_handle.write('Name\tSequence\tWarnings')

    if temperature_option:
        sequence_handle.write('\tTm (DNA/ssPNA)')

    sequence_handle.write('\n')

    for key in list(sequence_dict.keys()):
        sequence = pna_tools.revComp(sequence_dict[key])
        sequence_handle.write('%s\t%s\t' % (key, sequence))

        warnings = sol.checkWarnings(sequence)
        if warnings:
            first = True
            for warning in warnings:
                if first:
                    sequence_handle.write('%s' % warning)
                    first = False
                else:
                    sequence_handle.write(', %s' % warning)

        if temperature_option:
            dH, dS, dG, temp_dna, temp_pna = sol.calculateTmPNA(sequence)
            sequence_handle.write('\t%s' % round(temp_pna, 1))

        sequence_handle.write('\n')

    sequence_handle.close()
    return
