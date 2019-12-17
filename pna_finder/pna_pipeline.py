# Applies the arguments submitted to the main PNA Finder dialog to the internal PNA Finder functions and external
# alignment programs

from . import string_genes as sg
from .align_tools import *
from . import pna_tools
from . import pna_solubility as sol
import gffutils
import os
import subprocess


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
                sol_warn    -   int (optional)
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

    filebase = out_dir + '/' + submit_dict['id_list'].split('/')[-1].split('.')[0]

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

    bedfile = filebase + '.bed'
    bedfile_edited = filebase + '_edit.bed'
    fasta_outfile = filebase + '.fa'
    sequence_outfile = filebase + '_seq.out'

    # Run pna_tools programs
    pna_tools.findID(dbfile, id_list, bedfile, feature_types=feature_types)
    pna_tools.editBed(bedfile, bedfile_edited, window=window, sequence_length=length)

    # Check version and run BEDTools
    bash_check = ['%s/bash' % bash_dir, '--login', '-c', 'bedtools --version']
    version_call = subprocess.check_output(bash_check).decode('utf-8')
    bedtools_version = version_call.split(' ')[1].split('\n')[0]

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
    sequence_obj = open(sequence_outfile, 'w')
    sequence_obj.write('Name\tSequence\t')

    if submit_dict['sol_warn'] == 1:
        sequence_obj.write('Warnings\t')

    if submit_dict['STRING'] == 1:
        sequence_obj.write('Network Genes\tNetwork Nodes\tNetwork Edges\t')

    string_names = []

    sequence_obj.write('\n')

    if submit_dict['sol_warn'] == 0 and submit_dict['STRING'] == 0:
        pass
    else:
        for key in list(sequence_dict.keys()):
            name = key.rsplit('_', 1)[0]
            sequence = sol.rev_comp(sequence_dict[key])
            sequence_obj.write('%s\t%s\t' % (key, sequence))

            if submit_dict['sol_warn'] == 1:
                warnings = sol.checkWarnings(sequence)
                for warning in warnings:
                    sequence_obj.write('%s\t' % warning)

            if submit_dict['STRING'] == 1 and name not in string_names:
                net_genes, nodes, edges = sg.STRING_net(name, species, exclude_tm=True)
                sequence_obj.write('%s\t%s\t%s' % (net_genes, nodes, edges))
                string_names.append(name)

            sequence_obj.write('\n')

    sequence_obj.close()
    return


def find_off_targets(submit_dict):
    """
    Executes the PNA Finder Find Off-Targets function using arguments provided by the dictionary submit_dict:
    :param submit_dict:
        get_sequences submit_dict keys
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
                count       -   int (optional)
                self_target -   int (optional)
    (see manual.pdf in package root directory for more info)
    :return:
    """

    bash_dir = submit_dict['bash']
    shell_type = submit_dict['shell_type']

    # Retrieve job parameters
    window = (submit_dict['start'], submit_dict['end'])
    mismatch = submit_dict['mismatch']
    pna_length = submit_dict['length']
    index_option = submit_dict['index_option']
    count = submit_dict['count']
    self_target = submit_dict['self_target']

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

    filebase = out_dir + '/' + submit_dict['targets'].split('/')[-1].split('.')[0]

    samfile = filebase + '.sam'
    bamfile = filebase + '.bam'
    sbamfile = filebase + '.sorted.bam'
    bedfile = filebase + '.bed'
    outfile = filebase + '.out'
    if count:
        countfile = filebase + '.count'
    else:
        countfile = None

    # Parse and execute index option
    if index_option == 0:
        # Check version
        bash_check = ['%s/bash' % bash_dir, '--login', '-c', 'bowtie2-build --version']
        version_call = subprocess.check_output(bash_check).decode('utf-8')
        bt2_build_version = version_call.split('\n')[0].split()[-1]

        index_base = assembly.rsplit('.', 1)[0]
        bowtie2_build(assembly, index_base,
                      bash_sub=shell_type, bash_dir=bash_dir,
                      version=bt2_build_version,
                      err_file='%s/error_files/bowtie2-build.stderr' % out_dir)
    elif index_option == 1:
        index_base = assembly.rsplit('.', 2)[0]
        if index_base.split('.')[-1] == 'rev':
            index_base = index_base.rsplit('.', 1)[0]
    else:
        raise ValueError('Something has gone horribly wrong!')

    # Check version and run Bowtie2
    bash_check = ['%s/bash' % bash_dir, '--login', '-c', 'bowtie2 --version']
    version_call = subprocess.check_output(bash_check).decode('utf-8')
    bt2_version = version_call.split('\n')[0].split()[-1]

    flags = ['-f', '-N', str(mismatch), '-L', str(pna_length), '-a']
    if mismatch == 0:
        flags += ['--no-1mm-upfront']

    bowtie2(targets,
            index_base,
            samfile,
            flags,
            bash_sub=shell_type, bash_dir=bash_dir,
            err_file='%s/error_files/bowtie2.stderr' % out_dir,
            version=bt2_version)

    # Check version and run SAMTools
    bash_check = ['%s/bash' % bash_dir, '--login', '-c', 'samtools --version']
    version_call = subprocess.check_output(bash_check).decode('utf-8')
    samtools_version = version_call.split(' ')[1].split('\n')[0]

    samtools('view', samfile, outfile=bamfile, flags=['-b'], bash_sub=shell_type, bash_dir=bash_dir,
             err_file='%s/error_files/samtools_view.stderr' % out_dir,
             version=samtools_version)
    samtools('sort', bamfile, outfile=sbamfile, bash_sub=shell_type, bash_dir=bash_dir,
             err_file='%s/error_files/samtools_sort.stderr' % out_dir,
             version=samtools_version)
    samtools('index', sbamfile, bash_sub=shell_type, bash_dir=bash_dir,
             err_file='%s/error_files/samtools_index.stderr' % out_dir,
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
             outfile=bedfile,
             flags=flags,
             bash_sub=shell_type, bash_dir=bash_dir,
             err_file='%s/error_files/bedtools.stderr' % out_dir,
             version=bedtools_version)

    # Filter alignments
    pna_tools.process_bedwindow(bedfile,
                                outfile,
                                dist_filter=window[1],
                                OT_count_file=countfile,
                                self_OT=self_target,
                                feature_types=feature_types)

    return
