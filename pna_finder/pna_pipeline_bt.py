# These functions apply the arguments submitted to the main PNA Finder dialog to the internal PNA Finder functions
# (found within pna_tools.py and pna_solubility.py) and external alignment programs (Bowtie, SAMtools, bedtools)

# Import the necessary packages
from . import string_genes as sg
from Bio import SeqIO, Seq
from .align_tools import *
from . import pna_tools
from . import pna_solubility as sol
import gffutils
import os
import subprocess
import shutil
import warnings


def getSequences(submit_dict):
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
                string_id   -   int (optional)
            file keys
                id_list     -   str
                assembly    -   str
                gff         -   str
                output      -   str
            option keys
                full_search -   int (optional)
                gff_option  -   int
                STRING      -   int (optional)
                warnings    -   int (optional)
                temp_option -   int (optional)
    (see manual.pdf in package root directory for more info)
    :return:
    """

    # Retrieve bash shell location and shell type
    bash_dir = submit_dict['bash']
    shell_type = submit_dict['shell_type']

    # Retrieve job parameters
    window = (submit_dict['start'], submit_dict['end'])
    length = tuple(submit_dict['length'])

    feature_types = submit_dict['feature_type'].split(',')
    for ii in range(len(feature_types)):
        feature_types[ii] = feature_types[ii].strip()
    feature_types = tuple(feature_types)

    # Retrieve option selections
    full_search = submit_dict['full_search']
    gff_option = submit_dict['gff_option']
    warnings_option = submit_dict['warnings']
    temperature_option = submit_dict['temp_option']
    string_option = submit_dict['STRING']

    # Retrieve STRING analysis parameters
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

    # Initialize output folder
    out_dir = submit_dict['output']
    out_base = out_dir
    out_num = 0
    folder_check = True
    name_change = False

    # Create output folder, append number if already exists
    while folder_check:
        try:
            os.mkdir(out_dir)
            folder_check = False
            os.mkdir('%s/error_files' % out_dir)
        except WindowsError:
            name_change = True
            out_dir = out_base + ' (%s)' % str(out_num + 1)
            out_num += 1
    os.mkdir(out_dir + '/intermediate_files/')

    # Notify user of output folder name change
    if name_change:
        print('Directory %s already exists, job results written to directory %s' % (out_base, out_dir))

    # Parse and execute gff option
    if gff_option == 0:
        dbfile = annotation.rsplit('.', 1)[0] + '.db'
        gffutils.create_db(annotation, dbfn=dbfile, id_spec=["gene", "Name"],
                           force=True, keep_order=True,
                           merge_strategy='merge', sort_attribute_values=True)
    elif gff_option == 1:
        dbfile = annotation
    else:
        raise ValueError('Something has gone horribly wrong!')

    # Initialize output file names
    filebase = submit_dict['id_list'].split('/')[-1].split('.')[0]
    bedfile = out_dir + '/intermediate_files/' + filebase + '.bed'
    bedfile_edited = out_dir + '/intermediate_files/' + filebase + '.edit.bed'
    fasta_outfile = out_dir + '/' + filebase + '.fa'
    sequence_outfile = out_dir + '/' + filebase + '.out'
    string_outfile = None
    if string_option:
        string_outfile = out_dir + '/' + 'string.out'

    # Run pna_tools programs
    pna_tools.findID(dbfile, id_list, bedfile, feature_types=feature_types, id_type='id', full_search=full_search,
                     error_file='%s/error_files/findID.stderr' % out_dir)
    pna_tools.editBed(bedfile, bedfile_edited, window=window, sequence_length=length)

    # Check version for BEDTools
    bash_check = ['%s/bash' % bash_dir, '--login', '-c', 'bedtools --version']
    version_call = subprocess.check_output(bash_check).decode('utf-8')
    bedtools_version = version_call.split(' ')[1].split('\n')[0]

    # Run BEDTools getfasta to construct output FASTA file
    flags = ['-name', '-s']
    bedtools('getfasta',
             [assembly, bedfile_edited],
             outfile=fasta_outfile,
             flags=flags,
             bash_sub=shell_type, bash_dir=bash_dir,
             err_file='%s/error_files/bedtools.stderr' % out_dir,
             version=bedtools_version)

    # Run sequence solubility checks and STRING interaction finder
    sequence_dict = SeqIO.to_dict(SeqIO.parse(fasta_outfile, "fasta"))
    sequence_handle = open(sequence_outfile, 'w')
    sequence_handle.write('Name\tSequence')

    # Add headings for analysis options
    if warnings_option == 1:
        sequence_handle.write('\tWarnings')
    if temperature_option == 1:
        sequence_handle.write('\tTm (DNA/ssPNA)')
    sequence_handle.write('\n')

    # Initialize STRING file if option selected
    string_handle = None
    if submit_dict['STRING'] == 1:
        string_handle = open(string_outfile, 'w')
        string_handle.write('Gene\tNetwork Nodes\tNetwork Edges\tNetwork Genes\n')
    string_names = []

    # Run analysis and construct files for selected options
    if not warnings_option and not temperature_option and not string_option:
        pass
    else:
        for pna in list(sequence_dict.keys()):
            name = pna.rsplit('_', 1)[0]
            sequence = str(sequence_dict[pna].seq.reverse_complement())
            sequence_handle.write('%s\t%s\t' % (pna, sequence))

            # Check for solubility, complement warnings
            if warnings_option:
                warns = sol.checkWarnings(sequence)
                if warns:
                    sequence_handle.write(', '.join(warns))

            # Check melting temperature of PNA/ssDNA duplex
            if temperature_option:
                dH, dS, dG, temp_dna, temp_pna = sol.calculateTmPNA(sequence)
                sequence_handle.write('\t%s' % round(temp_pna, 1))

            sequence_handle.write('\n')

            # Run STRING protein interaction network analysis
            if string_option and name not in string_names:
                net_genes, nodes, edges = sg.STRING_net(name, species, exclude_tm=True)
                string_handle.write('%s\t%s\t%s' % (name, nodes, edges))
                try:
                    string_handle.write('\t%s' % ','.join(net_genes))
                except TypeError:
                    pass
                string_handle.write('\n')
                string_names.append(name)

    # Close analysis files
    sequence_handle.close()
    try:
        string_handle.close()
    except AttributeError:
        pass

    return


def findOffTargets(submit_dict):
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
                feature_type -  str
            file keys
                targets     -   str
                assembly    -   str
                gff         -   str
                output      -   str
            option keys
                index_option -  int
                count       -   int (optional)
                homology    -   int (optional)
                remove_files -  int (optional)
    (see manual.pdf in package root directory for more info)
    :return:
    """

    # Retrieve bash shell location and shell type
    bash_dir = submit_dict['bash']
    shell_type = submit_dict['shell_type']

    # Retrieve job parameters
    window = (submit_dict['start'], submit_dict['end'])
    mismatches = tuple(submit_dict['mismatch'])

    feature_types = submit_dict['feature_type'].split(',')
    for ii in range(len(feature_types)):
        feature_types[ii] = feature_types[ii].strip()

    # Retrieve option selections
    index_option = submit_dict['index_option']
    count = submit_dict['count']
    check_homology = submit_dict['homology']
    remove_files = submit_dict['remove_files']
    small_genome = False  # TODO: Add "-y" flag to Bowtie search for small genomes (check size automatically)

    # Retrieve job file paths
    targets = submit_dict['targets']
    assembly = submit_dict['assembly']
    annotation = submit_dict['gff']

    # Initialize output folder
    out_dir = submit_dict['output']
    out_base = out_dir
    out_num = 0
    folder_check = True
    name_change = False

    # Create output folder, append number if already exists
    while folder_check:
        try:
            os.mkdir(out_dir)
            folder_check = False
            os.mkdir('%s/error_files' % out_dir)
        except WindowsError:
            name_change = True
            out_dir = out_base + ' (%s)' % str(out_num + 1)
            out_num += 1
    os.mkdir(out_dir + '/intermediate_files/')

    # Notify user of output folder name change
    if name_change:
        print('Directory %s already exists, job results written to directory %s' % (out_base, out_dir))

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

    # Find PNA length
    targets_dict = SeqIO.to_dict(SeqIO.parse(targets, 'fasta'))
    length_dict = {}
    sam_dict = {}

    for pna in targets_dict:
        length = str(len(targets_dict[pna]))
        try:
            length_dict[length][pna] = targets_dict[pna]
            if '%s.rev' % pna not in length_dict:  # Add reverse sequence for N to 5' parallel alignment
                length_dict[length]['%s.rev' % pna] = targets_dict[pna][::-1]
        except KeyError:
            length_dict[length] = {}
            length_dict[length][pna] = targets_dict[pna]
            if '%s.rev' % pna not in length_dict:  # Add reverse sequence for N to 5' parallel alignment
                length_dict[length]['%s.rev' % pna] = targets_dict[pna][::-1]

    # Check Bowtie version
    bash_check = ['%s/bash' % bash_dir, '--login', '-c', 'bowtie --version']
    version_call = subprocess.check_output(bash_check).decode('utf-8')
    bt_version = version_call.split('\n')[0].split()[-1]

    # Initialize file naming variable
    filebase = submit_dict['targets'].split('/')[-1].split('.')[0]

    # Iterate on PNA length (for Bowtie run)
    for length in length_dict:
        if len(length_dict) > 1:
            print('Separating FASTA file for PNA length %s...' % length)

        # Construct FASTA file for given PNA length
        temp_targets = out_dir + '/intermediate_files/' + filebase + '.%smer.fa' % length
        SeqIO.write(length_dict[length].values(), temp_targets, 'fasta')
        temp_len = length

        # Add Bowtie flags to set FASTA file input, set seed length, set "all alignments" output
        flags = ['-f', '-l', str(length), '-a']
        if small_genome:
            flags += ['-y']

        # Iterate on number of mismatches
        for mismatch in mismatches:
            # Add Bowtie flag to set maximum mismatches in seed
            flags += ['-n', str(mismatch)]

            # Initialize SAM file naming variable, SAM output folder, error file folder
            sambase = out_dir + '/intermediate_files/%sMM/' % mismatch + filebase
            os.mkdir('%s/intermediate_files/%sMM/' % (out_dir, mismatch))
            os.mkdir('%s/error_files/%sMM/' % (out_dir, mismatch))

            # Set intermediate file names
            temp_sam = sambase + '.%smer.sam' % length
            temp_bam = sambase + '.%smer.bam' % length
            temp_sbam = sambase + '.%smer.sorted.bam' % length
            temp_bed = sambase + '.%smer.bed' % length
            out_bed = sambase + '.bed'

            # Initialize dictionary of temporary intermediate files, to be accessed below
            sam_dict[mismatch] = []
            sam_dict[mismatch] += [[temp_len, temp_sam, temp_bam, temp_sbam, temp_bed, out_bed]]

            # Run Bowtie
            bowtie(temp_targets,
                   index_base,
                   temp_sam,
                   flags,
                   bash_sub=shell_type, bash_dir=bash_dir,
                   err_file='%s/error_files/%sMM/bowtie.%smer.stderr' % (out_dir, mismatch, length),
                   version=bt_version)

    # Check SAMTools version
    bash_check = ['%s/bash' % bash_dir, '--login', '-c', 'samtools --version']
    version_call = subprocess.check_output(bash_check).decode('utf-8')
    samtools_version = version_call.split(' ')[1].split('\n')[0]

    # Initialize BED files for compilation
    bed_handle = None
    bed_dict = {}

    # Iterate through sam_dict entries
    for mismatch in sam_dict:
        first = True

        # Iterate through dictionaries for differing seed (target sequence) length
        for sam in sam_dict[mismatch]:
            # Retrieve dictionary entries, rename variables for clarity
            length = sam[0]
            samfile = sam[1]
            bamfile = sam[2]
            sbamfile = sam[3]
            temp_bed = sam[4]
            out_bed = sam[5]

            # Run SAMtools view, sort, and index
            samtools('view', samfile, outfile=bamfile, flags=['-b'], bash_sub=shell_type, bash_dir=bash_dir,
                     err_file='%s/error_files/%sMM/samtools_view.%smer.stderr' % (out_dir, mismatch, length),
                     version=samtools_version)
            samtools('sort', bamfile, outfile=sbamfile, bash_sub=shell_type, bash_dir=bash_dir,
                     err_file='%s/error_files/%sMM/samtools_sort.%smer.stderr' % (out_dir, mismatch, length),
                     version=samtools_version)
            samtools('index', sbamfile, bash_sub=shell_type, bash_dir=bash_dir,
                     err_file='%s/error_files/%sMM/samtools_index.%smer.stderr' % (out_dir, mismatch, length),
                     version=samtools_version)

            # Check bedtools version
            bash_check = ['%s/bash' % bash_dir, '--login', '-c', 'bedtools --version']
            version_call = subprocess.check_output(bash_check).decode('utf-8')
            bedtools_version = version_call.split(' ')[1].split('\n')[0]

            # Set bedtools function, input files, flags
            bedtools_function = 'window'
            infiles = [sbamfile, annotation]
            flags = ['-l', str(-1 * window[0]), '-r', str(-1 * window[0]), '-sw', '-sm', '-bed']

            # Run bedtools
            bedtools(bedtools_function,
                     infiles,
                     outfile=temp_bed,
                     flags=flags,
                     bash_sub=shell_type, bash_dir=bash_dir,
                     err_file='%s/error_files/%sMM/bedtools.%smer.stderr' % (out_dir, mismatch, length),
                     version=bedtools_version)

            # Compile BED file outputs (for differing seed length) for each mismatch quantity
            if len(sam_dict[mismatch]) > 1:
                if first:
                    bed_handle = open(out_bed, 'w')
                    bed_dict[mismatch] = out_bed
                    first = False
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

    # Filter BED file alignments for those within specified range of feature start
    first = True
    for mismatch in bed_dict:
        # Initialize output file names
        bedfile = bed_dict[mismatch]
        outfile = out_dir + '/' + filebase + '.%smm.out' % mismatch
        if count:
            countfile = out_dir + '/' + filebase + '.%smm.count' % mismatch
        else:
            countfile = None

        # Only run check_homology subfunction if the loop is on the first mismatch quantity
        if check_homology and first:
            homology_outfile = out_dir + '/' + 'homology.out'
            first = False
        else:
            check_homology = False
            homology_outfile = None

        # Run processBedWindow
        pna_tools.processBedWindow(bedfile,
                                   outfile,
                                   dist_filter=window[1],
                                   countfile=countfile,
                                   homology_outfile=homology_outfile,
                                   check_homology=check_homology,
                                   feature_types=feature_types)

    # Remove intermediate files if option specified (useful for large mammalian genome files)
    if remove_files:
        shutil.rmtree('%s/error_files' % out_dir)
        shutil.rmtree('%s/intermediate_files' % out_dir)

    return


def sequenceWarnings(submit_dict):
    """
    Executes the PNA Finder Check Sequence Warnings function using arguments provided by the dictionary submit_dict:
    :param submit_dict:
       check_warnings submit_dict keys
            file keys
                pna_list    -   str
                rnafold     -   str (optional)
                output      -   str
            option keys
                pna_option  -   int
                temp_option -   int (optional)
                rnafold_option - int (optional)
    (see manual.pdf in package root directory for more info)
    :return:
    """

    # Retrieve option selections
    temperature_option = submit_dict['temp_option']
    rnafold_option = submit_dict['rnafold_option']

    # Parse PNA option
    if submit_dict['pna_option'] == 0:
        pna_targets = submit_dict['pna_list']
        target_dict = SeqIO.to_dict(SeqIO.parse(pna_targets, 'fasta'))
    elif submit_dict['pna_option'] == 1:
        pna_list = submit_dict['pna_list']
        target_dict = dict()
        with open(pna_list) as pna_handle:
            for line in pna_handle:
                [pna, value] = line.rstrip().split('\t')
                target_dict[pna] = SeqIO.SeqRecord(seq=value, id=pna)
    else:
        raise ValueError('Something has gone horribly wrong!')

    # Initialize output folder
    out_dir = submit_dict['output']
    out_base = out_dir
    out_num = 0
    folder_check = True
    name_change = False

    # Create output folder, append number if already exists
    while folder_check:
        try:
            os.mkdir(out_dir)
            folder_check = False
            os.mkdir('%s/error_files' % out_dir)
        except WindowsError:
            name_change = True
            out_dir = out_base + ' (%s)' % str(out_num + 1)
            out_num += 1

    # Notify user of output folder name change
    if name_change:
        print('Directory %s already exists, job results written to directory %s' % (out_base, out_dir))

    # Initialize file names
    filebase = submit_dict['pna_list'].split('/')[-1].split('.')[0]
    sequence_outfile = out_dir + '/' + filebase + '.out'

    # Initialize output file, write header
    sequence_handle = open(sequence_outfile, 'w')
    sequence_handle.write('Name\tSequence\tTarget\tWarnings')

    if temperature_option:
        sequence_handle.write('\tTm (DNA/ssPNA)')

    # Perform setup steps for RNA folding analysis
    rna_dict = fold_dir = fold_handle = None
    if rnafold_option and 'rnafold' in submit_dict:
        sequence_handle.write('\t2° Structure Fraction')

        # Take RNA sequences from file, initialize dictionary of sequences
        rna_sequences = submit_dict['rnafold']
        rna_dict = SeqIO.to_dict(SeqIO.parse(rna_sequences, 'fasta'))

        # Initialize folding output folder and files
        fold_dir = out_dir + '/folding_files'
        os.mkdir(fold_dir)
        out_fold = fold_dir + '/fold.out'

        # Append folding column to header
        fold_handle = open(out_fold, 'w')
        fold_handle.write('Name\tSequence\tFold Plot\tEnergy\tTarget Fold Fraction\n')
    else:
        rnafold_option = 0
    sequence_handle.write('\n')

    # Iterate through PNA list
    for pna in list(target_dict.keys()):
        # Retrieve PNA target and sequence from target_dict
        target = str(target_dict[pna].seq)
        sequence = str(Seq.Seq(target).reverse_complement())
        sequence_handle.write('%s\t%s\t%s\t' % (pna, sequence, target))

        # Check for solubility, complement warnings
        warns = sol.checkWarnings(sequence)
        if warns:
            sequence_handle.write(', '.join(warns))

        # Check melting temperature of PNA/ssDNA duplex
        if temperature_option:
            dH, dS, dG, temp_dna, temp_pna = sol.calculateTmPNA(sequence)
            sequence_handle.write('\t%s' % round(temp_pna, 1))

        # Run RNA folding analysis
        if rnafold_option:
            try:
                rna = rna_dict[pna]

                # Look for PNA target sequence in RNA sequence, skip if not found
                start = str(rna.seq).find(target)
                if start == -1:
                    warnings.warn('Target sequence for PNA name "%s" not found in RNA, skipping...' % pna)
                    continue
                else:
                    end = start + len(target)
                    coordinates = (start, end)

                # Create temporary FASTA file to be used in the RNAfold program
                temp_fasta = fold_dir + '/%s_rna.fa' % pna
                SeqIO.write(rna, temp_fasta, 'fasta')

                # Run RNAfold, place outputs in fold_dict dictionary
                fold_dict = pna_tools.rnafold(temp_fasta, coordinates=coordinates, out_dir=fold_dir)

                # Write results to general output file
                sequence_handle.write('\t%s' % round(fold_dict[pna][3], 2))
                fold_handle.write('%s' % pna)
                for item in fold_dict[pna]:
                    fold_handle.write('\t%s' % item)
                fold_handle.write('\n')

                # Remove RNAfold temporary FASTA file
                os.remove(temp_fasta)

            except KeyError:
                warnings.warn('No RNA sequence matching PNA name "%s" found, skipping...' % pna)
        sequence_handle.write('\n')

    # Close analysis files
    sequence_handle.close()
    try:
        fold_handle.close()
    except AttributeError:
        pass

    return
