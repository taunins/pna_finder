import gffutils
from Bio import SeqIO
from collections import defaultdict
import warnings
import subprocess
import os


# Dictionary handling

def fastaToDict(infasta):
    """
    NOTE: This function is obsolete, and should be replaced by using SeqIO.to_dict on the SeqIO.parse iterator

    Takes FASTA file and returns a Python dictionary of IDs as keys and sequences as values
    :param infasta: File path for FASTA file
    :return: Dictionary of FASTA records
    """

    fasta_dict = defaultdict()

    for seq_record in SeqIO.parse(infasta, 'fasta'):
        fasta_dict[seq_record.id] = str(seq_record.seq)

    return fasta_dict


def dictToFasta(fasta_dict, filename, line_length=80):
    """
    NOTE: This function is obsolete, and should be replaced by using SeqIO.write on the SeqIO.parse iterator

    Takes a Python dictionary of IDs as keys and sequences as values and outputs a FASTA file
    :param fasta_dict: Dictionary of FASTA IDs as keys, FASTA sequences as values
    :param filename: Name of FASTA file that will be written with this dictionary
    :param line_length: Length of FASTA file line, 80 as default
    :return:
    """

    if type(fasta_dict) is dict and type(filename) is str:
        out_fasta = open(filename, "w")
        for key in list(fasta_dict.keys()):
            if line_length:
                line_start = 0
                line_end = line_length
                total_lines = (len(fasta_dict[key]) // line_length) + 1
                out_fasta.write('>%s\n' % key)
                for ii in range(total_lines):
                    out_fasta.write('%s\n' % fasta_dict[key][line_start:line_end])
                    line_start = line_end
                    line_end += line_length
            else:
                out_fasta.write('>%s\n%s\n' % (key, fasta_dict[key]))
    else:
        raise TypeError('dictToFasta accepts a dict input (arg 1) and str input file path (arg 2)')

    return


# Match ID to GFF/GTF file

def isAttribute(gene_id, feature, parent=None, child=None, hierarchy=None):
    """
    Takes a gene/protein ID and checks a CDS feature and gene (parent) feature, if supplied, for matching IDs
    :param gene_id: gene/protein ID
    :param feature: feature that the ID will be compared to
    :param parent: parent feature, if applicable
    :param child: child feature, if applicable
    :param hierarchy: a list that gives the order by which an annotation record's different attributes will be tried for
    use as the primary name for a given ID
    :return: True/False, name/None: First returned value is True or False depending on whether the gene/protein ID
    provided matches the feature provided. Second value is the name that has been assigned to the feature, if it is a
    match. If it is not a match, the second argument is None.
    """

    # Initialize dictionaries for attributes
    attribute_dict = {'feature': {'gene': None, 'protein_id': None, 'Dbxref': None,  # GFF attributes
                                  'gene_id': None, 'gene_name': None, 'transcript_id': None, 'tss_id': None},  # GTF attributes
                      'parent': {'gene': None, 'Name': None, 'gene_synonym': None, 'locus_tag': None, 'Dbxref': None,
                                 'gene_id': None, 'gene_name': None, 'transcript_id': None, 'tss_id': None}}

    if not hierarchy:  # Set default hierarchy for feature search
        hierarchy = ['feature.gene', 'parent.gene', 'parent.Name', 'parent.gene_synonym', 'feature.protein_id',
                     'parent.locus_tag', 'feature.Dbxref', 'parent.Dbxref',
                     'feature.gene_id', 'feature.gene_name', 'parent.gene_id', 'parent.gene_name',
                     'feature.transcript_id', 'feature.tss_id', 'parent.transcript_id', 'parent.tss_id']

    if parent and not child:
        pass
    elif child and not parent:  # Reduce to single case based on symmetry of feature type possibilities
        parent = feature
        feature = child
    elif child and parent:
        warnings.warn('isAttribute does not take both parent and child for single feature, using parent as default')
        pass
    else:  # Handle case where feature type cannot be inferred from parent/child features
        attribute_dict['feature'] = {'gene': None, 'protein_id': None, 'Dbxref': None, 'Name': None,
                                     'gene_synonym': None, 'locus_tag': None,
                                     'gene_id': None, 'gene_name': None, 'transcript_id': None, 'tss_id': None}

    # Iterate through feature and parent gene records to search for all identifiers, remove entry from hierarchy if not
    # found
    for attribute in list(attribute_dict['feature'].keys()):
        try:
            attribute_dict['feature'][attribute] = feature[attribute]
        except KeyError:
            try:  # Remove attributes not found in feature from hierarchy list, if they exist in list
                hierarchy.remove('feature.%s' % attribute)
            except ValueError:
                pass

    if parent:
        for attribute in list(attribute_dict['parent'].keys()):
            try:
                attribute_dict['parent'][attribute] = parent[attribute]
            except KeyError:
                hierarchy.remove('parent.%s' % attribute)  # Not necessary to use exception for this remove statement

    # Qualifiers that have been found are iterated through to look for matches
    for entry in hierarchy:
        try:
            [key_1, key_2] = entry.split('.')  # Split hierarchy list entries to retrieve nested dictionary values
            if key_2 == 'Dbxref' and attribute_dict[key_1][key_2]:  # Handle database cross-reference, None/list value
                identifier = [db_id.split(':')[1] for db_id in attribute_dict[key_1][key_2]]
            else:
                identifier = attribute_dict[key_1][key_2]

            if gene_id in identifier:  # If identifier found, name PNA by first entry in hierarchy, return True and name
                name = attribute_dict[hierarchy[0].split('.')[0]][hierarchy[0].split('.')[1]][0]
                return True, name
        except TypeError:  # Handle case in which attribute_dict Dbxref entry is not iterable, though this should not occur
            pass

    # If no match found, return False and None for name
    return False, None


def findID(gff_db, in_list, out_bed,
           feature_types=('CDS',),
           id_type='id',
           full_search=False,
           error_file=None):
    """
    Takes a path to a GFF3 file database (already created through gffutils.createdb function) and a file path for an ID
    list formatted as a single column list) of gene/protein identifiers, and outputs a .bed file of those records
    :param gff_db: A gffutils database that will be searched for features/attributes that match the IDs provided
    :param in_list: A single-column list of gene/protein IDs
    :param out_bed: The path and name of the output BED file
    :param feature_types: The types of features that should be examined. If left as None, the function searches all
    features
    :param id_type: A string variable with possible entries of 'id' (default) and 'position'. Using 'id' searches for
    specific gene identifiers (gene name, database key, etc.) while position searches for feature start index within the
    chromosome.
    :param full_search: A boolean variable that tells the function whether to look for genes using identifiers that are
    not the database key of the gffutils database. If 'True' is passed, the function will first look through database
    keys, then will search other annotation types for the gene/protein IDs for which no match was found.
    :param error_file: Optional output file that gives identifiers and their matching gene name.
    :return:
    """

    # Open gffutils database file
    db = gffutils.FeatureDB(gff_db)

    # Read gene/protein ID file into list
    with open(in_list) as f_handle:
        id_list = [line.rstrip() for line in f_handle]

    # Initialize output file, matches and not_found lists
    bed_handle = open(out_bed, "w")
    matches = []
    not_found = []

    # Perform initial search of database keys (regardless of whether full_search=True)
    if id_type == 'id':
        for gene_id in id_list:
            try:  # Retrieve feature information
                feature = db[gene_id]  # Exception statement will be raised here, if at all
                name = feature.id
                chromosome = feature.seqid
                start = feature.start
                end = feature.end
                score = feature.score
                strand = feature.strand
                matches.append([gene_id, name])

                # Write BED file
                bed_handle.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (chromosome, start, end, name, score, strand))
            except gffutils.exceptions.FeatureNotFoundError:
                # Add gene to not_found list if no matching database key found
                not_found.append(gene_id)
    elif id_type == 'position':  # Add all genes to not_found list if 'position' selected for id_type
        not_found = id_list
    else:
        raise ValueError('ID type not supported, input must be "id" or "position"')

    # Check for full search option and gene/protein IDs that have not been found, or for 'position' id_type
    if (full_search and not_found) or id_type == 'position':
        # Filter database to include only selected feature types
        if feature_types == ('',):
            features = db.all_features()
        else:
            features = db.all_features(featuretype=feature_types)

        # Iterate through database and not_found list
        for feature in features:
            for gene_id in not_found:
                condition = False
                name = None

                # Use isAttribute function to search for matches within genome annotations
                if id_type == 'id':
                    parent = None

                    # Look for parent and child features to expand annotation search
                    try:
                        parent = list(db.parents(feature['ID'][0]))[0]
                    except KeyError:
                        try:
                            parent = list(db.parents(feature['gene_id'][0]))[0]
                        except IndexError:
                            pass
                    except IndexError:
                        pass

                    child = None
                    try:
                        child = list(db.children(feature['ID'][0]))[0]
                    except KeyError:
                        try:
                            child = list(db.children(feature['gene_id'][0]))[0]
                        except IndexError:
                            pass
                    except IndexError:
                        pass

                    # Return True when feature match is found, with a feature name that will be used to name the PNA
                    condition, name = isAttribute(gene_id, feature, parent=parent, child=child)

                # Check for feature start position that matches gene/protein start indices
                elif id_type == 'position':  # TODO: Chromosome ID, more efficient search, genes with same start/end
                    try:
                        if feature.start == int(gene_id) or feature.end == int(gene_id):
                            condition = True
                            try:
                                name = feature['gene'][0]
                            except KeyError:
                                try:
                                    name = feature['gene_id'][0]
                                except KeyError:
                                    name = feature.id
                    except ValueError:
                        warnings.warn('If id_type is "position" list entries must be integers, skipping %s' % gene_id)
                        not_found.remove(gene_id)
                else:
                    raise ValueError('id_type must be either "id" or "position"')

                # Check if feature match is found, print message and write to BED file
                if condition:
                    print('%s matches record for %s' % (gene_id, name))
                    chromosome = feature.seqid
                    start = feature.start
                    end = feature.end
                    score = feature.score
                    strand = feature.strand

                    bed_handle.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (chromosome, start, end, name, score, strand))
                    not_found.remove(gene_id)
                    break

    # Print message if any gene IDs cannot be found
    if not_found:
        print('No matching record found for gene ID(s): %s' % ', '.join(not_found))
    else:
        print('All gene IDs successfully matched to features')

    bed_handle.close()

    # Write matches to error file
    if error_file:
        with open(error_file, 'w') as error_handle:
            error_handle.write('gene_id\tfeature_match_id\n')
            for pair in matches:
                error_handle.write('%s\n' % '\t'.join(pair))
            for gene_id in not_found:
                error_handle.write('%s\tNone\n' % gene_id)
        error_handle.close()

    return


# BEDTools file helpers

def editBed(in_bed, out_bed,
            window=(-5, -5),
            sequence_length=(12,)):
    """
    Takes input BED file (typically of targeted features/sequences produced by findID) and produces new BED file of
    enumerated sequences targeting the full window specified relative to the start codon
    :param in_bed: input BED file name
    :param out_bed: output BED file name
    :param window: the window of bases relative to the start codon for which antisense sequences will be designed
    :param sequence_length: the desired length of antisense sequences
    :return:
    """

    # Import input BED file
    with open(in_bed) as fhandle:
        bed_list = [line.rstrip().split('\t') for line in fhandle]

    # Initialize edited BED file
    bed_handle = open(out_bed, "w")

    # Iterate through lines of input BED file
    for record in bed_list:
        ref_id = record[0]
        name = record[3]
        score = record[4]
        strand = record[5]

        # Handle plus/minus strand for altering BED file start and end indices
        if strand == '+':
            feat_start = int(record[1])
            slide_set = list(range(window[0] - 1, window[1]))

            # Iterate through the user-specified window to retrieve all target sequence indices for given antisense
            # oligomer length within range
            for length in sequence_length:
                rec_num = 0
                for jj in slide_set:
                    start = feat_start + jj
                    end = start + length

                    if len(sequence_length) == 1:
                        feat_id = name + '_%s' % str(rec_num)  # Number the target sequence names
                        rec_num += 1
                    else:
                        feat_id = name + '_%s.%s' % (str(length), str(rec_num))     # Number by PNA length as well
                        rec_num += 1

                    # Write to new BED file
                    bed_handle.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (ref_id, str(start), str(end), feat_id, score, strand))

        elif strand == '-':
            feat_start = int(record[2])
            slide_set = list(range(window[0], window[1] + 1))

            # Iterate through the user-specified window to retrieve all target sequence indices for given antisense
            # oligomer length within range
            for length in sequence_length:
                rec_num = 0
                for jj in slide_set:
                    start = feat_start + -1 * jj
                    end = start - length

                    if len(sequence_length) == 1:
                        feat_id = name + '_%s' % str(rec_num)  # Number the target sequence names
                        rec_num += 1
                    else:
                        feat_id = name + '_%s.%s' % (str(length), str(rec_num))  # Number by PNA length as well
                        rec_num += 1

                    # Write to new BED file
                    bed_handle.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (ref_id, str(end), str(start), feat_id, score, strand))

    bed_handle.close()
    return


def processBedWindow(bedfile, outfile,
                     dist_filter=20,
                     countfile=None,
                     check_homology=False,
                     homology_outfile=None,
                     feature_types=None):
    """
    Processes the output from a BEDTools window function to determine which feature-overlapping alignments are likely
    to cause antisense gene expression/mRNA translation inhibition.
    :param bedfile: Input BED file that will be analyzed
    :param outfile: Output .out file that will display filtered and processed results
    :param dist_filter: The distance downstream from the given feature start where expression/translation inhibition is
    still expected when an antisense molecule binds
    :param countfile: Off-target count file that tabulates the number of inhibitory off-targets for each antisense
    sequence
    :param check_homology: Allows user to check for matches between target gene name (designated by the name in the
    FASTA file prior to the underscore) and a given off-target annotation
    :param homology_outfile: Output file
    :param feature_types: The type of feature for which the user wants to find off-targets
    :return:
    """

    # Initialize output file and write header
    out_handle = open(outfile, "w")
    out_handle.write('PNA ID\tChromosome ID\tAlignment Strand\tAlignment Start\tFeature Start\tRelative Position\t'
                     'Feature Info\n')

    # Handle feature checking function
    feature_check = True
    if feature_types is None:
        feature_types = ['CDS']
    elif feature_types == ['']:
        # Skip feature filtering if feature_type variable left blank
        feature_check = False

    # Initialize dictionaries for off-target counting and duplicates, homology list
    ot_count = {}
    count_duplicates = {}
    potential_homology = []

    # Initialize off-target output and homology output handle variables
    ot_handle = None
    homology_handle = None
    line_number = 2

    # Check whether to filter feature alignments from downstream
    if dist_filter == '':
        distance_pass = True
    else:
        distance_pass = False

    # Open BED file, iterate line by line
    with open(bedfile, 'r') as bed_handle:
        for line in bed_handle:
            # Split BED file line into columns, retrieve feature type for subsequent filter step
            record = line.split('\t')
            feature_type = record[8]

            # Skip loop if alignment feature type does not match input feature types
            if feature_type in feature_types:
                pass
            elif feature_check:
                continue

            # Retrieve relevant BED file fields for output file
            pna_name = record[3]
            ref_id = record[0]
            strand = record[5]
            feat_info = record[14].split('\n')[0]   # Remove line break from feature info

            # Check alignment strand, use appropriate feature start and end to calculate alignment distance
            if strand == '+':
                align_start = record[1]
                feat_start = record[9]
                distance = int(align_start) - int(feat_start) + 1
            elif strand == '-':
                align_start = record[2]
                feat_start = record[10]
                distance = int(feat_start) - int(align_start)
            else:
                raise IndexError('Error reading genome strand id')

            # Check whether feature is within distance filter limit
            if distance_pass or distance < dist_filter:
                # Compile information for output file
                out_line = [pna_name, ref_id, strand, align_start, feat_start, str(distance), feat_info]

                # Keep count of number of off-targets for each PNA
                if countfile:
                    if pna_name not in list(ot_count.keys()):
                        ot_count[pna_name] = 1
                        count_duplicates[pna_name] = {}     # Initialize new dictionary within count_duplicates dictionary for PNA name
                        count_duplicates[pna_name][ref_id] = [feat_start]   # Add feature start index (within list) to chromosome/reference dictionary entry
                    elif ref_id not in list(count_duplicates[pna_name].keys()):
                        count_duplicates[pna_name][ref_id] = []     # Add empty list to new chromosome/reference entry for existing PNA name dictionary

                    # Count off-target alignment for new feature start. Do not count if it has the same start index as
                    # another feature on the same chromosome/reference (this avoids multiple counting of human alternate RNA
                    # transcripts)
                    if feat_start not in count_duplicates[pna_name][ref_id]:
                        ot_count[pna_name] += 1
                        count_duplicates[pna_name][ref_id] += [feat_start]

                # Check for potential PNA target feature homology in strain by using a simple substring search
                if check_homology:  # TODO: Improve homology check by parsing feature IDs, allowing keyword inputs
                    if pna_name.split('_')[0].upper() in feat_info.upper():
                        potential_homology += [[pna_name, line_number, feat_info]]

                # Write to outfile
                for entry in out_line:
                    out_handle.write('%s\t' % str(entry))
                out_handle.write('\n')

                line_number += 1

    out_handle.close()

    # Write to off-target count file
    if countfile:
        ot_handle = open(countfile, 'w')
        for key in list(ot_count.keys()):
            ot_handle.write('%s\t%s\n' % (key, ot_count[key]))

    try:
        ot_handle.close()
    except AttributeError:
        pass

    # Write homology output to file or print
    if check_homology and homology_outfile:
        homology_handle = open(homology_outfile, 'w')
        homology_handle.write('PNA ID\tOutfile Line Number\tFeature Info\n')
        for entry in potential_homology:
            homology_handle.write('%s\t%s\t%s\n' % (entry[0], entry[1], entry[2]))
    elif check_homology and not homology_outfile:
        print('Potential alignments to homologous features:')
        for entry in potential_homology:
            print('%s\t%s\n' % (entry[0], entry[2]))

    try:
        homology_handle.close()
    except AttributeError:
        pass

    return


# Miscellaneous functions

def revComp(sequence):
    """
    NOTE: This function is obsolete, and should be replaced by using reverse_complement() on a Bio.Seq object

    Returns the reverse complement of a DNA or RNA sequence
    :param sequence: DNA/RNA sequence
    :return: Reverse complement sequence string
    """

    # Standardize entry to be uppercase
    try:
        sequence = sequence.upper()
    except AttributeError:
        raise TypeError('rev_comp takes a nucleobase sequence string')

    # Initialize output variable
    rc_sequence = ''

    # Check for mixed DNA and RNA sequence, set nucleic acid type
    if 'T' in sequence and 'U' in sequence:
        raise ValueError('Mixed RNA and DNA sequences')
    elif 'U' in sequence:
        seq_type = 'RNA'
    else:
        seq_type = 'DNA'

    # Iterate backwards through input sequence, append complement to reverse complement output
    for letter in sequence[::-1]:
        if letter == 'A':
            if seq_type == 'DNA':
                rc_sequence += 'T'
            elif seq_type == 'RNA':
                rc_sequence += 'U'
        elif letter == 'T' or letter == 'U':
            rc_sequence += 'A'
        elif letter == 'C':
            rc_sequence += 'G'
        elif letter == 'G':
            rc_sequence += 'C'
        else:
            raise ValueError('Non-RNA or non-DNA base included')

    return rc_sequence


def rnafold(fasta, outfile=None, coordinates=None, out_dir=None):
    """
    Returns either an tab-delimited output file or a dictionary of sequence names keyed to lists of the following
    values: [0] mRNA sequence, [1] folding plot, [2] free energy of folding, [3] (optional) secondary structure fraction
    of sub-sequence (as designated by coordinate list input)
    :param fasta: mRNA sequence file
    :param outfile: (optional) path to output file
    :param coordinates: (optional) coordinates of sub-sequence for which fractional folding will be returned
    :param out_dir: (optional) directory in which to place the .ps and .ss RNA folding files
    :return: Dictionary with FASTA file id as the key, and a list of [sequence, fold_plot, energy] as the value
    """

    # Change working directory to the input out directory if provided. This will place intermediate folding files into
    # this directory
    if out_dir:
        try:
            os.chdir(out_dir)
        except FileNotFoundError:
            warnings.warn('Output directory does not exist, file outputs placed in current working directory')

    # Pass FASTA file to RNAfold program in command line, process output
    fold_output = subprocess.check_output(['RNAfold', fasta]).decode('utf-8')
    out_dict = {}
    out_list = fold_output.split('\r\n')

    # Process output list and place variables in output dictionary
    for ii in range(len(out_list)):
        try:
            if out_list[ii][0] == '>':
                name = out_list[ii].split(' ')[0][1:]
                sequence = out_list[ii + 1]
                fold_plot = out_list[ii + 2].split(' ')[0]
                energy = float(out_list[ii + 2].split(' ', 1)[1][1:-1].strip())

                out_dict[name] = [sequence, fold_plot, energy]
                if coordinates:
                    try:
                        target_fold = fold_plot[coordinates[0]:coordinates[1]]
                        fraction = (target_fold.count('(') + target_fold.count(')')) / len(target_fold)
                        out_dict[name] += [fraction]
                    except TypeError or IndexError:
                        Warning('Input variable "coordinates" must be a list of two integers')
        except IndexError:
            pass

    # Write output to file if provided
    if outfile:
        with open(outfile, 'w') as out_handle:
            out_handle.write('Name\tSequence\tFold Plot\tEnergy')
            if coordinates:
                out_handle.write('\tTarget Fold Fraction')
            out_handle.write('\n')

            for name in out_dict.keys():
                out_handle.write('%s' % name)
                for element in out_dict[name]:
                    out_handle.write('\t%s' % element)
                out_handle.write('\n')

    return out_dict


def truncateFasta(fasta, gff_db, out_fasta, out_gff, out_info=None,
                  feature_types=('CDS',),
                  window=(-50, 50)):
    """
    Shortens FASTA file to relevant window of bases for PNA/ASO off-target searching
    :param fasta:
    :param gff_db:
    :param out_fasta:
    :param out_gff:
    :param out_info:
    :param feature_types:
    :param window:
    :return:
    """

    fasta_dict = fastaToDict(fasta)
    truncate_dict = {}

    db = gffutils.FeatureDB(gff_db)

    out_gff = open(out_gff, "w")
    out_info = open(out_info, "w")

    if out_info:
        info = (fasta, ', '.join(feature_types), str(window))
        info_text = 'Truncation settings: %s\nFeature Types: %s\nStart Range: %s' % info
        out_info.write(info_text)
        out_info.close()

    if feature_types == ['']:
        db_list = list(db.all_features())
    else:
        db_list = list(db.all_features(featuretype=feature_types))

    prev_feat_info = None
    for feature in db_list:
        feat_info = (feature.seqid, feature.start, feature.end, feature.strand)
        if feat_info == prev_feat_info:
            print('%s is a duplicate feature, skipping...' % feature.id)
            continue

        if feature.strand == '+':
            start = feature.start + window[0]
            end = feature.start + window[1]
            feature_start = -1 * window[0]
            feature_end = window[1] - window[0]
            if feature_start < 0:
                feature_start = 1
        elif feature.strand == '-':
            start = feature.end + window[0]
            end = feature.end + window[1]
            feature_start = 1
            feature_end = -1 * window[0]
            if feature_start < 0:
                feature_start = 0
        else:
            warnings.warn('Skipping feature %s, no strand designation' % feature.id)
            continue

        truncate_dict[feature.seqid + '_%s' % feature.id] = fasta_dict[feature.seqid][start:end]

        feature_list = str(feature).split('\t')
        feature_list[0] = feature.seqid + '_%s' % feature.id

        feature_list[3] = str(feature_start)
        feature_list[4] = str(feature_end)

        out_gff.write('\t'.join(feature_list) + '\n')

        prev_feat_info = feat_info

    out_gff.close()
    dictToFasta(truncate_dict, out_fasta)
    return


def checkHomology(outfile, homology_outfile=None):
    """
    Function to check homology outside of the processBedWindow function
    :return:
    """

    potential_homology = []
    line_number = 2
    h_handle = None

    with open(outfile) as fhandle:
        for line in fhandle:
            record = line.split('\t')
            pna_name = record[0]
            gene_target = pna_name.split('_')[0]
            feat_info = record[6]

            if gene_target.upper() in feat_info.upper():
                potential_homology += [[pna_name, line_number, feat_info]]

            line_number += 1

    if homology_outfile:
        h_handle = open(homology_outfile, 'w')
        for item in potential_homology:
            h_handle.write('%s\t%s\t%s\n' % (item[0], item[1], item[2]))
    else:
        print('Potential alignments to homologous features:')
        for item in potential_homology:
            print('%s\t%s\t%s\n' % (item[0], item[1], item[2]))

    try:
        h_handle.close()
    except AttributeError:
        pass

    return
