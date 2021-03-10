import gffutils
from Bio import SeqIO
from collections import defaultdict
import warnings
import subprocess


# Dictionary handling

def fastaToDict(infasta):
    """
    Takes FASTA file and returns a Python dictionary of IDs as keys and sequences as values
    :param infasta: File path for FASTA file
    :return: fastaDict: a dictionary of FASTA records
    """

    fasta_dict = defaultdict()

    for seq_record in SeqIO.parse(infasta, 'fasta'):
        fasta_dict[seq_record.id] = str(seq_record.seq)

    return fasta_dict


def dictToFasta(fasta_dict, filename, line_length=80):
    """
    Takes a Python dictionary of IDs as keys and sequences as values and outputs a FASTA file
    :param fasta_dict: Dictionary of FASTA IDs as keys, FASTA sequences as values
    :param filename: Name of FASTA file that will be written with this dictionary
    :param line_length:
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


# Match ID to GFF/GTF file

def isAttribute(gene_id, feature, parent=None, child=None, hierarchy=None):
    """
    Takes a gene/protein ID and checks a CDS feature and gene (parent) feature, if supplied, for matching IDs
    :param gene_id: gene/protein ID
    :param feature: feature that the ID will be compared to
    :param parent: parent feature, if applicable
    :param child:
    :param hierarchy: a list that gives the order by which an annotation record's different attributes will be tried for
    use as the primary name for a given ID
    :return: True/False, name/None: First returned value is True or False depending on whether the gene/protein ID
    provided matches the feature provided. Second value is the name that has been assigned to the feature, if it is a
    match. If it is not a match, the second argument is None.
    """

    # Initialize dictionaries for attributes
    attribute_dict = {'feature': {'gene': None, 'protein_id': None, 'Dbxref': None,  # GFF attributes
                                  'gene_id': None, 'gene_name': None, 'transcript_id': None, 'tss_id': None},
                      # GTF attributes
                      'parent': {'gene': None, 'Name': None, 'gene_synonym': None, 'locus_tag': None, 'Dbxref': None,
                                 'gene_id': None, 'gene_name': None, 'transcript_id': None, 'tss_id': None}}

    if not hierarchy:  # Set default hierarchy
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

    # Iterate through feature and parent gene records to search for all identifiers
    for key in list(attribute_dict['feature'].keys()):
        try:
            attribute_dict['feature'][key] = feature[key]
        except KeyError:
            try:
                hierarchy.remove('feature.%s' % key)
            except ValueError:
                pass

    if parent:
        for key in list(attribute_dict['parent'].keys()):
            try:
                attribute_dict['parent'][key] = parent[key]
            except KeyError:
                hierarchy.remove('parent.%s' % key)

    # Qualifiers that have been found are iterated through to look for matches
    for entry in hierarchy:
        try:
            [key_1, key_2] = entry.split('.')
            if key_2 == 'Dbxref' and attribute_dict[key_1][key_2]:
                identifier = [db_id.split(':')[1] for db_id in attribute_dict[key_1][key_2]]
            else:
                identifier = attribute_dict[key_1][key_2]

            if gene_id in identifier:
                name = attribute_dict[hierarchy[0].split('.')[0]][hierarchy[0].split('.')[1]][0]
                return True, name
        except TypeError:
            pass

    return False, None


def findID(gff_db, in_list, out_bed, feature_types=('CDS',), id_type='id', full_search=False, error_file=None):
    """
    Takes a path to a GFF3 file database (already created through gffutils.createdb function) and a file path for an ID
    list formatted as a single column list) of gene/protein identifiers, and outputs a .bed file of those records
    :param gff_db: A gffutils database that will be searched for features/attributes that match the IDs provided
    :param in_list: A single-column list of gene/protein IDs
    :param out_bed: The path and name of the output BED file
    :param feature_types: The types of features that should be examined. If left as None, the function searches all
    features
    :param id_type:
    :param full_search:
    :return:
    """

    db = gffutils.FeatureDB(gff_db)

    with open(in_list) as f_handle:
        id_list = [line.rstrip() for line in f_handle]

    outfile = open(out_bed, "w")
    not_found = []
    matches = []

    if id_type == 'id':
        for gene_id in id_list:
            try:
                feature = db[gene_id]
                name = feature.id
                chromosome = feature.seqid
                start = feature.start
                end = feature.end
                score = feature.score
                strand = feature.strand
                matches.append([gene_id, name])

                outfile.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (chromosome, start, end, name, score, strand))
            except gffutils.exceptions.FeatureNotFoundError:
                not_found.append(gene_id)
    else:
        not_found = id_list

    if (full_search and not_found) or id_type == 'position':
        if feature_types == ('',):
            features = list(db.all_features())
        else:
            features = list(db.all_features(featuretype=feature_types))

        for gene_id in not_found:
            for feature in features:
                condition = False
                name = None
                if id_type == 'id':
                    parent = None
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

                    condition, name = isAttribute(gene_id, feature, parent=parent, child=child)

                elif id_type == 'position':
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
                else:
                    raise ValueError('id_type must be either "id" or "position"')

                if condition:
                    print('%s matches record for %s' % (gene_id, name))
                    chromosome = feature.seqid
                    start = feature.start
                    end = feature.end
                    score = feature.score
                    strand = feature.strand

                    outfile.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (chromosome, start, end, name, score, strand))
                    not_found.remove(gene_id)
                    break

    if not_found:
        print('No matching record found for gene ID(s): %s' % ', '.join(not_found))
    else:
        print('All gene IDs successfully matched to features')

    outfile.close()

    if error_file:
        with open(error_file, 'w') as error_handle:
            error_handle.write('gene_id\tfeature_match_id\n')
            for pair in matches:
                error_handle.write('%s\n' % '\t'.join(pair))
            for gene_id in not_found:
                error_handle.write('%s\tNone\n' % gene_id)
        error_handle.close()


# BEDTools file helpers

def editBed(in_bed, out_bed, window=(-5, -5), sequence_length=12):
    """
    Takes input BED file (typically of targeted features/sequences produced by findID) and produces new BED file of
    enumerated sequences targeting the full window specified relative to the start codon
    :param in_bed: input BED file name
    :param out_bed: output BED file name
    :param window: the window of bases relative to the start codon for which antisense sequences will be designed
    :param sequence_length: the desired length of antisense sequences
    :return:
    """

    # Import bed file
    with open(in_bed) as fhandle:
        bed_list = [line.rstrip().split('\t') for line in fhandle]

    outfile = open(out_bed, "w")

    for record in bed_list:
        rec_num = 0
        ref_id = record[0]
        name = record[3]
        score = record[4]
        strand = record[5]

        if strand == '+':
            feat_start = int(record[1])
            slide_set = list(range(window[0] - 1, window[1]))

            for jj in slide_set:
                start = feat_start + jj
                end = start + sequence_length
                feat_id = name + '_%s' % str(rec_num)

                outfile.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (ref_id, str(start), str(end), feat_id, score, strand))
                rec_num += 1

        elif strand == '-':
            feat_start = int(record[2])
            slide_set = list(range(window[0], window[1] + 1))

            for jj in slide_set:
                start = feat_start + -1 * jj
                end = start - sequence_length
                feat_id = name + '_%s' % str(rec_num)

                outfile.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (ref_id, str(end), str(start), feat_id, score, strand))
                rec_num += 1

    outfile.close()


def processBedWindow(infile, outfile,
                     dist_filter=20,
                     ot_count_file=None,
                     check_homology=False,
                     homology_outfile=None,
                     feature_types=None):
    """
    Processes the output from a BEDTools window function to determine which feature-overlapping alignments are likely
    to cause antisense gene expression/mRNA translation inhibition.
    :param infile: Input BED file that will be analyzed
    :param outfile: Output .out file that will display filtered and processed results
    :param dist_filter: The distance downstream from the given feature start where expression/translation inhibition is
    still expected when an antisense molecule binds
    :param ot_count_file: Off-target count file that tabulates the number of inhibitory off-targets for each antisense
    sequence
    :param check_homology:
    :param homology_outfile:
    :param feature_types: The type of feature for which the user wants to find off-targets
    :return:
    """

    out = open(outfile, "w")
    out.write('PNA ID\tChromosome ID\tAlignment Strand\tAlignment Start\tFeature Start\tAlignment to STC\tFeature '
              'Info\n')
    ot_count = {}

    feature_check = True
    if feature_types is None:
        feature_types = ['CDS']
    elif feature_types == ['']:
        feature_check = False

    count_duplicates = {}
    potential_homology = []

    ot_handle = None
    h_handle = None
    line_number = 2

    if dist_filter == '':
        distance_pass = True
    else:
        distance_pass = False

    log_mismatch = False  # TODO
    with open(infile) as fhandle:
        for line in fhandle:
            record = line.split('\t')
            feature_type = record[8]

            # check for only input feature types, if feature types input not left blank
            if feature_type in feature_types:
                pass
            elif feature_check:
                continue
            else:
                pass

            pna_name = record[3]
            ref_id = record[0]
            strand = record[5]
            feat_info = record[14].split('\n')[0]

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

            if distance_pass or distance < dist_filter:
                out_line = [pna_name, ref_id, strand, align_start, feat_start, str(distance), feat_info]

                # keep count of number of off-targets for each PNA
                if pna_name not in list(ot_count.keys()):
                    ot_count[pna_name] = 1
                    count_duplicates[pna_name] = {}
                    count_duplicates[pna_name][ref_id] = [feat_start]
                elif ref_id not in list(count_duplicates[pna_name].keys()):
                    count_duplicates[pna_name][ref_id] = []

                if feat_start not in count_duplicates[pna_name][ref_id]:
                    ot_count[pna_name] += 1
                    count_duplicates[pna_name][ref_id] += [feat_start]

                if check_homology:  # TODO: Improve this
                    if pna_name.split('_')[0].upper() in feat_info.upper():
                        potential_homology += [[pna_name, line_number, feat_info]]

                # write outfile
                for item in out_line:
                    out.write('%s\t' % str(item))
                out.write('\n')

                line_number += 1

    if ot_count_file:
        ot_handle = open(ot_count_file, 'w')
        for key in list(ot_count.keys()):
            ot_handle.write('%s\t%s\n' % (key, ot_count[key]))
    else:
        for key in list(ot_count.keys()):
            print(('%s\t%s' % (key, ot_count[key])))

    if check_homology and homology_outfile:
        h_handle = open(homology_outfile, 'w')
        h_handle.write('PNA ID\tOutfile Line Number\tFeature Info\n')
        for item in potential_homology:
            h_handle.write('%s\t%s\t%s\n' % (item[0], item[1], item[2]))
    elif check_homology and not homology_outfile:
        print('Potential alignments to homologous features:')
        for item in potential_homology:
            print('%s\t%s\n' % (item[0], item[2]))

    out.close()

    try:
        ot_handle.close()
    except AttributeError:
        pass

    try:
        h_handle.close()
    except AttributeError:
        pass

    return


# Miscellaneous functions

def revComp(sequence):
    """

    :param sequence:
    :return:
    """

    try:
        sequence = sequence.upper()
    except AttributeError:
        raise TypeError('rev_comp takes a nucleobase sequence string')

    rc_sequence = ''

    if 'T' in sequence and 'U' in sequence:
        raise ValueError('Mixed RNA and DNA sequences')
    elif 'U' in sequence:
        seq_type = 'RNA'
    else:
        seq_type = 'DNA'

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


def rnafold(fasta, outfile=None, coordinates=None):
    """
    Returns either an tab-delimited output file or a dictionary of sequence names keyed to lists of the following
    values: [0] mRNA sequence, [1] folding plot, [2] free energy of folding, [3] (optional) secondary structure fraction
    of sub-sequence (as designated by coordinate list input)
    :param fasta: mRNA sequence file
    :param outfile: (optional) path to output file
    :param coordinates: (optional) coordinates of sub-sequence for which fractional folding will be returned
    :return:
    """

    fold_output = subprocess.check_output(['RNAfold', fasta]).decode('utf-8')
    out_dict = {}
    out_list = fold_output.split('\r\n')

    for ii in range(len(out_list)):
        try:
            if out_list[ii][0] == '>':
                name = out_list[ii].split(' ')[0][1:]
                sequence = out_list[ii + 1]
                fold_plot = out_list[ii + 2].split(' ')[0]
                energy = float(out_list[ii + 2].split(' ')[1][1:-1])

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
    else:
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
