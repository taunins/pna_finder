import gffutils
from Bio import SeqIO
from collections import defaultdict
import warnings


def bashToWindows(path):
    """
    Provides workaround for Windows bash shell path handling, allowing Python subprocess.call functions to work on
    bash shells in Windows
    :param path: Cygwin or Windows bash path (e.g. str starting with '/cygdrive/c/...' or '/mnt/c/...'
    :return: Windows path (e.g. str starting with 'C:/...')
    """

    if type(path) is str:
        path_list = path.split('/')
        if path_list[1] == 'cygdrive':
            winpath = path_list[2].upper() + ':/' + '/'.join(path_list[3:])
            return winpath
        elif path_list[1] == 'mnt':
            winpath = path_list[2].upper() + ':/' + '/'.join(path_list[3:])
            return winpath
        else:
            raise ValueError('bashToWindows takes Windows bash shell file path, e.g. /cygdrive/c/... or /mnt/c/...')
    else:
        raise TypeError('bashToWindows takes string path input')


def windowsToBash(path, shell_type):
    """
    Provides workaround for Windows bash shell path handling, allowing Python subprocess.call functions to work on
    bash shells in Windows
    :param path: Windows file path
    :param shell_type: 'cygwin' or 'winbash'
    :return:
    """

    if type(path) is str:
        if '/' in path:
            path_list = path.split('/')
        elif '\\' in path:
            path_list = path.split('\\')
        else:
            raise ValueError('windowsToBash takes Windows file path, e.g. C:/...')

        if shell_type == 'cygwin':
            try:
                bashpath = '/cygdrive/' + path_list[0][0].lower() + '/' + '/'.join(path_list[1:])
                return bashpath
            except IndexError:
                raise ValueError('windowsToBash takes Windows file path, e.g. C:/...')
        elif shell_type == 'winbash':
            try:
                bashpath = '/mnt/' + path_list[0][0].lower() + '/' + '/'.join(path_list[1:])
                return bashpath
            except IndexError:
                raise ValueError('windowsToBash takes Windows file path, e.g. C:/...')
        else:
            raise ValueError('windowsToBash takes Windows file path, e.g. C:/...')
    else:
        raise TypeError('windowsToBash takes string path input')


def fastaToDict(infasta):
    """
    Takes FASTA file and returns a Python dictionary of IDs as keys and sequences as values
    :param infasta: File path for FASTA file
    :return: fastaDict: a dictionary of FASTA records
    """

    fastaDict = defaultdict(list)

    for seq_record in SeqIO.parse(infasta, 'fasta'):
        fastaDict[seq_record.id] = str(seq_record.seq)

    return fastaDict


def dictToFasta(inDict, filename):
    """
    Takes a Python dictionary of IDs as keys and sequences as values and outputs a FASTA file
    :param inDict: Dictionary of FASTA IDs as keys, FASTA sequences as values
    :param filename: Name of FASTA file that will be written with this dictionary
    :return:
    """

    if type(inDict) is dict and type(filename) is str:
        outfasta = open(filename, "w")
        for key in list(inDict.keys()):
            outfasta.write('>%s\n%s\n' % (key, inDict[key]))
    else:
        raise TypeError('dictToFasta accepts a dict input (arg 1) and str input file path (arg 2)')


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
    att_dict = defaultdict()
    att_dict['feature'] = {'gene': None, 'protein_id': None, 'Dbxref': None,    # GFF attributes
                           'gene_id': None, 'gene_name': None, 'transcript_id': None, 'tss_id': None}   # GTF attributes
    att_dict['parent'] = {'gene': None, 'Name': None, 'gene_synonym': None, 'locus_tag': None, 'Dbxref': None,
                          'gene_id': None, 'gene_name': None, 'transcript_id': None, 'tss_id': None}

    if not hierarchy:   # Set default hierarchy
        hierarchy = ['feature.gene', 'parent.gene', 'parent.Name', 'parent.gene_synonym', 'feature.protein_id',
                     'parent.locus_tag', 'feature.Dbxref', 'parent.Dbxref',
                     'feature.gene_id', 'feature.gene_name', 'feature.transcript_id', 'feature.tss_id',
                     'parent.gene_id', 'parent.gene_name', 'parent.transcript_id', 'parent.tss_id']

    if parent and not child:
        pass
    elif child and not parent:  # Reduce to single case based on symmetry of feature type possibilities
        parent = feature
        feature = child
    elif child and parent:
        raise ValueError('isAttribute does not take both parent and child for single feature record')
    else:   # Handle case where feature type cannot be inferred from parent/child features
        att_dict['feature'] = {'gene': None, 'protein_id': None, 'Dbxref': None, 'Name': None, 'gene_synonym': None,
                               'locus_tag': None,
                               'gene_id': None, 'gene_name': None, 'transcript_id': None, 'tss_id': None}
        att_dict['parent'] = {}

    # Iterate through feature and parent gene records to search for all identifiers
    for key in list(att_dict['feature'].keys()):
        try:
            att_dict['feature'][key] = feature[key]
        except KeyError:
            hierarchy.remove('feature.%s' % key)

    if parent:
        for key in list(att_dict['parent'].keys()):
            try:
                att_dict['parent'][key] = parent[key]
            except KeyError:
                hierarchy.remove('parent.%s' % key)

    # Qualifiers that have been found are iterated through to look for matches
    for entry in hierarchy:
        try:
            for identifier in att_dict[entry.split('.')[0]][entry.split('.')[1]]:
                if gene_id in identifier:
                    name = att_dict[hierarchy[0].split('.')[0]][hierarchy[0].split('.')[1]][0]
                    return True, name
        except TypeError:
            pass

    return False, None


def findID(gff_db, in_list, out_bed, feature_types=None):
    """
    Takes a path to a GFF3 file database (already created through gffutils.createdb function) and a file path for an ID
    list formatted as a single column list) of gene/protein identifiers, and outputs a .bed file of those records
    :param gff_db: A gffutils database that will be searched for features/attributes that match the IDs provided
    :param in_list: A single-column list of gene/protein IDs
    :param out_bed: The path and name of the output BED file
    :param feature_types: The types of features that should be examined. If left as None, the function searches all
    features
    :return:
    """

    if feature_types is None:
        feature_types = ['CDS']

    db = gffutils.FeatureDB(gff_db)
    with open(in_list) as fhandle:
        id_list = [line.rstrip() for line in fhandle]

    outfile = open(out_bed, "w")

    if feature_types == ['']:
        db_list = list(db.all_features())
    else:
        db_list = list(db.all_features(featuretype=feature_types))

    for gene_id in id_list:
        for feature in db_list:

            parent = None
            try:
                parent = list(db.parents(feature['ID'][0]))[0]
            except KeyError:
                try:
                    parent = list(db.parents(feature['gene_id'][0]))[0]
                except:
                    pass
            except:
                pass

            child = None
            try:
                child = list(db.children(feature['ID'][0]))[0]
            except KeyError:
                try:
                    child = list(db.children(feature['gene_id'][0]))[0]
                except:
                    pass
            except:
                pass

            condition, name = isAttribute(gene_id, feature, parent=parent, child=child)
            if condition:
                print('%s matches record for %s' % (gene_id, name))
                chrom = feature.seqid
                start = feature.start
                end = feature.end
                strand = feature.strand

                outfile.write('%s\t%s\t%s\t%s\t%s\n' % (chrom, start, end, name, strand))
                db_list.remove(feature)
                break
            else:
                db_list = db_list[1:] + [db_list[0]]

        else:
            print('No matching record found for gene ID %s' % gene_id)

    outfile.close()


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
        refID = record[0]
        name = record[3]
        strand = record[4]

        if strand == '+':
            feat_start = int(record[1])
            slide_set = list(range(window[0] - 1, window[1]))

            for jj in slide_set:
                start = feat_start + jj
                end = start + sequence_length
                id = name + '_%s' % str(rec_num)

                outfile.write('%s\t%s\t%s\t%s\t%s\n' % (refID, str(start), str(end), id, strand))
                rec_num += 1

        elif strand == '-':
            feat_start = int(record[2])
            slide_set = list(range(window[0], window[1] + 1))

            for jj in slide_set:
                start = feat_start + -1 * jj
                end = start - sequence_length
                id = name + '_%s' % str(rec_num)

                outfile.write('%s\t%s\t%s\t%s\t%s\n' % (refID, str(end), str(start), id, strand))
                rec_num += 1

    outfile.close()


# FIND OFF TARGETS TOOLS

def process_bedwindow(infile, outfile,
                      dist_filter=20,
                      OT_count_file=None,
                      self_OT=False,
                      feature_types=None):
    """
    Processes the output from a BEDTools window function to determine which feature-overlapping alignments are likely
    to cause antisense gene expression/mRNA translation inhibition.
    :param infile: Input BED file that will be analyzed
    :param outfile: Output .out file that will display filtered and processed results
    :param dist_filter: The distance downstream from the given feature start where expression/translation inhibition is
    still expected when an antisense molecule binds
    :param OT_count_file: Off-target count file that tabulates the number of inhibitory off-targets for each antisense
    sequence
    :param self_OT: True/False depending on whether the off-target genome is that which the antisense sequence
    originates from
    :param feature_types: The type of feature for which the user wants to find off-targets
    :return:
    """

    out = open(outfile, "w")
    out.write('PNA ID\tChromosome ID\tAlignment Strand\tAlignment Start\tFeature Start\tAlignment to STC\tFeature '
              'Info\n')
    OT_count = defaultdict(list)

    feature_check = True
    if feature_types is None:
        feature_types = ['CDS']
    elif feature_types == ['']:
        feature_check = False

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

            seq_name = record[3]
            ref_id = record[0]
            strand = record[5]
            feat_info = record[14].split('\n')[0]

            if strand == '+':
                align_start = record[1]
                feat_start = record[9]
                STC_dist = int(align_start) - int(feat_start) + 1
            elif strand == '-':
                align_start = record[2]
                feat_start = record[10]
                STC_dist = int(feat_start) - int(align_start)
            else:
                raise IndexError('Error reading genome strand id')

            outline = [seq_name, ref_id, strand, align_start, feat_start, str(STC_dist), feat_info]

            if STC_dist < dist_filter:
                # keep count of number of off-targets for each PNA
                if seq_name not in list(OT_count.keys()):
                    if self_OT:
                        OT_count[seq_name] = 0
                    elif not self_OT:
                        OT_count[seq_name] = 1
                    else:
                        warnings.warn('self_OT must be True or False, defaulting to False...')
                        self_OT = False
                else:
                    OT_count[seq_name] += 1

                # write outfile
                for item in outline:
                    out.write('%s\t' % str(item))
                out.write('\n')

    if OT_count_file:
        ot_handle = open(OT_count_file, 'w')
        for key in list(OT_count.keys()):
            ot_handle.write('%s\t%s\n' % (key, OT_count[key]))
    else:
        for key in list(OT_count.keys()):
            print(('%s\t%s' % (key, OT_count[key])))

    out.close()