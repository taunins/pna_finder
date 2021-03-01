import subprocess


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


def bowtie_build(fasta, refpath, err_file='bowtie-build.stderr',
                 version=None,
                 bash_sub=False, bash_dir='/',):
    """
    Construct bowtie indices (.ebwt files)
    :param fasta: Path to FASTA containing reference sequences
    :param refpath: Path to index base (e.g. C:/directory/e_coli_MG1655)
    :param err_file: File path to where output is written
    :param version: Enforces bowtie version number (e.g., '1.2.3')
    :param bash_sub: False/cygwin/winbash
    :param bash_dir: Location of bash.exe, if bash_sub not False
    :return:
    """

    try:
        if not bash_sub:
            version_call = subprocess.check_output(['bowtie-build', '--version']).decode('utf-8')
        else:
            bash_check = ['%s/bash' % bash_dir, '--login', '-c', 'bowtie-build --version']
            version_call = subprocess.check_output(bash_check).decode('utf-8')
    except OSError:
        raise RuntimeError('bowtie-build not found; check if it is installed and in $PATH\n')

    local_version = version_call.split('\n')[0].split()[-1]
    assert version == local_version, 'bowtie-build version incompatibility %s != %s' % (version, local_version)

    # Change file formats for winbash option
    if bash_sub == 'winbash':
        fasta = windowsToBash(fasta, 'winbash')
        refpath = windowsToBash(refpath, 'winbash')

    # Make sure spaces in files are escaped for bash
    if ' ' in fasta:
        fasta = '"' + fasta + '"'

    if ' ' in refpath:
        refpath = '"' + refpath + '"'

    # Initialize error/output file
    err_handle = open(err_file, "w")

    # Call bowtie-build
    if not bash_sub:
        subprocess.call(['bowtie-build', '-f', fasta, refpath],
                        stdout=err_handle,
                        stderr=subprocess.STDOUT)
    else:
        bash_call = ['%s/bash' % bash_dir, '--login', '-c', 'bowtie-build -f %s %s' % (fasta, refpath)]
        subprocess.call(bash_call,
                        stdout=err_handle,
                        stderr=subprocess.STDOUT)

    err_handle.close()


def bowtie(fasta, refpath, samfile, flags, err_file='bowtie.stderr',
           version=None,
           bash_sub=False, bash_dir='/'):
    """
    Call bowtie2-align
    :param fasta: fasta file for alignment
    :param refpath: Path to bowtie index files (.bt2)
    :param samfile: path to samfile output
    :param flags: A list of bowtie flags such as ['-f', '-N', '0', '-L', '12', '-a']
    :param err_file: File path to where output is written
    :param version: Enforces bowtie version number
    :param bash_sub: False/cygwin/winbash
    :param bash_dir: Location of bash.exe, if bash_sub not False
    :return:
    """

    # Check that we have access to bowtie
    try:
        if not bash_sub:
            version_call = subprocess.check_output(['bowtie', '--version']).decode('utf-8')
        else:
            bash_check = ['%s/bash' % bash_dir, '--login', '-c', 'bowtie --version']
            version_call = subprocess.check_output(bash_check).decode('utf-8')
    except OSError:
        raise RuntimeError('bowtie not found; check if it is installed and in $PATH\n')

    # Check that version is the expected version
    local_version = version_call.split('\n')[0].split()[-1]
    assert version == local_version, 'bowtie version incompatibility %s != %s' % (version, local_version)

    # Change file formats for winbash option
    if bash_sub == 'winbash':
        fasta = windowsToBash(fasta, 'winbash')
        refpath = windowsToBash(refpath, 'winbash')
        samfile = windowsToBash(samfile, 'winbash')

    # Make sure spaces in files are escaped for bash
    if ' ' in fasta:
        fasta = '"' + fasta + '"'

    if ' ' in refpath:
        refpath = '"' + refpath + '"'

    if ' ' in samfile:
        samfile = '"' + samfile + '"'

    # Initialize error/output file
    err_handle = open(err_file, "w")

    # Call bowtie2
    bowtie_args = ['bowtie', refpath, fasta, '-S', samfile] + list(flags)

    if not bash_sub:
        subprocess.call(bowtie_args,
                        stdout=err_handle,
                        stderr=subprocess.STDOUT)
    else:
        bash_call = ['%s/bash' % bash_dir, '--login', '-c'] + [' '.join(bowtie_args)]
        subprocess.call(bash_call,
                        stdout=err_handle,
                        stderr=subprocess.STDOUT)

    err_handle.close()


def bowtie2_build(fasta, refpath, err_file='bowtie2-build.stderr',
                  version=None,
                  bash_sub=False, bash_dir='/',):
    """
    Construct bowtie2 indices (.bt2 files)
    :param fasta: Path to FASTA containing reference sequences
    :param refpath: Path to index base (e.g. C:/directory/e_coli_MG1655)
    :param err_file: File path to where output is written
    :param version: Enforces bowtie2 version number (e.g., '2.2.1')
    :param bash_sub: False/cygwin/winbash
    :param bash_dir: Location of bash.exe, if bash_sub not False
    :return:
    """

    try:
        if not bash_sub:
            version_call = subprocess.check_output(['bowtie2-build', '--version']).decode('utf-8')
        else:
            bash_check = ['%s/bash' % bash_dir, '--login', '-c', 'bowtie2-build --version']
            version_call = subprocess.check_output(bash_check).decode('utf-8')
    except OSError:
        raise RuntimeError('bowtie2-build not found; check if it is installed and in $PATH\n')

    local_version = version_call.split('\n')[0].split()[-1]
    assert version == local_version, 'bowtie2-build version incompatibility %s != %s' % (version, local_version)

    # Change file formats for winbash option
    if bash_sub == 'winbash':
        fasta = windowsToBash(fasta, 'winbash')
        refpath = windowsToBash(refpath, 'winbash')

    # Make sure spaces in files are escaped for bash
    if ' ' in fasta:
        fasta = '"' + fasta + '"'

    if ' ' in refpath:
        refpath = '"' + refpath + '"'

    # Initialize error/output file
    err_handle = open(err_file, "w")

    # Call bowtie2-build
    if not bash_sub:
        subprocess.call(['bowtie2-build', '-f', fasta, refpath],
                        stdout=err_handle,
                        stderr=subprocess.STDOUT)
    else:
        bash_call = ['%s/bash' % bash_dir, '--login', '-c', 'bowtie2-build -f %s %s' % (fasta, refpath)]
        subprocess.call(bash_call,
                        stdout=err_handle,
                        stderr=subprocess.STDOUT)

    err_handle.close()


def bowtie2(fasta, refpath, samfile, flags, err_file='bowtie2.stderr',
            version=None,
            bash_sub=False, bash_dir='/'):
    """
    Call bowtie2-align
    :param fasta: fasta file for alignment
    :param refpath: Path to bowtie2 index files (.bt2)
    :param samfile: path to samfile output
    :param flags: A list of bowtie2 flags such as ['-f', '-N', '0', '-L', '12', '-a']
    :param err_file: File path to where output is written
    :param version: Enforces bowtie2 version number
    :param bash_sub: False/cygwin/winbash
    :param bash_dir: Location of bash.exe, if bash_sub not False
    :return:
    """

    # Check that we have access to bowtie2
    try:
        if not bash_sub:
            version_call = subprocess.check_output(['bowtie2', '--version']).decode('utf-8')
        else:
            bash_check = ['%s/bash' % bash_dir, '--login', '-c', 'bowtie2 --version']
            version_call = subprocess.check_output(bash_check).decode('utf-8')
    except OSError:
        raise RuntimeError('bowtie2 not found; check if it is installed and in $PATH\n')

    # Check that version is the expected version
    local_version = version_call.split('\n')[0].split()[-1]
    assert version == local_version, 'bowtie2 version incompatibility %s != %s' % (version, local_version)

    # Change file formats for winbash option
    if bash_sub == 'winbash':
        fasta = windowsToBash(fasta, 'winbash')
        refpath = windowsToBash(refpath, 'winbash')
        samfile = windowsToBash(samfile, 'winbash')

    # Make sure spaces in files are escaped for bash
    if ' ' in fasta:
        fasta = '"' + fasta + '"'

    if ' ' in refpath:
        refpath = '"' + refpath + '"'

    if ' ' in samfile:
        samfile = '"' + samfile + '"'

    # Initialize error/output file
    err_handle = open(err_file, "w")

    # Call bowtie2
    bowtie_args = ['bowtie2', '-x', refpath, '-U', fasta, '-S', samfile] + list(flags)

    if not bash_sub:
        subprocess.call(bowtie_args,
                        stdout=err_handle,
                        stderr=subprocess.STDOUT)
    else:
        bash_call = ['%s/bash' % bash_dir, '--login', '-c'] + [' '.join(bowtie_args)]
        subprocess.call(bash_call,
                        stdout=err_handle,
                        stderr=subprocess.STDOUT)

    err_handle.close()


def samtools(function, infile, outfile=None, flags=(), err_file='samtools.stderr',
             version=None,
             bash_sub=False, bash_dir='/'):
    """
    Call samtools functions view, sort, or index
    :param function: view, sort, index
    :param infile: file path to input file
    :param outfile: file path to output file if needed (e.g. outfile=None for index function)
    :param flags: list of additional flags that can be passed to the given function
    :param err_file: samtools output/error file
    :param version: Enforces samtools version number
    :param bash_sub: False/cygwin/winbash
    :param bash_dir: Location of bash.exe, if bash_sub not False
    :return:
    """

    # Check that we have access to samtools
    try:
        if not bash_sub:
            version_call = subprocess.check_output(['samtools', '--version']).decode('utf-8')
        else:
            bash_check = ['%s/bash' % bash_dir, '--login', '-c', 'samtools --version']
            version_call = subprocess.check_output(bash_check).decode('utf-8')
    except OSError:
        raise RuntimeError('samtools not found; check if it is installed and in $PATH\n')

    # Check that version is the expected version
    local_version = version_call.split(' ')[1].split('\n')[0]
    assert version == local_version, 'samtools version incompatibility %s != %s' % (version, local_version)

    # Check that chosen function is supported
    if function not in ['view', 'sort', 'index']:
        raise Exception('Function %s not currently supported by python samtools wrapper' % function)

    # Change file formats for winbash option
    if bash_sub == 'winbash':
        infile = windowsToBash(infile, 'winbash')

        if outfile:
            outfile = windowsToBash(outfile, 'winbash')

    # Make sure spaces in files are escaped for bash
    if ' ' in infile:
        infile = '"' + infile + '"'

    if outfile and ' ' in outfile:
        outfile = '"' + outfile + '"'

    # Initialize arguments
    samtools_args = ['samtools', function, infile] + list(flags)

    if outfile:
        samtools_args += ['-o', outfile]

    # Initialize error file
    err_handle = open(err_file, "w")

    # Call samtools
    if not bash_sub:
        subprocess.call(samtools_args,
                        stdout=err_handle,
                        stderr=subprocess.STDOUT)
    else:
        bash_call = ['%s/bash' % bash_dir, '--login', '-c'] + [' '.join(samtools_args)]
        subprocess.call(bash_call,
                        stdout=err_handle,
                        stderr=subprocess.STDOUT)

    err_handle.close()


def bedtools(function, infiles, outfile=None, flags=(), err_file='bedtools.stderr',
             version=None,
             bash_sub=False, bash_dir='/'):
    """
    Call bedtools functions window, getfasta
        window arguments: infiles = [bed/bam file, genome gff file]
                          outfile = output bedfile
                          flags = standard flags are ['-l', '#', '-r', '0', '-sw', '-sm', '-bed']: '-l' and '-r' and the
                                  numbers following them specify the left and right side extension of the window from
                                  either edge of the feature. '-l' corresponds to the start site and '-r' corresponds
                                  to the end of the feature when '-sw' is specified, which tells bedtools window to
                                  define left and right based on strand. '-sm' specifies that overlaps will only be
                                  reported when on the same strand. '-bed' specifies bedfile output.
        getfasta arguments: infiles = [genome fasta, bedfile with genome coordinates]
                            outfile = output fasta file
                            flags = standard flags are ['-name', '-s']: '-name' specifies that the output fasta
                                    identifier will be pulled from the name column (column 3) of the bedfile, '-s'
                                    specifies that strandedness will be respected in extracting fasta sequences.
    :param function: str of function name
    :param infiles: list of str for file paths to input files
    :param outfile: str file path to output file if needed
    :param flags: list of additional flags that can be passed to the given function
    :param err_file: bedtools output/error file
    :param version: Enforces bedtools version number
    :param bash_sub: False/cygwin/winbash
    :param bash_dir: Location of bash.exe, if bash_sub not False
    :return:
    """

    # Check access to samtools
    try:
        if not bash_sub:
            version_call = subprocess.check_output(['bedtools', '--version']).decode('utf-8')
        else:
            bash_check = ['%s/bash' % bash_dir, '--login', '-c', 'bedtools --version']
            version_call = subprocess.check_output(bash_check).decode('utf-8')
    except OSError:
        raise RuntimeError('bedtools not found; check if it is installed and in $PATH\n')

    # Check that version is the expected version
    local_version = version_call.split(' ')[1].split('\n')[0]
    assert version == local_version, 'bedtools version incompatibility %s != %s' % (version, local_version)

    # Check that chosen function is supported
    if function not in ['window', 'getfasta']:
        raise Exception('Function %s not currently supported by python bedtools wrapper' % function)

    # Change file formats for shell option
    if bash_sub:
        for n, infile in enumerate(infiles):
            infiles[n] = windowsToBash(infile, bash_sub)

    # Make sure spaces in files are escaped for bash
    for n, infile in enumerate(infiles):
        if ' ' in infile:
            infiles[n] = '"' + infile + '"'

    # Initialize error file
    err_handle = open(err_file, 'w')

    # Assemble arguments list, call function
    bedtools_args = ['bedtools', function]

    if function == 'window':
        # Check for correct input number files, and for file type of first input. Takes 1 bam/bed, 1 gff file
        if len(infiles) != 2:
            raise IndexError('bedtools window requires 2 input files')

        if infiles[0].replace('"', '').split('.')[-1] == 'bed':
            bedtools_args += ['-a', infiles[0]]
        elif infiles[0].replace('"', '').split('.')[-1] == 'bam':
            bedtools_args += ['-abam', infiles[0]]
        else:
            raise NameError('bedtools window takes only .bed or .bam files as input')

        bedtools_args += ['-b', infiles[1]] + list(flags)

        # Initialize output file
        out_handle = open(outfile, 'w')

        # Call bedtools
        if not bash_sub:
            subprocess.call(bedtools_args,
                            stdout=out_handle,
                            stderr=err_handle)
        else:
            bash_call = ['%s/bash' % bash_dir, '--login', '-c'] + [' '.join(bedtools_args)]
            subprocess.call(bash_call,
                            stdout=out_handle,
                            stderr=err_handle)

        out_handle.close()

    elif function == 'getfasta':
        # Check for correct input number files. Takes 1 fasta, 1 bed file
        # Standard flags for extractSeq function: -name -s
        if len(infiles) != 2:
            raise IndexError('bedtools getfasta requires 2 input files')

        bedtools_args += ['-fi', infiles[0], '-bed', infiles[1]] + list(flags)

        # Initialize output file
        if outfile:
            if bash_sub:
                outfile = windowsToBash(outfile, bash_sub)

            outfile = '"' + outfile + '"'
            bedtools_args += ['-fo', outfile]
        else:
            raise IndexError('getfasta requires an outfile')

        # Call bedtools
        if not bash_sub:
            subprocess.call(bedtools_args,
                            stderr=err_handle)
        else:
            bash_call = ['%s/bash' % bash_dir, '--login', '-c'] + [' '.join(bedtools_args)]
            subprocess.call(bash_call,
                            stderr=err_handle)

    else:
        print('%s function not currently supported' % function)

    err_handle.close()
