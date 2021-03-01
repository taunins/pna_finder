import datetime
from tkinter import *
import tkinter.filedialog
import functools
import os
from collections import OrderedDict
from . import pna_pipeline_bt
import pna_finder_dev


class pna_app:
    def __init__(self, master, bash_path='/', home_path='/', shell_type='cygwin'):
        """
        Presents user with choice of the two PNA Finder functions: Get Sequences and Find Off-Targets.
        Initializes tkinter frame and file search directories.
        (see manual.pdf in package root directory for more info)
        :param master:
        :param bash_path:
        :param home_path:
        """

        self.master = master
        self.frame_welcome = Frame(self.master)
        self.frame_welcome.pack()
        self.bash_path = bash_path
        self.home_path = home_path
        self.shell_type = shell_type

        Label(
            self.frame_welcome,
            text='Welcome to the PNA Finder Toolbox', font=('Arial', 12, 'bold')
        ).grid(columnspan=2, padx=5, pady=5)

        Label(
            self.frame_welcome,
            text='Select toolbox function:', font=('Arial', 11, 'italic')
        ).grid(rowspan=3, row=1, column=0, padx=5, pady=5)

        self.function = IntVar()
        self.function.set(0)

        Radiobutton(
            self.frame_welcome,
            text='Get Sequences',
            variable=self.function, value=0
        ).grid(row=1, column=1, sticky=W, padx=5, pady=0)

        Radiobutton(
            self.frame_welcome,
            text='Find Off-Targets',
            variable=self.function, value=1
        ).grid(row=2, column=1, sticky=W, padx=5, pady=0)

        Radiobutton(
            self.frame_welcome,
            text='Check Sequence Warnings',
            variable=self.function, value=2
        ).grid(row=3, column=1, sticky=W, padx=5, pady=0)

        Button(
            self.frame_welcome,
            text="OK", command=self.run_pnaFinder
        ).grid(row=10, column=0, sticky=W, padx=5, pady=10)

        Button(
            self.frame_welcome,
            text="QUIT", fg="red", command=self.frame_welcome.quit,
        ).grid(row=10, column=1, sticky=E, padx=5, pady=10)

    def run_pnaFinder(self):
        """
        Initializes dictionaries for PNA Finder GUI labels, parameters, and filenames. Runs the PNA Finder toolbox
        function that the user has selected.
        :return:
        """

        top = self.top = Toplevel(self.master)
        self.master.withdraw()

        self.labels = {}
        self.flabels = {}
        self.parameters = {}
        self.files = {}

        # JOB NAME BOX
        self.labels['job_name'] = Label(
            self.top,
            text='Job Name:', font=('Arial', 11, 'bold'))
        self.labels['job_name'].grid(row=0, column=0, sticky=W, padx=5, pady=5)

        self.job_name = StringVar()
        self.job_name.set('')
        Entry(self.top, width=20, textvariable=self.job_name).grid(row=0, columnspan=3, padx=20, pady=5)

        if self.function.get() == 0:
            top.title('PNA Finder - Get Sequences')
            self.function_name = 'get_sequences'
            self.app_gs()
        elif self.function.get() == 1:
            top.title('PNA Finder - Find Off-Targets')
            self.function_name = 'find_off_targets'
            self.app_ot()
        elif self.function.get() == 2:
            top.title('PNA Finder - Check Sequence Warnings')
            self.function_name = 'sequence_warnings'
            self.app_warnings()

    def app_gs(self):
        """
        Creates a dialog box that may be filled in with the Get Sequences function parameters and file selections
            get_sequences Inputs:
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

        # PNA SEQUENCE PARAMETERS
        self.labels['start'] = Label(  # Prompt for window start index
            self.top,
            text='Sequence Window Start:')
        self.labels['start'].grid(row=1, column=0, sticky=W, padx=5, pady=5)

        self.labels['end'] = Label(  # Prompt for window end index
            self.top,
            text='Sequence Window End:')
        self.labels['end'].grid(row=2, column=0, sticky=W, padx=5, pady=5)

        self.labels['length'] = Label(  # Prompt for PNA length
            self.top,
            text='PNA Length:')
        self.labels['length'].grid(row=3, column=0, sticky=W, padx=5, pady=5)

        self.labels['feature_type'] = Label(  # Prompt for PNA length
            self.top,
            text='Annotation Record Types:')
        self.labels['feature_type'].grid(row=4, column=0, sticky=W, padx=5, pady=5)

        self.parameters['start'] = StringVar()
        self.parameters['start'].set('-5')

        self.window_start = Entry(self.top, width=5, textvariable=self.parameters['start'])
        self.window_start.grid(row=1, column=2, sticky=E, padx=5, pady=5)

        self.parameters['end'] = StringVar()
        self.parameters['end'].set('-5')

        self.window_end = Entry(self.top, width=5, textvariable=self.parameters['end'])
        self.window_end.grid(row=2, column=2, sticky=E, padx=5, pady=5)

        self.parameters['length'] = StringVar()
        self.parameters['length'].set('12')

        self.pna_length = Entry(self.top, width=5, textvariable=self.parameters['length'])
        self.pna_length.grid(row=3, column=2, sticky=E, padx=5, pady=5)

        self.parameters['feature_type'] = StringVar()
        self.parameters['feature_type'].set('CDS')

        self.pna_length = Entry(self.top, width=10, textvariable=self.parameters['feature_type'])
        self.pna_length.grid(row=4, column=2, sticky=E, padx=5, pady=5)

        # ID LIST FILE SELECTION
        self.labels['id_list'] = Label(
            self.top,
            font=('Arial', 9, 'bold'),
            text='Select Gene ID List'
        )
        self.labels['id_list'].grid(row=7, column=0, sticky=W, padx=5, pady=5)

        id_list_button = Button(  # Button to select gene ID list file
            self.top,
            text="Select",
            command=functools.partial(self.browse, self.top, 'id_list',
                                      title='Select Gene ID List',
                                      filetypes=(("text files", "*.txt"), ("all files", "*.*")),
                                      row=7, column=1))
        id_list_button.grid(row=7, column=2, sticky=E, padx=5, pady=5)

        # GENOME ASSEMBLY FILE SELECTION
        self.labels['assembly'] = Label(
            self.top,
            font=('Arial', 9, 'bold'),
            text='Select Genome Assembly File'
        )
        self.labels['assembly'].grid(row=8, column=0, sticky=W, padx=5, pady=5)

        assembly_button = Button(  # Button to select genome assembly file
            self.top,
            text="Select",
            command=functools.partial(self.browse, self.top, 'assembly',
                                      title='Select genome assembly FASTA File',
                                      filetypes=(("fasta files", "*.fa *.fna"), ("all files", "*.*")),
                                      row=8, column=1))
        assembly_button.grid(row=8, column=2, sticky=E, padx=5, pady=5)

        # ANNOTATION FILE SELECTION
        Label(
            self.top,
            text='Select Annotation File Option:', font=("Arial", 11, "italic")
        ).grid(row=9, rowspan=2, sticky=N, pady=(20, 0))

        self.gff_option = IntVar()
        self.gff_option.set(0)

        Radiobutton(
            self.top,
            text='Build new gffutils database from GFF file',
            variable=self.gff_option, value=0
        ).grid(row=9, column=1, columnspan=2, sticky=W, padx=5, pady=(10, 0))

        Radiobutton(
            self.top,
            text='Select existing gffutils database file',
            variable=self.gff_option, value=1
        ).grid(row=10, column=1, columnspan=2, sticky=W, padx=5, pady=0)

        self.labels['gff'] = Label(
            self.top,
            font=('Arial', 9, 'bold'),
            text='Select Annotation File'
        )
        self.labels['gff'].grid(row=11, column=0, sticky=W, padx=5, pady=5)

        gff_button = Button(
            self.top,
            text="Select",
            command=functools.partial(self.browse, self.top, 'gff',
                                      title='Select Annotation File',
                                      row=11, column=1))
        gff_button.grid(row=11, column=2, sticky=E, padx=5, pady=5)

        # OUTPUT DIRECTORY SELECTION
        self.labels['output'] = Label(
            self.top,
            font=('Arial', 9, 'bold'),
            text='Select Output Directory'
        )
        self.labels['output'].grid(row=12, column=0, sticky=W, padx=5, pady=5)

        output_button = Button(
            self.top,
            text="Select",
            command=functools.partial(self.browse_folder, self.top, 'output',
                                      title='Select Output Directory',
                                      row=12, column=1))
        output_button.grid(row=12, column=2, sticky=E, padx=5, pady=5)

        # STRING ANALYSIS, MELT TEMPERATURE, AND SOLUBILITY WARNINGS OPTIONS
        self.STRING = IntVar()
        self.STRING.set(0)

        self.parameters['string_id'] = StringVar()
        self.parameters['string_id'].set('')

        self.warnings = IntVar()
        self.warnings.set(1)

        self.temp_option = IntVar()
        self.temp_option.set(1)

        self.labels['string_id'] = Label(
            self.top,
            text='Species NCBI Taxonomy ID', font=("Arial", 8, "bold"))
        self.labels['string_id'].grid(row=13, column=1, columnspan=2, sticky=E, pady=(10, 0))

        Checkbutton(
            self.top,
            text='Run STRING analysis',
            variable=self.STRING
        ).grid(row=14, column=0, sticky=W, padx=5, pady=(0, 5))

        self.string_id = Entry(self.top, width=5, textvariable=self.parameters['string_id'])
        self.string_id.grid(row=14, column=2, sticky=E, padx=5, pady=0)

        Checkbutton(
            self.top,
            text='Run sequence warnings analysis',
            variable=self.warnings
        ).grid(row=15, column=0, sticky=W, padx=5, pady=5)

        Checkbutton(
            self.top,
            text='Calculate single-stranded PNA/DNA melt temperature',
            variable=self.temp_option
        ).grid(row=16, column=0, columnspan=3, sticky=W, padx=5, pady=5)

        # JOB SUBMISSION
        Button(
            self.top,
            text="SUBMIT", command=self.submit
        ).grid(row=17, column=0, sticky=W, padx=5, pady=10)

        # WINDOW QUIT
        Button(
            self.top,
            text="QUIT", fg="red", command=self.quit_app,
        ).grid(row=17, column=2, sticky=E, padx=5, pady=10)

    def app_ot(self):
        """
        Creates a dialog box that may be filled in with the Find Off Targets function parameters and file selections
            find_off_targets Inputs:
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
                    homology    -   int (optional)
                    remove_files -  int (optional)
        (see manual.pdf in package root directory for more info)
        :return:
        """

        # OFF TARGET SEARCH PARAMETERS
        self.labels['start'] = Label(  # Prompt for window start index
            self.top,
            text='Off-Target Window Start:')
        self.labels['start'].grid(row=1, column=0, sticky=W, padx=5, pady=5)

        self.labels['end'] = Label(  # Prompt for window end index
            self.top,
            text='Off-Target Window End:')
        self.labels['end'].grid(row=2, column=0, sticky=W, padx=5, pady=5)

        self.labels['mismatch'] = Label(  # Prompt for number of mismatches allowed
            self.top,
            text='Maximum Mismatches:')
        self.labels['mismatch'].grid(row=3, column=0, sticky=W, padx=5, pady=5)

        self.labels['feature_type'] = Label(  # Prompt for GFF record type
            self.top,
            text='Annotation Record Types:')
        self.labels['feature_type'].grid(row=5, column=0, sticky=W, padx=5, pady=5)

        self.parameters['start'] = StringVar()
        self.parameters['start'].set('-20')

        self.window_start = Entry(self.top, width=5, textvariable=self.parameters['start'])
        self.window_start.grid(row=1, column=2, sticky=E, padx=5, pady=5)

        self.parameters['end'] = StringVar()
        self.parameters['end'].set('20')

        self.window_end = Entry(self.top, width=5, textvariable=self.parameters['end'])
        self.window_end.grid(row=2, column=2, sticky=E, padx=5, pady=5)

        self.parameters['mismatch'] = StringVar()
        self.parameters['mismatch'].set('0')

        self.mismatch = Entry(self.top, width=5, textvariable=self.parameters['mismatch'])
        self.mismatch.grid(row=3, column=2, sticky=E, padx=5, pady=5)

        self.parameters['feature_type'] = StringVar()
        self.parameters['feature_type'].set('CDS')

        self.feature_type = Entry(self.top, width=10, textvariable=self.parameters['feature_type'])
        self.feature_type.grid(row=5, column=2, sticky=E, padx=5, pady=5)

        # PNA FASTA FILE SELECTION
        self.labels['targets'] = Label(
            self.top,
            font=('Arial', 9, 'bold'),
            text='Select PNA Targets File'
        )
        self.labels['targets'].grid(row=6, column=0, sticky=W, padx=5, pady=5)

        targets_button = Button(
            self.top,
            text="Select",
            command=functools.partial(self.browse, self.top, 'targets',
                                      title='Select PNA targets FASTA file',
                                      filetypes=(("fasta files", "*.fa *.fna"), ("all files", "*.*")),
                                      row=6, column=1))
        targets_button.grid(row=6, column=2, sticky=E, padx=5, pady=5)

        # GENOME ASSEMBLY FILE SELECTION
        Label(
            self.top,
            text='Select Index File Option:', font=("Arial", 11, "italic")
        ).grid(row=7, rowspan=2, sticky=N, pady=(20, 0))

        self.index_option = IntVar()
        self.index_option.set(0)

        Radiobutton(
            self.top,
            text='Build new index from FASTA file',
            variable=self.index_option, value=0
        ).grid(row=8, column=1, sticky=W, padx=5, pady=(5, 0))

        Radiobutton(
            self.top,
            text='Select existing index (any of the 6 files)',
            variable=self.index_option, value=1
        ).grid(row=9, column=1, sticky=W, padx=5, pady=(0, 5))

        self.labels['assembly'] = Label(
            self.top,
            font=('Arial', 9, 'bold'),
            text='Select Assembly/Index File'
        )
        self.labels['assembly'].grid(row=10, column=0, sticky=W, padx=5, pady=5)

        assembly_button = Button(
            self.top,
            text="Select",
            command=functools.partial(self.browse, self.top, 'assembly',
                                      title='Select Assembly/Index File',
                                      row=10, column=1))
        assembly_button.grid(row=10, column=2, sticky=E, padx=5, pady=5)

        # ANNOTATION FILE SELECTION
        self.labels['gff'] = Label(
            self.top,
            font=('Arial', 9, 'bold'),
            text='Select Annotation File'
        )
        self.labels['gff'].grid(row=11, column=0, sticky=W, padx=5, pady=5)

        gff_button = Button(
            self.top,
            text="Select",
            command=functools.partial(self.browse, self.top, 'gff',
                                      title='Select Annotation File',
                                      filetypes=(("gff files", "*.gff *.gtf"), ("all files", "*.*")),
                                      row=11, column=1))
        gff_button.grid(row=11, column=2, sticky=E, padx=5, pady=5)

        # OUTPUT DIRECTORY SELECTION
        self.labels['output'] = Label(
            self.top,
            font=('Arial', 9, 'bold'),
            text='Select Output Directory'
        )
        self.labels['output'].grid(row=12, column=0, sticky=W, padx=5, pady=5)

        output_button = Button(
            self.top,
            text="Select",
            command=functools.partial(self.browse_folder, self.top, 'output',
                                      title='Select Output Directory',
                                      row=12, column=1))
        output_button.grid(row=12, column=2, sticky=E, padx=5, pady=5)

        # FILE OUTPUT OPTIONS
        self.count = IntVar()
        self.count.set(1)

        self.homology = IntVar()
        self.homology.set(0)

        self.remove_files = IntVar()
        self.remove_files.set(0)

        Checkbutton(
            self.top,
            text='Create off-target count file',
            variable=self.count
        ).grid(row=13, column=0, sticky=W, padx=5, pady=5)

        Checkbutton(
            self.top,
            text='Check for potential feature homology',
            variable=self.homology
        ).grid(row=14, column=0, columnspan=2, sticky=W, padx=5, pady=5)

        Checkbutton(
            self.top,
            text='Remove intermediate and error files after run',
            variable=self.remove_files
        ).grid(row=15, column=0, columnspan=2, sticky=W, padx=5, pady=5)

        # JOB SUBMISSION
        Button(
            self.top,
            text="SUBMIT", command=self.submit
        ).grid(row=16, column=0, sticky=W, padx=5, pady=10)

        # WINDOW QUIT
        Button(
            self.top,
            text="QUIT", fg="red", command=self.quit_app,
        ).grid(row=16, column=2, sticky=E, padx=5, pady=10)

    def app_warnings(self):
        """
        Creates a dialog box that may be filled in with the Check Warnings file selections
            sequence_warnings Inputs:
                self.files dictionary
                    pna_list    -   str
                    output      -   str
        (see manual.pdf in package root directory for more info)
        :return:
        :return:
        """

        # PNA FILE SELECTION
        Label(
            self.top,
            text='Select PNA Upload Filetype:', font=("Arial", 11, "italic")
        ).grid(row=1, rowspan=2, sticky=N, pady=(20, 0))

        self.pna_option = IntVar()
        self.pna_option.set(0)

        Radiobutton(
            self.top,
            text='Upload PNA target FASTA file',
            variable=self.pna_option, value=0
        ).grid(row=1, column=1, columnspan=2, sticky=W, padx=5, pady=(10, 0))

        Radiobutton(
            self.top,
            text='Upload PNA sequences table',
            variable=self.pna_option, value=1
        ).grid(row=2, column=1, columnspan=2, sticky=W, padx=5, pady=0)

        self.labels['pna_list'] = Label(
            self.top,
            font=('Arial', 9, 'bold'),
            text='Select PNA File'
        )
        self.labels['pna_list'].grid(row=3, column=0, sticky=W, padx=5, pady=5)

        pna_button = Button(
            self.top,
            text="Select",
            command=functools.partial(self.browse, self.top, 'pna_list',
                                      title='Select PNA File',
                                      row=3, column=1))
        pna_button.grid(row=3, column=2, sticky=E, padx=5, pady=5)

        # OUTPUT DIRECTORY SELECTION
        self.labels['output'] = Label(
            self.top,
            font=('Arial', 9, 'bold'),
            text='Select Output Directory'
        )
        self.labels['output'].grid(row=4, column=0, sticky=W, padx=5, pady=5)

        output_button = Button(
            self.top,
            text="Select",
            command=functools.partial(self.browse_folder, self.top, 'output',
                                      title='Select Output Directory',
                                      row=4, column=1))
        output_button.grid(row=4, column=2, sticky=E, padx=5, pady=5)

        # MELT TEMPERATURE CALCULATION OPTION
        self.temp_option = IntVar()
        self.temp_option.set(1)

        Checkbutton(
            self.top,
            text='Calculate single-stranded PNA/DNA melt temperature',
            variable=self.temp_option
        ).grid(row=5, column=0, columnspan=3, sticky=W, padx=5, pady=5)

        # JOB SUBMISSION
        Button(
            self.top,
            text="SUBMIT", command=self.submit
        ).grid(row=6, column=0, sticky=W, padx=5, pady=10)

        # WINDOW QUIT
        Button(
            self.top,
            text="QUIT", fg="red", command=self.quit_app,
        ).grid(row=6, column=2, sticky=E, padx=5, pady=10)

    def browse(self, frame, filekey, title='Select file', filetypes=(("all files", "*.*"),),
               row=0, column=0):
        """
        Function prompts user to select a file for input into the PNA Finder toolbox, and records choice in GUI display
        :param frame: Tkinter active frame
        :param filekey: Key for the self.files dictionary, corresponding to the file being selected
        :param title: File browser window title
        :param filetypes: File browser possible filetype selections
        :param row: Row location of file label within the Tkinter frame
        :param column: Column location of file label within the Tkinter frame
        :return:
        """

        try:
            if self.flabels['%s_display' % filekey].get():
                label_pass = True
            else:
                label_pass = False
        except KeyError:
            self.flabels['%s_display' % filekey] = StringVar()
            label_pass = False

        if filekey == 'gff' and self.function.get() == 0:
            if self.gff_option.get() == 0:
                filetypes = (("gff/gtf files", "*.gff *gtf"), ("all files", "*.*"))
            elif self.gff_option.get() == 1:
                filetypes = (("db files", "*.db"), ("all files", "*.*"))
            else:
                raise ValueError('Something has gone horribly wrong!')

        elif filekey == 'assembly' and self.function.get() == 1:
            if not self.index_option.get():
                filetypes = (("fasta files", "*.fa *.fna"), ("all files", "*.*"))
                title = 'Select genome assembly FASTA file'
            elif self.index_option.get():
                filetypes = (("Bowtie index files", "*.ebwt"), ("all files", "*.*"))
                title = 'Select genome assembly BT index file'
            else:
                raise ValueError('Something has gone horribly wrong!')

        elif filekey == 'pna_list':
            if not self.pna_option.get():
                filetypes = (("fasta files", "*.fa *.fna"), ("all files", "*.*"))
                title = 'Select PNA target FASTA file'
            if self.pna_option.get():
                filetypes = (("text files", "*.txt"), ("all files", "*.*"))
                title = 'Select PNA sequence table file'

        f_out = str(tkinter.filedialog.askopenfilename(initialdir=self.home_path, title=title,
                                                       filetypes=filetypes))  # TODO: try askopenfilenames

        if f_out:
            self.files[filekey] = f_out
            file_loc_list = self.files[filekey].split('/')
            file_type = file_loc_list[-1].split('.')[-1]

            if len(file_loc_list[-1]) > 15:
                if len(file_loc_list) > 2:
                    if len(file_loc_list[-1].rsplit('.', 1)[0]) > 15:
                        temp_label = file_loc_list[0] + '/.../' + file_loc_list[-1].rsplit('.', 1)[0][
                                                                  0:15] + '... .' + file_type
                    else:
                        temp_label = file_loc_list[0] + '/.../' + file_loc_list[-1]
                    self.flabels['%s_display' % filekey].set(temp_label)
                else:
                    temp_label = file_loc_list[0] + '/' + file_loc_list[-1][0:15] + '... .' + file_type
                    self.flabels['%s_display' % filekey].set(temp_label)
            else:
                if len(file_loc_list) > 2:
                    temp_label = file_loc_list[0] + '/.../' + file_loc_list[-1]
                    self.flabels['%s_display' % filekey].set(temp_label)
                else:
                    self.flabels['%s_display' % filekey].set(self.files[filekey])

            if label_pass:
                pass
            else:
                Label(
                    frame,
                    textvariable=self.flabels['%s_display' % filekey]
                ).grid(row=row, column=column, sticky=E, padx=5, pady=5)

        return

    def browse_folder(self, frame, filekey, title='Select folder', row=0, column=0):
        """
        Function prompts user to select a folder for input into the PNA Finder toolbox, and records choice in GUI display
        :param frame: Tkinter active frame
        :param filekey: Key for the self.files dictionary, corresponding to the folder being selected
        :param title: File browser window title
        :param row: Row location of file label within the Tkinter frame
        :param column: Column location of file label within the Tkinter frame
        :return:
        """

        try:
            if self.flabels['%s_display' % filekey]:
                label_pass = True
            else:
                label_pass = False
        except KeyError:
            self.flabels['%s_display' % filekey] = StringVar()
            label_pass = False

        d_out = str(tkinter.filedialog.askdirectory(initialdir=self.home_path, title=title))

        if d_out:
            self.files[filekey] = d_out
            file_loc_list = self.files[filekey].split('/')

            if len(file_loc_list[-1]) > 15:
                if len(file_loc_list) > 2:
                    self.flabels['%s_display' % filekey].set(file_loc_list[0] + '/.../' + file_loc_list[-1][0:15])
                else:
                    self.flabels['%s_display' % filekey].set(file_loc_list[0] + '/' + file_loc_list[-1][0:15])
            else:
                if len(file_loc_list) > 2:
                    self.flabels['%s_display' % filekey].set(file_loc_list[0] + '/.../' + file_loc_list[-1])
                else:
                    self.flabels['%s_display' % filekey].set(self.files[filekey])

            if label_pass:
                pass
            else:
                Label(
                    frame,
                    textvariable=self.flabels['%s_display' % filekey]
                ).grid(row=row, column=column, sticky=E, padx=5, pady=5)

        return

    def check_inputs(self):
        """
        Function checks the inputs for PNA tools to determine whether the necessary inputs are submitted in the
        correct format
            get_sequences Inputs:
                self.parameters dictionary:
                    start       -   int
                    end         -   int
                    pna_length  -   int
                    feature_type -   str
                    string_id   -   int (optional)
                self.files dictionary
                    id_list     -   str
                    assembly    -   str
                    gff         -   str
                    output      -   str
                self.STRING     -   int (optional)
                self.warnings   -   int (optional)

            find_off_targets Inputs:
                self.parameters dictionary:
                    start       -   int
                    end         -   int
                    mismatch    -   int
                    feature_type -  str
                self.files dictionary
                    targets     -   str
                    assembly    -   str
                    gff         -   str
                    output      -   str
                self.count      -   int (optional)
                self.homology   -   int (optional)
                self.remove_files - int (optional)

            sequence_warnings Inputs
                pna_option  -   int
                pna_list    -   str
                temp_option -   int (optional)
        """

        check = True
        check_warnings = OrderedDict()
        check_messages = OrderedDict()
        submit_dict = OrderedDict()
        p_dict = f_dict = None

        if self.job_name.get() == '':
            now = datetime.datetime.now()
            self.job_name.set('%s_%s' % (str(now.date()), self.function_name))
        else:
            for character in '<>:"/\|?*':
                if character in self.job_name.get():
                    check = False
                    check_warnings['job_name'] = r'Job name must be valid directory name (Remove characters: <>:"/\|?*)'

        if self.function.get() == 0:
            p_dict = OrderedDict([('start', 'Sequence Window Start'),
                                  ('end', 'Sequence Window End'),
                                  ('length', 'PNA Length'),
                                  ('feature_type', 'Annotation Record Types'),
                                  ('string_id', 'STRING Species ID')])

            f_dict = OrderedDict([('id_list', 'Gene ID List'),
                                  ('assembly', 'Genome Assembly File'),
                                  ('gff', 'Annotation File'),
                                  ('output', 'Output Directory')])

            submit_dict['gff_option'] = self.gff_option.get()
            submit_dict['STRING'] = self.STRING.get()
            submit_dict['warnings'] = self.warnings.get()
            submit_dict['temp_option'] = self.temp_option.get()

        elif self.function.get() == 1:
            p_dict = OrderedDict([('start', 'Off-Target Window Start'),
                                  ('end', 'Off-Target Window End'),
                                  ('mismatch', 'Maximum Mismatches'),
                                  ('feature_type', 'Annotation Record Types')])

            f_dict = OrderedDict([('targets', 'PNA Targets File'),
                                  ('assembly', 'Assembly/Index File'),
                                  ('gff', 'Annotation File'),
                                  ('output', 'Output Directory')])

            submit_dict['index_option'] = self.index_option.get()
            submit_dict['count'] = self.count.get()
            submit_dict['homology'] = self.homology.get()
            submit_dict['remove_files'] = self.remove_files.get()

        elif self.function.get() == 2:
            p_dict = OrderedDict([])
            f_dict = OrderedDict([('pna_list', 'PNA Sequences File'),
                                  ('output', 'Output Directory')])

            submit_dict['pna_option'] = self.pna_option.get()
            submit_dict['temp_option'] = self.temp_option.get()

        for p_key in p_dict:
            try:  # Entry data check
                submit_dict[p_key] = int(self.parameters[p_key].get().strip())
                if self.function.get() == 0 and p_key == 'string_id' and not self.STRING.get():
                    check_messages[p_key] = 'STRING ID entered but STRING analysis option not selected'

                elif p_key == 'end' and submit_dict[p_key] < submit_dict['start']:
                    check = False
                    check_warnings[p_key] = 'Window Start must be less than or equal to Window End'

            except ValueError:
                if self.function.get() == 0:
                    if p_key == 'string_id' and not self.STRING.get():
                        continue
                elif self.function.get() == 1:
                    if p_key == 'end' and self.parameters[p_key].get() == '':
                        submit_dict[p_key] = self.parameters[p_key].get()
                        continue
                    elif p_key == 'mismatch':
                        try:
                            mismatch_list = []
                            for mismatch in self.parameters[p_key].get().split(','):
                                mismatch_list += [int(mismatch.strip())]
                            submit_dict[p_key] = mismatch_list
                            continue
                        except ValueError:
                            check = False
                            check_warnings[p_key] = 'Mismatch must be integer or comma-separated list of integers'
                            continue

                if p_key == 'feature_type':
                    submit_dict[p_key] = self.parameters[p_key].get().strip()
                    continue

                check = False
                if self.parameters[p_key].get() == '':
                    check_warnings[p_key] = 'No value entered for %s' % p_dict[p_key]
                else:
                    check_warnings[p_key] = 'Must enter integer value for %s' % p_dict[p_key]

        for f_key in f_dict:
            try:
                if os.path.exists(self.files[f_key]):
                    if f_key == 'output':
                        submit_dict[f_key] = str(self.files[f_key].rstrip()) + '/' + str(self.job_name.get().rstrip())

                    elif self.function.get() == 0 and f_key == 'gff' and \
                            ((self.gff_option.get() and self.files[f_key].rsplit('.', 1)[1] != 'db') or
                             (not self.gff_option.get() and self.files[f_key].rsplit('.', 1)[1] != 'gff')):
                        check = False
                        check_warnings[f_key] = 'Wrong filetype selected for Annotation File Option'

                    elif self.function.get() and f_key == 'assembly' and \
                            ((self.index_option.get() and self.files[f_key].rsplit('.', 1)[1] != 'ebwt' or
                              (not self.index_option.get() and self.files[f_key].rsplit('.', 1)[1] not in ['fa', 'fna',
                                                                                                           'fasta']))):
                        check = False
                        check_warnings[f_key] = 'Wrong filetype selected for Index File Option'

                    else:
                        submit_dict[f_key] = str(self.files[f_key].rstrip())

                else:
                    check = False
                    check_warnings[f_key] = 'Filepath for %s does not exist' % f_dict[f_key]

            except KeyError:
                check = False
                check_warnings[f_key] = ('Must select %s' % f_dict[f_key])

        return check, submit_dict, check_warnings, check_messages

    def submit(self):
        """
        Calls the self.check_inputs() function and displays warning/error messages if input criteria are not satisfied.
        If input criteria are satisfied, the function submits the job to the external pna_pipeline.py script according
        to the variable self.function.
        :return:
        """

        check, submit_dict, check_warnings, check_messages = self.check_inputs()

        for l_key in self.labels:
            self.labels[l_key].config(fg='black')

        if not check:
            try:
                self.submit_warning.destroy()
            except AttributeError:
                pass

            warning_text = ''

            for w_key in check_warnings:
                self.labels[w_key].config(fg='red')
                warning_text += check_warnings[w_key] + '\n'

            self.submit_warning = Label(
                self.top,
                text=warning_text, fg='red')
            self.submit_warning.grid(row=20, columnspan=3, padx=5, pady=5)
            return

        else:
            submit_dict['home'] = self.home_path
            submit_dict['bash'] = self.bash_path
            submit_dict['shell_type'] = self.shell_type

            try:
                self.submit_warning.destroy()
            except AttributeError:
                pass

            if check_messages:
                message_text = ''

                for m_key in check_messages:
                    self.labels[m_key].config(fg='#FF6D00')
                    message_text += check_messages[m_key] + '\n'

                self.submit_warning = Label(
                    self.top,
                    text=message_text, fg='#FF6D00')
                self.submit_warning.grid(row=20, column=0, columnspan=3, padx=5, pady=5)

            start_time = datetime.datetime.now()
            print('Started at %s' % start_time)

            if self.function.get() == 0:
                pna_pipeline_bt.get_sequences(submit_dict)
            elif self.function.get() == 1:
                pna_pipeline_bt.find_off_targets(submit_dict)
            elif self.function.get() == 2:
                pna_pipeline_bt.sequence_warnings(submit_dict)

            end_time = datetime.datetime.now()
            print('Finished at %s' % end_time)
            print('Run execution time: %s' % (end_time - start_time))

    def quit_app(self):
        self.top.destroy()
        self.master.deiconify()


def startup():
    """
    Checks for proper start_info.txt file, to ensure preset home and bash directories.
    (see manual.pdf in package root directory for more info)
    :return:
    """

    bash_path = home_path = shell_type = False
    warning_text = 'The file "start_info.txt" has been deleted or moved.'

    try:
        start_handle = open('%s/start_info.txt' % pna_finder_dev.__path__[0], 'r')
        warning_text = 'The file "start_info.txt" is improperly built.'

        start_dict = {}
        for line in start_handle:
            path_dict = line.rstrip().split('\t')
            start_dict[path_dict[0]] = path_dict[1]

        home_path = start_dict['home']
        bash_path = start_dict['bash']
        shell_type = start_dict['shell']

    except IOError:
        start_warning = Tk()
        start_warning.title('PNA Finder Toolbox')

        Label(
            start_warning,
            text='ERROR', fg='red',
            font=('Arial', 12, 'bold')
        ).grid(row=0, padx=5, pady=5)

        Label(
            start_warning,
            text='%s\nRun start_info.py to rebuild this file.' % warning_text,
            font=('Arial', 10, 'italic')
        ).grid(row=1, padx=5, pady=5)

        Button(
            start_warning,
            text="QUIT", fg='red', command=start_warning.quit
        ).grid(row=2, padx=5, pady=10)

        start_warning.mainloop()
        start_warning.destroy()

    return bash_path, home_path, shell_type


def run():
    """
    Runs the PNA Finder toolbox GUI
    :return:
    """
    bash_path, home_path, shell_type = startup()

    if bash_path and home_path:
        root = Tk()
        root.title('PNA Finder Toolbox')

        pna_app(root, bash_path=bash_path, home_path=home_path, shell_type=shell_type)

        root.mainloop()
        root.destroy()
