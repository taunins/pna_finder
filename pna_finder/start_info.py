# Initial setup module for the PNA Finder toolbox

from tkinter import *
import tkinter.filedialog
import functools
import os


class start_app:
    def __init__(self, master):
        """
        Initializes the PNA Finder start_info app, in order to set home and bash directories.
        (see manual.pdf in package root directory for more info)
        :param master: tkinter GUI frame
        """
        self.master = master
        self.master.title('Choose toolbox startup settings')

        # Dictionaries for filenames and GUI label
        self.labels = {}
        self.flabels = {}
        self.files = {}

        # Variable initialization
        self.submit_warning = None

        # Header
        self.labels['header'] = Label(
            self.master,
            text=15*' ' + 'PNA Finder Startup Settings' + 15*' ',
            font=('Arial', 11, 'bold')
        )
        self.labels['header'].grid(row=0, columnspan=3)

        # Home directory
        self.labels['home'] = Label(
            self.master,
            font=('Arial', 9, 'bold'),
            text='Home directory'
        )
        self.labels['home'].grid(row=1, column=0, sticky=W, padx=5, pady=5)

        home_button = Button(
            self.master,
            text="Select",
            command=functools.partial(self.browse_folder, self.master, 'home',
                                      title='Select home directory',
                                      row=1, column=1))
        home_button.grid(row=1, column=2, sticky=E, padx=5, pady=5)

        # Bash directory
        self.labels['bash'] = Label(
            self.master,
            text='Bash shell directory',
            font=('Arial', 9, 'bold'),
        )
        self.labels['bash'].grid(row=2, column=0, sticky=W, padx=5, pady=5)

        bash_button = Button(
            self.master,
            text="Select",
            command=functools.partial(self.browse_folder, self.master, 'bash',
                                      title='Select Bash shell directory',
                                      row=2, column=1))
        bash_button.grid(row=2, column=2, sticky=E, padx=5, pady=5)

        # Bash shell type selection
        Label(
            self.master,
            text='Bash shell type:', font=("Arial", 9, "bold")
        ).grid(row=3, column=0, sticky=W, padx=5, pady=5)

        self.bash_option = IntVar()
        self.bash_option.set(0)

        Radiobutton(
            self.master,
            text='Cygwin',
            variable=self.bash_option, value=0
        ).grid(row=3, column=1, sticky=W, padx=5, pady=5)

        Radiobutton(
            self.master,
            text='Windows Bash',
            variable=self.bash_option, value=1
        ).grid(row=3, column=2, sticky=W, padx=5, pady=5)

        # OK and QUIT buttons
        Button(
            self.master,
            text="OK", command=self.writeStartInfo
        ).grid(row=10, column=0, sticky=W, padx=5, pady=10)

        Button(
            self.master,
            text="QUIT", fg="red", command=self.master.quit,
        ).grid(row=10, column=2, sticky=E, padx=5, pady=10)

    def browse_folder(self, frame, filekey, title='Select folder', row=0, column=0):
        """
        Function prompts user to select a folder for input into the start_info dialog, and records choice in GUI display
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

        d_out = str(tkinter.filedialog.askdirectory(initialdir='/', title=title))

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

    def writeStartInfo(self):
        """
        Function writes user-selected start folders to start_info.txt file
        :return:
        """
        if 'home' in list(self.files.keys()) and 'bash' in list(self.files.keys()):
            info_handle = open('%s/start_info.txt' % os.path.dirname(__file__), 'w')
            if self.bash_option.get() == 0:
                bash_option = 'cygwin'
            else:
                bash_option = 'winbash'

            info_handle.write('home\t%s\nbash\t%s\nshell\t%s' % (self.files['home'], self.files['bash'], bash_option))

            info_handle.close()
            self.master.quit()
        else:
            message_text = 'Both "home" and "bash" directories must be selected.'
            self.submit_warning = Label(
                self.master,
                text=message_text, fg='#FF6D00')
            self.submit_warning.grid(row=15, columnspan=3, padx=5, pady=5)


def run():
    """
    Function runs start_info
    :return:
    """
    root = Tk()
    start_app(root)

    root.mainloop()
    root.destroy()
