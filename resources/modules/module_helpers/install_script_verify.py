#!/usr/bin/env python
"""
Script that verifies the installation scripts.

Can be run from command line. Stores the logfile
inside the folder where it is run from.

python install_script_verify.py <FILENAME>

"""
import os
import csv
import argparse
import subprocess
import urllib.request

class ModuleVeri:
    """
    A class made to verify the module mirrors in the 
    installation scripts.

    ...

    Attributes
    ----------
    veri_file : str
    Name of the file that is to be verified.
    
    vars : dict
    The variables found inside the installation script
    and their values.

    Methods
    -------
    create_bash()
    Creates a temporary bash file.
    
    check_template()
    Checks if the installation script adheres to
    the installation template
    
    get_var_value(varname)
    Gets the value of a certain variable inside the
    bash file.

    verify_url(url)
    Pings a certain URL to verify if the URL exists
    or not.

    append_to_csv(list_dict)
    Appends the list_dict to the logfile.

    """
    def __init__(self, veri_file):
        self.veri_file = veri_file
        if not os.path.isfile('verify_modules.csv'):
            with open('verify_modules.csv', 'w') as file:
                writer = csv.writer(file)
                writer.writerow(['FILENAME', 'FOLLOW_TEMPLATE',
                                 'ARCHIVE_URL', 'URL_WORK'])
        with open(self.veri_file, 'r') as f:
            self.data = f.readlines()
        self.vars = {}

    def create_bash(self):
        """
        Parses the installation scripts and generates a smaller 
        BASH file which is used to extract the variables from.

        """        
        temp_bash = []
        flag = False
        for line in self.data:
            temp_bash.append(line)
            if 'ARCHIVE_URL' in line.strip()[0:11]:
                temp_bash.append(line)
                flag = True
                break
        if flag:
            with open('temp.sh', 'w') as f:
                for line in temp_bash:
                    f.write(line)
        else:
            raise Exception("Cannot generate temporary bash file")

    def check_template(self):
        """
        Checks if the installation script follows the installation
        script template.

        Output
        ------
        data: boolean
        Returns True or False depending if the installation script
        follows the template.

        """
        flag = False
        for line in self.data:
            if '$MODULE_INSTALL_SCRIPT_DIR/install_module.sh' in line:
                flag = True
        return flag

    def get_var_value(self, varname):
        """
        Locally runs the temporary bash script generated to get
        variable values.
        
        Parameters
        ----------
        varname: str
        Variable name.

        Returns 
        -------
        data: str
        Returns the value of the variable.

        """
        CMD = 'echo $(source {}; echo ${})'.format('temp.sh', varname)
        p = subprocess.Popen(CMD, stdout=subprocess.PIPE, shell=True, executable='/bin/bash')
        return p.stdout.readlines()[0].strip().decode('utf-8')

    def verify_url(self, url):
        """
        Checks if the ARCHIVE_URL is a valid mirror or not.

        Parameters
        ----------
        url: str
        The URL which is to be verified.

        Returns 
        -------
        data: bool/str
        Returns the response.code if the response has failed,
        or a boolean value of True.

        """
        try:
            response = urllib.request.urlopen(url, timeout=10)
        except Exception as e:
            return e
        if response.code == 200:
            return True
        return response.code

    def append_to_csv(self, list_dict):
        """
        Appends the installation script verification results
        to the CSV.

        Parameters
        ----------
        list_dict: Consists of a dict containing different keys
        which are logged 

        """
        keys = list_dict[0].keys()
        with open('verify_modules.csv', 'a') as output_file:
            dict_writer = csv.DictWriter(output_file, keys)
            dict_writer.writerows(list_dict)

if __name__=='__main__':
    parser = argparse.ArgumentParser(description='Single script to scrape metadata of all\
        the software found inside the stack folder')
    parser.add_argument('path', metavar='P', type=str, nargs=1,
        help='Path to the installation script to be checked')

    args = parser.parse_args()
    sw_name = args.path[0]

    modver_ = ModuleVeri(sw_name)
    modver_.create_bash()
    url_ = modver_.get_var_value('ARCHIVE_URL')
    arr_ = []
    dict_ = {'FILENAME': sw_name,
             'FOLLOW_TEMPLATE': modver_.check_template(),
             'ARCHIVE_URL': url_,
             'URL_WORK': modver_.verify_url(url_)}

    arr_.append(dict_)
    modver_.append_to_csv(arr_)