#!/usr/bin/env python
"""
Script that verifies the installation scripts.

"""
import os
import sys
import csv
import subprocess
import urllib.request

class moduleVeri:
    def __init__(self, veri_file):
        self.veri_file = veri_file
        if not os.path.isfile('module_helpers/verify_modules.csv'):
            with open('module_helpers/verify_modules.csv', 'w') as file:
                writer = csv.writer(file)
                writer.writerow(['FILENAME', 'FOLLOW_TEMPLATE',
                                 'ARCHIVE_URL', 'URL_WORK?'])
        self.file_path = 'module_helpers/verify_modules.csv'
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
        with open('module_helpers/verify_modules.csv', 'a') as output_file:
            dict_writer = csv.DictWriter(output_file, keys)
            dict_writer.writerows(list_dict)

sw_name = sys.argv[1]
modVer = moduleVeri(sw_name)
modVer.create_bash()
url_ = modVer.get_var_value('ARCHIVE_URL')
arr_ = []
dict_ = {'FILENAME': sw_name,
         'FOLLOW_TEMPLATE': modVer.check_template(),
         'ARCHIVE_URL': url_,
         'URL WORK?': modVer.verify_url(url_)}

arr_.append(dict_)
modVer.append_to_csv(arr_)