"""
A set of web scraping scripts.
Focused on scraping software metadata from various different websites.
logfile variables is the destination where the logfile will be stored.
"""

import time
import os
import requests
import lxml.html
from bs4 import BeautifulSoup

logfile = os.getcwd()

class searchLogger:
    """
    The function handles logging done when scrapers are run.
    """
    def __init__(self):
        self.log_format = '{time} - {type} - {message}\n'
        self.logfile_path = os.path.join(logfile, 'logfile.txt')
        if not os.path.isfile(self.logfile_path):
            with open(self.logfile_path, 'w') as file:
                curr_time = time.ctime(time.time())
                file.write(f'LOG FILE GENERATED ON {curr_time}\n')
        
    def append_log(self, dict_data):
        """
        Appends a single line to the Logfile.

        Parameters
        ----------
        dict_data : {'time': time,
                    'type': type of the message,
                    'message': content of the log}
            Contains the log message.
        """
        curr_time = time.ctime(time.time())
        dict_data.update({'time': curr_time})
        with open(self.logfile_path, 'a') as file:
            file.write(self.log_format.format(**dict_data))
    
    def reset_log(self):
        """
        Removes the old logfile, and creates a new one.
        """
        with open(self.logfile_path, 'w') as file:
                curr_time = time.ctime(time.time())
                file.write(f'LOG FILE GENERATED ON {curr_time}\n')

class searchBioconda:
    """
    Scraper made to scrape Bioconda data from - 
    https://bioconda.github.io/conda-package_index.html
    """
    def __init__(self):
        self.keys_to_ignore = ['Links',
                               'Recipe',
                               'DEPENDENCIES',
                               'VERSIONS_AVAILABLE',
                               'Documentation',
                               'Developer docs',
                               'Homepage']
        root_source = 'https://bioconda.github.io/conda-package_index.html'
        page_source = requests.get(root_source)
        page_soup = BeautifulSoup(page_source.text, features='lxml')
        self.xpath_dict = {'INFO': '/html/body/div[1]/div[2]/div/div/div/dl[1]/dd/p',
             'VERSIONS_AVAILABLE': '/html/body/div[1]/div[2]/div/div/div/dl[2]/dd/dl/dd[1]',
             'DEPENDENCIES': '/html/body/div[1]/div[2]/div/div/div/dl[2]/dd/dl/dd[2]'}
        self.all_keys = ['CHANNEL_LINK', 'NAME', 'INFO', 'VERSIONS_AVAILABLE', 
                               'Homepage', 'Documentation', 'Developer docs', 'License']
        self.check_if_exist = list(set(self.all_keys) - set(self.keys_to_ignore))
        sws_ = page_soup.find('table', class_='indextable modindextable').findAll('a')    
        sws_dict = {}
        for elem in sws_:
            sws_dict[elem.text.strip()] = 'https://bioconda.github.io/' + elem.get('href')
        self.db = sws_dict
        self.logger = searchLogger()
        self.logger.append_log({
            'type': 'INITIALIZATION',
            'message': 'Database scraped from Bioconda Index'
            })

    def return_db(self):
        """
        Returns the database set scraped from - 
        https://bioconda.github.io/conda-package_index.html
        
        Returns 
        -------
        data: dict
        A dict with keys as software name, and values as
        links.

        """
        return self.db
    
    def fetch_software(self,  name):
        """
        Fetches the software metadata from their channel page.
        
        Parameters
        ----------
        name: str
        Name of the software that is to be extracted.

        Returns 
        -------
        data: dict
        A dict with keys as software name, and values as
        links.

        """
        name = name.lower()
        template_dict = {}
        if name in list(self.db.keys()):
            sw_link = self.db[name]
            template_dict['CHANNEL_LINK'] = sw_link
            template_dict['NAME'] = name
            source = requests.get(sw_link)
            bs4_source = BeautifulSoup(source.text, features='lxml')
            lxml_source = lxml.html.fromstring(source.content)
            temp_box = bs4_source.find('dl', class_='field-list simple')
            for key in list(self.xpath_dict.keys()):
                try:
                    template_dict[key] = lxml_source.xpath(self.xpath_dict[key])[0].text_content().strip().replace('\n', ',')
                except Exception as e:
                    self.logger.append_log({
                        'type': 'Exception @ {}'.format(name),
                        'message': e
                    })
                    print("Exception encountered! {}, exception - {}".format(key, e))
            keys = lxml_source.xpath('/html/body/div[1]/div[2]/div/div/div/dl[1]/dd/dl/dt')
            values = lxml_source.xpath('/html/body/div[1]/div[2]/div/div/div/dl[1]/dd/dl/dd')
            for index, key in enumerate(keys):
                template_dict[key.text_content()] = values[index].text_content().strip().replace('\n', ',')
            for key in self.keys_to_ignore:
                template_dict.pop(key, None)
            self.logger.append_log({
                'type': 'Message @ {}'.format(name),
                'message': 'Extracted the following {}'.format(list(template_dict.keys()))
            })
            keys_not_found = list(set(self.check_if_exist) - set(template_dict.keys()))
            if keys_not_found:
                self.logger.append_log({
                    'type': 'Not Found @ {}'.format(name),
                    'message': 'Did not find the following {}'.format(list(keys_not_found))
                })
            return template_dict
        return False
    
    def isBioconda(self, name):
        """
        Checks if a certain software is available in Bioconda Index.
        
        Parameters
        ----------
        name: str
        Name of the software that is to be checked.

        Returns 
        -------
        data: boolean
        Returns a True or False depending on the result.

        """
        if name.lower() in list(self.db.keys()):
            return True
        return False


class searchPyPi:
    """
    Scraper made to scrape metadata from PyPi. 
    """
    def __init__(self):
        self.keys_to_ignore = ['MAINTAINER', 'AUTHOR']
        self.search_source = 'https://pypi.org/search/?q={}'
        self.field_paths = {'INFO': '/html/body/main/div[2]/div/div',
                           'MAINTAINER': '/html/body/main/div[3]/div/div/div[1]/div[5]/span/a/span[2]',
                           'AUTHOR': '/html/body/main/div[3]/div/div/div[1]/div[4]/p/a',
                           'HOMEPAGE': '/html/body/main/div[3]/div/div/div[1]/div[2]/ul/li/a', 
                           'LICENSE': '/html/body/main/div[3]/div/div/div[1]/div[4]/p[1]'}
        self.search_xpath = '/html/body/main/div/div/div[2]/form/div[3]/ul/li/a'
        self.keys_to_ignore = ['AUTHOR', 'MAINTAINER', 'HOMEPAGE']
    
    def searchPackage(self, name):
        """
        Checks if a certain software is available on PyPi.
        
        Parameters
        ----------
        name: str
        Name of the software that is to be searched.

        Returns 
        -------
        data: dict
        Returns a dict containing the software as the key
        and URL as the value..

        """
        queried_page = self.search_source.format(name)
        page_source = requests.get(queried_page)
        lxml_html = lxml.html.fromstring(page_source.content)
        search_results = lxml_html.xpath(self.search_xpath)
        result_dict = {}
        for result in search_results:
            sw_name = result.text_content().strip().split(' ')[0]
            sw_name = sw_name.replace('\n', '').lower()
            if sw_name == name: 
                return 'https://pypi.org' + result.get('href')
        return None
    
    def getMetadata(self, name, link):
        """
        Extracts the software metadata from the PyPi page.
        
        Parameters
        ----------
        name: str
        Name of the software.

        link: str
        URL of the PyPi page.

        Returns 
        -------
        data: dict
        Returns a dict containing the software metadata.

        """
        queried_page = requests.get(link)
        lxml_html = lxml.html.fromstring(queried_page.content)
        sw_dict = {'NAME': name,
                  'CHANNEL_LINK': link}
        for field in list(self.field_paths.keys()):
            try:
                if field in self.keys_to_ignore:
                    continue
                sw_dict[field] = lxml_html.xpath(self.field_paths[field])[0].text_content().strip()                
                if field is 'HOMEPAGE':
                    sw_dict[field] = lxml_html.xpath(self.field_paths[field])[0].get('href')
            except Exception as e:
                print("Exception encountered! {}, exception - {}".format(field, e))
        return sw_dict