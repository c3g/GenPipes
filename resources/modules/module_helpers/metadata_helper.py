#!/usr/bin/env python

"""
Driver Script to get a software's metadata.
The script is called by the installation scripts available.
Change c_path to destination you wish to store the JSONs.

Can be run directly from command line with a command line argument

python metadata_helper.py <SOFTWARE_NAME> <PATH_TO_SOFTSTACK>
"""

import pythonSearcher as pySrch
import json
import sys
import os


def run_mdh(sw_name, path=False):

	if path:
		c_path = path
	else:
		c_path = os.getcwd()

	logger = pySrch.searchLogger()

	logger.append_log({
		'type': 'READY',
		'message': f'Scraping {sw_name}'
		})

	template_dict = {
	      "CHANNEL_LINK": None,
	      "NAME": sw_name,
	      "INFO": None,
	      "License": None
	}

	bioconda_search = pySrch.searchBioconda()
	pypi_search= pySrch.searchPyPi()

	bconda_data = bioconda_search.fetch_software(sw_name)
	pypi_url = pypi_search.searchPackage(sw_name)

	pypi = False
	not_found = False

	if bconda_data:
		data = bconda_data
	elif pypi_url:
		data = pypi_search.getMetadata(sw_name, pypi_url)
		pypi = True
	else:
		data = template_dict
		not_found = True 

	f_name = '.metadata.json'

	if pypi:
		f_name = 'VERIFY_' + f_name
	if not_found:
		f_name = 'NOTFOUND_' + f_name

	json_path = os.path.join((c_path), sw_name)
	json_filepath = os.path.join(json_path, f_name)

	if not os.path.exists(json_path):
		os.makedirs(json_path)

	with open(json_filepath, 'w') as f:
		json.dump(data, f, indent=6)

	print(f'JSON inside {json_filepath}')

if __name__ == '__main__':
	sw_name = sys.argv[1]
	if len(sys.argv) == 3:
		path_ = sys.argv[2]
	run_mdh(sw_name, path_)
