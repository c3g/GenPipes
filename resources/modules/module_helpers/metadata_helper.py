#!/usr/bin/env python

"""
Driver Script to get a software's metadata.
The script is called by the installation scripts available.
Change c_path to destination you wish to store the JSONs.

Can be run directly from command line with a command line argument

python metadata_helper.py <SOFTWARE_NAME> <PATH_TO_SOFTSTACK>
"""

import os
import json
import argparse
import pythonSearcher as pySrch


def run_mdh(sw_name, path=False, overwrite=False):
	"""
	The pythonSearcher.py driver file.

	Parameters
	----------
	sw_name: str
	Name of the software to be extracted

	path: str
	Path to the software stack root folder

	overwrite: boolean
	True/False flag. If overwrite is True, then jsons
	are overwritten.

	Returns 
	-------
	data: str
	Returns the value of the variable.

	"""
	if path:
		c_path = path
	else:
		c_path = os.getcwd()

	logger = pySrch.SearchLogger()

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

	bioconda_search = pySrch.SearchBioconda()
	pypi_search= pySrch.SearchPyPi()

	bconda_data = bioconda_search.fetch_software(sw_name)
	pypi_url = pypi_search.search_package(sw_name)

	pypi = False
	not_found = False

	if bconda_data:
		data = bconda_data
	elif pypi_url:
		data = pypi_search.get_metadata(sw_name, pypi_url)
		pypi = True
	else:
		data = template_dict
		not_found = True 

	def_name = '.metadata.json'

	f_name = '.metadata.json'

	if pypi:
		f_name = 'VERIFY_' + f_name
	if not_found:
		f_name = 'NOTFOUND_' + f_name

	json_path = os.path.join((c_path), sw_name)
	check_json_path = os.path.join(json_path, def_name)
	json_filepath = os.path.join(json_path, f_name)

	if not os.path.exists(json_path):
		os.makedirs(json_path)

	if overwrite:
		with open(json_filepath, 'w') as f:
			json.dump(data, f, indent=6)
		print(f'JSON inside {json_filepath}')

	elif os.path.isfile(check_json_path):
		print(f'JSON skipped since exists {check_json_path}')
		pass

	else:
		with open(json_filepath, 'w') as f:
			json.dump(data, f, indent=6)
		print(f'JSON inside {json_filepath}')

	
if __name__ == '__main__':
	

	parser = argparse.ArgumentParser(description='Script to scrape metadata of a\
		package')
	parser.add_argument('name', metavar='N', type=str, nargs=1,
		help='Name of the software')
	parser.add_argument('path', metavar='P', type=str, nargs=1,
		help='Path to the software stack')
	parser.add_argument('-o', '--overwrite',
		help='Use the flag only if you want to overwrite previous JSONs',
		action='store_true')
	
	args = parser.parse_args()

	sw_name = args.name[0]
	path = args.path[0]

	run_mdh(sw_name, path=path, overwrite=args.overwrite)