#!/usr/bin/env python
'''
This is a single script to scrape metadata of all the software
within a specific folder. 

Pass the path as a command line argument.

Eg - python fullstack_generator.py <PATH>

The JSONs are stored inside the module folders inside the given
path.
'''
import os
import argparse
import metadata_helper as mdh

def generate_all_jsons(traverse_this=None, overwrite=False):
	if traverse_this is None:
		traverse_this = os.getcwd()
	all_soft = os.listdir(traverse_this)
	for soft in all_soft:
		print("Scraping {}".format(soft))
		mdh.run_mdh(soft, path=traverse_this, overwrite=overwrite)

if __name__ == '__main__':	
	parser = argparse.ArgumentParser(description='Single script to scrape metadata of all\
		the software found inside the stack folder')
	parser.add_argument('path', metavar='P', type=str, nargs=1,
		help='Path to the software stack')
	parser.add_argument('-o', '--overwrite',
		help='Use the flag only if you want to overwrite previous JSONs',
		action='store_true')
	parser.add_argument('-v', '--verbose',
		help='Use verbose flag to print the directories as they are crawled',
		action='store_true')
	
	args = parser.parse_args()
	path = args.path[0]
	if args.verbose:
		print(f'Traversing {path}')
	
	generate_all_jsons(traverse_this=path, overwrite=args.overwrite)



