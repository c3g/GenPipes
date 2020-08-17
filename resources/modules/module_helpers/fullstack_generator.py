'''
This is a single script to scrape metadata of all the software
within a specific folder. 

Pass the path as a command line argument.

Eg - 

python fullstack_generator.py <PATH>

The JSONs are stored inside the module folders inside the given
path.
'''
import os
import sys
import metadata_helper as mdh

def generate_all_jsons(traverse_this=None):
	if traverse_this is None:
		traverse_this = os.getcwd()
	all_soft = os.listdir(traverse_this)
	for soft in all_soft:
		print("Scraping {}".format(soft))
		mdh.run_mdh(soft, path=traverse_this)

if __name__ == '__main__':
	path = sys.argv[1]
	print(f'Traversing {path}')
	generate_all_jsons(traverse_this=path)
