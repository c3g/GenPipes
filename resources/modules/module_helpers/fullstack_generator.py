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

def generate_all_jsons(traverse_this=None, overwrite=False):
	if traverse_this is None:
		traverse_this = os.getcwd()
	all_soft = os.listdir(traverse_this)
	for soft in all_soft:
		print("Scraping {}".format(soft))
		mdh.run_mdh(soft, path=traverse_this, overwrite=False)

if __name__ == '__main__':
	try:
		sys.argv[1]
	except Exception as IndexError:
		print("Are you sure you're using it with proper arguments?\n \
			Try using -h or --help for help".replace('			', ''))
		sys.exit()
	if sys.argv[1] == '-h' or '--help':
		print("Driver Script to generate .metadata.json for full stack\n\
		Usage:\n\
		python fullstack_generator.py <PATH> <OVERWRITE_FLAG=True/False\
		(False by Default)".replace('		', ''))
		sys.exit()
	path = sys.argv[1]
	overwrite_flag = False
	if len(sys.argv) == 3:
		overwrite_flag = sys.argv[2]
	print(f'Traversing {path}')
	generate_all_jsons(traverse_this=path, overwrite=overwrite_flag)
