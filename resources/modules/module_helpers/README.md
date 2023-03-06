# C3G

## Description

The following set of scripts were written under the project, which's agenda was to generate an automated method to maintain the current software stack. 

**Link to PR**: https://bitbucket.org/mugqic/genpipes/pull-requests/181/soft-jsondb-gsoc2020/

### Independent Scripts

#### pythonSearcher.py

**INFO**: This script looks up the software's name on the internet, scrapes the important data, and then the scraped data is converted into a .metadata.json which is stored inside the software's folder inside the stack.

**USAGE**: 
```
pythonSearcher.py is a Python Class. To use it, the user is required to import it like a Python Package. 
```

#### json_combine.py

**INFO**: This script combines all the JSONs found in the software stack folders into a single combined_json.json

**USAGE**: 
```
python json_combine.py <PATH_TO_SOFT_STACK> <PATH_TO_STORE_JSON> ( Can be called through CLI )
```

#### install_script_verify.py

**INFO**: It verifies the installation bash scripts by going through their mirrors and verifying them. It then generates a logfile based on the data collected.

**USAGE**: 
```
python install_script_verify.py <PATH_TO_SCRIPT_NAME> ( Can be called through CLI )
```

### Driver Files

#### fullstack_generator.py

**INFO**: This script looks through the defined directory, and then generates metadata jsons for all the software folders found.

**USAGE**: 
```
python fullstack_generator.py <PATH> ( Can be called through CLI )
```

#### metadata_helper.py

**INFO**: Driver script to run pythonSearcher.py for a specific software. This has been exclusively written to be used in the bash scripts.

**USAGE**: 
```
python metadata_helper.py <SOFTWARE_NAME> <PATH_TO_SOFTSTACK> ( Can be called through CLI )
```

### FLOW

The following code scripts are made to run in a certain way. Below is a detailed pipeline of how the scripts should be executed.

1. **fullstack_generator.py** can be run on a software stack folder, to scrape metadata of all the packages inside the mentioned folder.
1. **json_combine.py** is to be run to generate a combined json file should it be needed.
1. **metadata_helper.py** can be run directly from the command line. It is made to be mentioned in the installation scripts such that when a certain installation script is run, it also scrapes the data and generates a .metadata.json for the software. It can be run directly as well, should it be needed.

Please note that *install_script_verify.py* is not a requirement, but can be run to verify whether an installation script's ARCHIVE_URL is functional or not.
