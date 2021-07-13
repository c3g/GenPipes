#################################################################
##################### Setup of dependencies	#####################
#################################################################

## 'MUGQIC_INSTALL_HOME_DEV' for development, 'MUGQIC_INSTALL_HOME' for production (don't write '$' before!)
if [[ ${1:-} == MUGQIC_INSTALL_HOME ]]
then
  INSTALL_HOME=MUGQIC_INSTALL_HOME
elif [[ ${1:-} == GENPIPES_INSTALL_HOME ]]
then
  INSTALL_HOME=GENPIPES_INSTALL_HOME
else
  INSTALL_HOME=MUGQIC_INSTALL_HOME_DEV
  NOPATCH=1
  NOWRAP=1
fi

## Indirection call to use $INSTALL_HOME value as variable name
INSTALL_DIR=${!INSTALL_HOME}/software/$SOFTWARE
ARCHIVE_DIR=${!INSTALL_HOME}/archive

WD="$INSTALL_HOME/software/magic_infinium/"
mkdir -p $WD
cd $WD
mkdir -p software
cd software

## IAAP-CLI linux (login at Illumina.com required) [module]
mkdir -p iaap
cd iaap
wget "https://www.dropbox.com/s/q53buu94ds2z4cc/iaap-cli-linux-x64-1.1.0-sha.80d7e5b3d9c1fdfc2e99b472a90652fd3848bbc7.tar.gz?dl=1" -O "iaap-cli-linux-x64-1.1.0-tar.gz"
tar -xvf iaap-cli-linux-x64-1.1.0-tar.gz

## BeadArrayFiles [module]
cd /lb/project/C3G/projects/auld_infinium_dev/software
mkdir -p BeadArrayFiles
cd BeadArrayFiles
wget https://github.com/Illumina/BeadArrayFiles/archive/1.3.3.tar.gz
tar -xvf 1.3.3.tar.gz
cd BeadArrayFiles-1.3.3/
module load mugqic/python/3.7.3
python setup.py install --user

# custom scripts
wget https://www.dropbox.com/s/81m24491n5iieko/dna_summary_magic.py?dl=1  -O "dna_summary_magic.py" . # TODO: mugqic tools?
wget https://www.dropbox.com/s/7ivq8hstdq6gwnc/gtc_control_report.py?dl=1 -O "gtc_control_report.py" . # TODO: mugqic tools?

## pysam [python module]
module load mugqic/python/3.7.3
#cd /lb/project/C3G/projects/auld_infinium_dev/software
pip install --install-option="--prefix=$WD/pysam"  --install-option="--ignore-installed" pysam
#mkdir -p pysam
#cd pysam
#wget https://files.pythonhosted.org/packages/99/5a/fc440eb5fffb5346e61a38b49991aa552e4b8b31e8493a101d2833ed1e19/pysam-0.16.0.1.tar.gz
#tar -xvf  pysam-0.16.0.1.tar.gz
#cd  pysam-0.16.0.1
#python setup.py install --user

# pyvcf [python module]
module load mugqic/python/3.7.3
cd /lb/project/C3G/projects/auld_infinium_dev/software
mkdir -p pyvcf
cd pyvcf
wget "https://files.pythonhosted.org/packages/20/b6/36bfb1760f6983788d916096193fc14c83cce512c7787c93380e09458c09/PyVCF-0.6.8.tar.gz"
tar -xvf  "PyVCF-0.6.8.tar.gz"
cd "PyVCF-0.6.8"
python setup.py install --user

## Illumina GTC2VCF [module]
cd /lb/project/C3G/projects/auld_infinium_dev/software
mkdir -p GTCtoVCF
cd GTCtoVCF
wget https://github.com/Illumina/GTCtoVCF/archive/1.2.1.tar.gz
tar -xvf 1.2.1.tar.gz 
cd GTCtoVCF-1.2.1
