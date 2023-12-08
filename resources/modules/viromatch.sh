#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=viromatch
VERSION=3.2.3.1
ARCHIVE=$SOFTWARE-$VERSION.tar.gz
ARCHIVE_URL=https://github.com/nicolargo/$SOFTWARE/archive/refs/tags/v$VERSION.tar.gz
SOFTWARE_DIR=${SOFTWARE^}-$VERSION
PYTHON_VERSION=3.9.1
PYTHON_SHORT_VERSION=${PYTHON_VERSION:0:3}
NOWRAP=1
NOPATCH=1

build() {
  cd $INSTALL_DOWNLOAD
  module load mugqic/python/$PYTHON_VERSION
  pip install --prefix=$INSTALL_DIR/$SOFTWARE_DIR --ignore-installed ${SOFTWARE}[action,browser,cloud,cpuinfo,docker,export,folders,gpu,graph,ip,raid,snmp,web,wifi,ratelimiter,datrie,appdirs,ConfigArgParse]==${VERSION}
  pip install --prefix=$INSTALL_DIR/$SOFTWARE_DIR --ignore-installed snakemake==5.25.0
  ln -s $(which python) $INSTALL_DIR/$SOFTWARE_DIR/bin/python
  ln -s $(which python3) $INSTALL_DIR/$SOFTWARE_DIR/bin/python3
  # Make sure all shebang are ok as it's a docker software
  sed -i -e 's@#! /usr/bin/python3.7@#!/usr/bin/env python3.7@g' $INSTALL_DIR/$SOFTWARE_DIR/bin/*.py
  # Remove hardcoded path to smk file
  sed -i -e 's@/usr/lib/python3.7@/cvmfs/soft.mugqic/CentOS6/software/viromatch/viromatch-master44edca0@g' $INSTALL_DIR/$SOFTWARE_DIR/viromatch/lib/pipeline/viromatch.py
  # Make viromatch able to import itself
  sed -i -e 's@#!/usr/bin/env python3.7@#!/usr/bin/env python3.7\n\nimport sys\nsys.path.append("/cvmfs/soft.mugqic/root/software/viromatch/viromatch-master44edca0")\n@g' $INSTALL_DIR/$SOFTWARE_DIR/bin/taxid2lineage.py
  sed -i -e 's@#!/usr/bin/env python3.7@#!/usr/bin/env python3.7\n\nimport sys\nsys.path.append("/cvmfs/soft.mugqic/root/software/viromatch/viromatch-master44edca0")\n@g' $INSTALL_DIR/$SOFTWARE_DIR/bin/viromatch.py
  sed -i -e 's@#!/usr/bin/env python3.7@#!/usr/bin/env python3.7\n\nimport sys\nsys.path.append("/cvmfs/soft.mugqic/root/software/viromatch/viromatch-master44edca0")\n@g' $INSTALL_DIR/$SOFTWARE_DIR/bin/tsv_taxonomy_count_prep.py
  sed -i -e 's@#!/usr/bin/env python3.7@#!/usr/bin/env python3.7\n\nimport sys\nsys.path.append("/cvmfs/soft.mugqic/root/software/viromatch/viromatch-master44edca0")\n@g' $INSTALL_DIR/$SOFTWARE_DIR/bin/sam_taxonomy_count_prep.py
  sed -i -e 's@#!/usr/bin/env python3.7@#!/usr/bin/env python3.7\n\nimport sys\nsys.path.append("/cvmfs/soft.mugqic/root/software/viromatch/viromatch-master44edca0")\n@g' $INSTALL_DIR/$SOFTWARE_DIR/bin/best_hit_filter_tsv.py
  sed -i -e 's@#!/usr/bin/env python3.7@#!/usr/bin/env python3.7\n\nimport sys\nsys.path.append("/cvmfs/soft.mugqic/root/software/viromatch/viromatch-master44edca0")\n@g' $INSTALL_DIR/$SOFTWARE_DIR/bin/best_hit_filter_sam.py
}

module_file() {
echo "\
#%Module1.0
proc ModulesHelp { } {
  puts stderr \"\tMUGQIC - $SOFTWARE \"
}
module-whatis \"$SOFTWARE\"

set             root                $INSTALL_DIR/$SOFTWARE_DIR
prepend-path    PATH                \$root/bin
prepend-path    PATH                /cvmfs/soft.mugqic/CentOS6/software/bwa/bwa-0.7.17
prepend-path    PATH                /cvmfs/soft.mugqic/CentOS6/software/bwa/bwa-0.7.17/bwakit
prepend-path    PATH                /cvmfs/soft.mugqic/CentOS6/software/samtools/samtools-1.9/bin
prepend-path    LIBRARY_PATH        /cvmfs/soft.mugqic/CentOS6/software/samtools/samtools-1.9/lib
prepend-path    LD_LIBRARY_PATH     /cvmfs/soft.mugqic/CentOS6/software/samtools/samtools-1.9/lib
prepend-path    PATH                /cvmfs/soft.mugqic/CentOS6/software/seqtk/seqtk-1.2
prepend-path    PATH                /cvmfs/soft.mugqic/CentOS6/software/diamond/diamond-2.1.4/bin
prepend-path    PATH                /cvmfs/soft.mugqic/CentOS6/software/vsearch/vsearch-1.11.1
prepend-path    PATH                /cvmfs/soft.mugqic/CentOS6/software/fqtrim/fqtrim-0.9.7/bin
prepend-path    PYTHONPATH          $PYTHONPATH
prepend-path    PYTHONPATH          \$root/lib/python${PYTHON_SHORT_VERSION}
prepend-path    PYTHONPATH          \$root/lib/python${PYTHON_SHORT_VERSION}/site-packages
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@