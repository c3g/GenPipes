#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=biodiscml_report
VERSION=1.0
SOFTWARE_DIR=${SOFTWARE}-${VERSION}
# 3.10.4 (mugqic/python's default) has a broken pip on at least Narval
# (PermissionError reading pip's own vendored files from CVMFS -- a
# pre-existing issue with that specific module on that cluster, confirmed
# reproducible independently of this script). 3.11.1 was confirmed working
# there instead.
PYTHON_VERSION=3.11.1
PYTHON_SHORT_VERSION=${PYTHON_VERSION%.*}

# This module isn't a single piece of software: it pins the Python
# dependencies of the biodiscml pipeline's `report` step (see
# genpipes/bfx/biodiscml.py's `report()` and
# genpipes/pipelines/biodiscml/biodiscml_report.ipynb) on top of
# mugqic/python, from PyPI -- the same way resources/modules/python_lib.sh
# installs mugqic/python's own libraries. Versions are the ones validated end
# to end against real BioDiscML outputs (see biodiscml_report/requirements.txt
# for the unpinned minimums).
#
# There's no single upstream tarball for a pinned-deps module like this, but
# install_module.sh always downloads one before calling build() --
# nbconvert's own source archive (one of the pinned packages) is used as that
# placeholder, like cnvkit.sh points at cnvkit's tarball while installing
# from PyPI.
ARCHIVE=nbconvert-7.17.1.tar.gz
ARCHIVE_URL=https://github.com/jupyter/nbconvert/archive/refs/tags/v7.17.1.tar.gz
PACKAGES="numpy==2.4.6 pandas==3.0.6 scipy==1.17.1 scikit-learn==1.9.1 plotly==7.1.0 jupyter==1.1.1 nbconvert==7.17.1"

build() {
  cd $INSTALL_DOWNLOAD

  module load mugqic/python/$PYTHON_VERSION

  # On Digital Research Alliance clusters, the default environment points pip
  # at the Alliance's own wheelhouse (PIP_CONFIG_FILE) and puts a
  # _manylinux.py on PYTHONPATH that disables all manylinux wheels. Both are
  # meant for the Alliance's own Python builds: their wheels need a newer
  # glibc than mugqic/python runs with, and fail to import from it
  # (e.g. "GLIBC_2.34 not found"). mugqic/python needs PyPI's manylinux
  # wheels instead, so bypass both for this install. No effect elsewhere.
  # Installed in a single pip call so dependencies are resolved together.
  env -u PIP_CONFIG_FILE -u PYTHONPATH \
    pip install \
      --index-url https://pypi.org/simple \
      --only-binary=:all: \
      --prefix=$INSTALL_DIR/$SOFTWARE_DIR \
      --ignore-installed \
      $PACKAGES

  ln -s $(which python) $INSTALL_DIR/$SOFTWARE_DIR/bin/python
  ln -s $(which python3) $INSTALL_DIR/$SOFTWARE_DIR/bin/python3
}

module_file() {
echo "\
#%Module1.0
proc ModulesHelp { } {
  puts stderr \"\tMUGQIC - $SOFTWARE : Python environment (jupyter/plotly/scipy/scikit-learn) for the biodiscml pipeline's report step\"
}
module-whatis \"$SOFTWARE\"

set             root                $INSTALL_DIR/$SOFTWARE_DIR
prepend-path    PATH                \$root/bin
prepend-path    PYTHONPATH          \$root/lib/python${PYTHON_SHORT_VERSION}
prepend-path    PYTHONPATH          \$root/lib/python${PYTHON_SHORT_VERSION}/site-packages
# pip --prefix puts Jupyter's data files (nbconvert templates, kernel specs)
# and config under \$root, which Jupyter doesn't search by default
prepend-path    JUPYTER_PATH        \$root/share/jupyter
prepend-path    JUPYTER_CONFIG_PATH \$root/etc/jupyter
prereq          mugqic/python/$PYTHON_VERSION
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
