#!/bin/bash

THIS_SCRIPT=$(basename "$0")
# Default GenPipes in a Container version to be used.
GIAC_VERSION=v4.0.0

usage() {
  echo "script usage: $THIS_SCRIPT -h [-v GiaC_version]"
  echo "Usage:"
  echo " -h                               Display this help message."
  echo " -v <GiaC_version>                GenPipes in a Container released version to be used."
  exit 1
  }

while getopts 'h:v:' OPTION; do
  case "$OPTION" in
  v)
    GIAC_VERSION="$OPTARG"
    ;;
  h)
    usage
    ;;
  ?)
    usage
    ;;
  esac
done

# Get Genpipes In A Container image. Cf. https://github.com/c3g/genpipes_in_a_container
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
CONTAINER_DIR=$SCRIPT_DIR/../../resources/container
mkdir "$CONTAINER_DIR"
wget -c http://github.com/c3g/genpipes_in_a_container/releases/download/${GIAC_VERSION}/wrapper_genpipes.tgz -O - | tar -xz -C ${CONTAINER_DIR} --strip-components=1
echo "Installing singularity image and configuration file in $CONTAINER_DIR"
if test -f "${CONTAINER_DIR}/etc/wrapper.conf"; then
  rel=$(realpath --relative-to=$PWD ${CONTAINER_DIR}/etc/wrapper.conf)
  read -r -p  "./$rel exist, do you want to overwrite it y[n]:?" response
   case "$response" in
     [yY][eE][sS]|[yY])
          echo "Updating config file."
          mv  "${CONTAINER_DIR}/etc/wrapper.conf.tpl" "${CONTAINER_DIR}/etc/wrapper.conf"
          ;;
     *)
          echo "Keeping old config file."
          ;;
    esac
else
   mv -i "${CONTAINER_DIR}/etc/wrapper.conf.tpl" "${CONTAINER_DIR}/etc/wrapper.conf"
fi
