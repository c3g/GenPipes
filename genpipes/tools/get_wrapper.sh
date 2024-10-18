#!/bin/bash
# Get Genpipes In A Container image
GIAC_VERSION=v2.1.0
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
CONTAINER_DIR=$SCRIPT_DIR/../../resources/container
mkdir "$CONTAINER_DIR"
wget -c http://github.com/c3g/genpipes_in_a_container/releases/download/${GIAC_VERSION}/wrapper_genpipes.tgz -O - | tar -xz -C ${CONTAINER_DIR} --strip-components=1
echo Installing singularity image and configuration file
if test -f "${CONTAINER_DIR}/etc/wrapper.conf"; then
  rel=$(realpath --relative-to=$PWD ${CONTAINER_DIR}/etc/wrapper.conf)
  echo "./$rel exist,"
  read -r -p  "  do you want to overwrite it y[n]:?" response

   case "$response" in
     [yY][eE][sS]|[yY])
          echo updating confing
          mv  "${CONTAINER_DIR}/etc/wrapper.conf.tpl" "${CONTAINER_DIR}/etc/wrapper.conf"
          ;;
     *)
          echo keeping old config
          ;;
    esac
else
   mv -i "${CONTAINER_DIR}/etc/wrapper.conf.tpl" "${CONTAINER_DIR}/etc/wrapper.conf"
fi


