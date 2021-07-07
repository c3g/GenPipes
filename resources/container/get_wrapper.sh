#!/bin/bash
# Get Genpipes In A Container image
GIAC_VERSION=v2.0.1
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
wget -c http://github.com/c3g/genpipes_in_a_container/releases/download/${GIAC_VERSION}/wrapper_genpipes.tgz -O - | tar -xz -C ${SCRIPT_DIR}   --strip-components=1
echo Installing singulaity image and configuration file
if test -f "${SCRIPT_DIR}/etc/wrapper.conf"; then
  rel=$(realpath --relative-to=$PWD ${SCRIPT_DIR}/etc/wrapper.conf)
  echo "./$rel exist,"
  read -r -p  "  do you want to overwrite it y[n]:?" response

   case "$response" in
     [yY][eE][sS]|[yY])
          echo updating confing
          mv  ${SCRIPT_DIR}/etc/wrapper.conf.tpl ${SCRIPT_DIR}/etc/wrapper.conf
          ;;
     *)
          echo keeping old config
          ;;
    esac
else
  echo mv -i ${SCRIPT_DIR}/etc/wrapper.conf.tpl ${SCRIPT_DIR}/etc/wrapper.conf
fi


