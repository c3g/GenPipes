#!/bin/bash
# Get Genpipes In A Container image
GIAC_VERSION=v2.0.0.rc
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
wget -c http://github.com/c3g/genpipes_in_a_container/releases/download/${GIAC_VERSION}/wrapper_genpipes.tgz -O - | tar -xz -C ${SCRIPT_DIR}   --strip-components=1
echo Installing singulaity image and configuration file
echo If ${SCRIPT_DIR}/etc/wrapper.conf is already configured, you should probably not overwirte it
mv -i ${SCRIPT_DIR}/etc/wrapper.conf.tpl ${SCRIPT_DIR}/etc/wrapper.conf

