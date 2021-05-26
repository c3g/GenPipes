#!/bin/bash
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
wget -c http://github.com/c3g/genpipes_in_a_container/releases/download/1.1.0/wrapper_genpipes.tgz -O - | tar -xz -C ${SCRIPT_DIR}   --strip-components=1
echo Installing singulaity image and configuration file
mv -i ${SCRIPT_DIR}/etc/wrapper.conf.tpl ${SCRIPT_DIR}/etc/wrapper.conf

