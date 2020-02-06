#!/bin/bash
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
wget -c http://github.com/c3g/genpipes_in_a_container/releases/download/1.0.3/wrapper_genpipes.tgz -O - | tar -xz -C $SCRIPT_DIR   --strip-components=1

