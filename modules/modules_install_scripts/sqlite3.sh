
###################
################### sqlite3
###################
VERSION="3071502"
INSTALL_PATH="$MUGQIC_INSTALL_HOME/software/sqlite-shell-linux-x86-"$VERSION
mkdir -p $INSTALL_PATH
FILE="sqlite-shell-linux-x86-"$VERSION".zip" #http://www.sqlite.org/sqlite-shell-linux-x86-3071502.zip
wget http://www.sqlite.org/$FILE
unzip $FILE
mv sqlite3  $INSTALL_PATH

# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - sqlite3 \"
}
module-whatis \"MUGQIC - sqlite3 shell \"
                      
set             root               \$::env(MUGQIC_INSTALL_HOME)/software/sqlite-shell-linux-x86-$VERSION
prepend-path    PATH               \$root
" > $VERSION

# version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"
" > .version

mkdir -p $MUGQIC_INSTALL_HOME/modulefiles/mugqic/sqlite3
mv .version $VERSION $MUGQIC_INSTALL_HOME/modulefiles/mugqic/sqlite3/


