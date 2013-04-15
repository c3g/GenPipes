

###################
################### Ray
###################
#module load compat-openmpi-x86_64 # 4abacus
#module load compat-openmpi-psm-x86_64
VERSION="2.1.0"
INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/ray/Ray-v$VERSION
NAME=Ray-v$VERSION
wget http://sourceforge.net/projects/denovoassembler/files/$NAME.tar.bz2/download
tar -xvf $NAME.tar.bz2
cd $NAME
# ...... TODO: meh
make PREFIX=ray-build
make install
ls ray-build
mpiexec -n 1 ray-build/Ray -o test -p test_1.fastq test_2.fastq -k 31
#mkdir -p $INSTALL_PATH
Then, type
== Options ==
You can provide compilation options to the Makefile.
MPICXX                  The path to the C++ compiler wrapper (usually called mpicxx)
PREFIX                  Where to install stuff
MAXKMERLENGTH           maximum k-mer length, default is MAXKMERLENGTH=32
FORCE_PACKING           save memory by not aligning addresses, default is FORCE_PACKING=n
ASSERT                  run assertions too, default is ASSERT=n
For other options, read the Makefile header.




