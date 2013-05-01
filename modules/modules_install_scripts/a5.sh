


###################
################### A5
###################
VERSION="20120518"
INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/a5/
mkdir -p $INSTALL_PATH
wget http://ngopt.googlecode.com/files/ngopt_a5pipeline_linux-x64_$VERSION.tar.gz
tar -xvf ngopt_a5pipeline_linux-x64_$VERSION.tar.gz
# sh ngopt_a5pipeline_linux-x64_$VERSION/test.a5.sh # only tested on login node


# TODO manual : modify the perl script to support long reads... -l + -r with empty file
	# PATCH: create dumy emptyfile	
	#my $idba_cmd = "idba -r $reads -o $WD/$outbase --mink ".IDBA_MIN_K." --maxk $maxrdlen";
	#my @dummy = `echo "" > empty.fa`;
	#my $idba_cmd = "idba -r empty.fa -l $reads -o $WD/$outbase --mink ".IDBA_MIN_K." --maxk $maxrdlen";
	

mv ngopt_a5pipeline_linux-x64_$VERSION $INSTALL_PATH

# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - A5 microbial genomes assembly pipeline \"
}
module-whatis \"MUGQIC - A5  \"
 
conflict        mugqic/bwa mugqic/fish mugqic/idba mugqic/samtools mugqic/sga mugqic/sspace mugqic/tagdust mugqic/bowtie mugqic/SSPACE             
set             root               \$::env(MUGQIC_INSTALL_HOME)/software/a5/ngopt_a5pipeline_linux-x64_$VERSION/bin
prepend-path    PATH               \$root
prepend-path    PATH               \$root/SSPACE
" > $VERSION

# version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"
" > .version

mkdir -p $MUGQIC_INSTALL_HOME/modulefiles/mugqic/a5
mv .version $VERSION $MUGQIC_INSTALL_HOME/modulefiles/mugqic/a5/





