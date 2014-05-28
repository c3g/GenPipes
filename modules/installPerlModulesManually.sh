#!/bin/bash

mkdir -p packages
cd packages

## Define install path here:
INSTALLPATH=$MUGQIC_INSTALL_HOME/software/perl/perl-5.18.2

## Download modules here. Try to use <perl -MCPAN -e 'install This::Module'>. If it doesnt work, install it by hand.
wget http://search.cpan.org/CPAN/authors/id/G/GR/GRANTM/XML-Simple-2.20.tar.gz
wget http://search.cpan.org/CPAN/authors/id/F/FA/FANGLY/Statistics-R-0.32.tar.gz
wget http://search.cpan.org/CPAN/authors/id/R/RO/ROBM/Cache-FastMmap-1.40.tar.gz
wget http://search.cpan.org/CPAN/authors/id/N/NW/NWCLARK/Devel-Size-0.79.tar.gz
wget http://search.cpan.org/CPAN/authors/id/P/PE/PEREINAR/File-Which-0.05.tar.gz
wget http://search.cpan.org/CPAN/authors/id/L/LE/LEONT/Module-Build-0.4205.tar.gz

## Define modules here in the MODULE array
MODULE=()


MODULE[${#MODULE[@]}]=XML-Simple-2.20
MODULE[${#MODULE[@]}]=Statistics-R-0.32
MODULE[${#MODULE[@]}]=Cache-FastMmap-1.40
MODULE[${#MODULE[@]}]=Devel-Size-0.79
MODULE[${#MODULE[@]}]=File-Which-0.05
MODULE[${#MODULE[@]}]=Module-Build-0.4205

## Loop the MODULE array and install modules.
for element in "${MODULE[@]}"; do

    echo "Installing module $element"
	tar -xvf $element.tar.gz
	cd $element
	perl Makefile.PL PREFIX=$INSTALLPATH
	make
	make install
	cd ..
 
done

cd ..

chmod -R 775 $MUGQIC_INSTALL_HOME/software/perl
