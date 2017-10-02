#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

# Perl module where to install libs
#PERL_MODULE=mugqic/perl/5.10.1
#PERL_MODULE=mugqic/perl/5.18.2
PERL_MODULE=mugqic/perl/5.22.1
module load $PERL_MODULE

# Install Perl modules from CPAN

# MakeMaker's prompt function will always return the default without waiting for user input.
export PERL_MM_USE_DEFAULT=1
# Install module prerequisites automatically
export PERL_EXTUTILS_AUTOINSTALL="--defaultdeps"

for MODULE in \
WWW::Curl \
; do
$PERL_HOME/bin/perl -MCPAN -e"CPAN::Shell->force(qw(install $MODULE))"
# Test if module is properly installed
$PERL_HOME/bin/perl -e "use $MODULE"
done

# Add permissions
chmod -R ug+rwX,o+rX-w $PERL_HOME
