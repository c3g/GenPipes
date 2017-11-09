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
DBD::mysql \
DBD::SQLite \
DBI \
DBI::SQL::Nano \
Devel::Size \
Error \
File::Slurp \
File::Which \
File::ShareDir \
File::Spec::Link \
Filesys::Df \
forks \
GD \
Graph \
Graph::Directed \
GraphViz \
HTML::TableExtract \
Inline \
IO::Scalar \
Parallel::ForkManager \
Parse::Range \
PDF::API2 \
PDF::Create \
Perl::Unsafe::Signals \
PerlIO::gzip \
PostScript::TextBlock \
Proc::ParallelLoop \
SOAP::Lite \
Sort::Naturally \
Spreadsheet::ParseExcel \
Statistics::Descriptive \
Statistics::R \
String::Approx \
SVG \
SVG::Graph \
Test::CPAN::Changes \
Test::CheckManifest \
Test::CPAN::Meta \
Test::CPAN::Meta::JSON \
Test::Pod \
Test::Pod::Coverage \
Test::TrailingSpace \
Text::CSV \
Text::CSV_XS \
Text::CSV::Encoded \
Text::NSP::Measures::2D::Fisher::twotailed \
RPC::PlClient \
Time::HiRes \
XML::DOM \
XML::DOM::XPath \
XML::LibXML \
XML::Parser::PerlSAX \
XML::SAX::Writer \
XML::Simple \
XML::Twig \
XML::Writer \
YAML \
DB_File \
Set::IntervalTree \
; do
$PERL_HOME/bin/perl -MCPAN -e"CPAN::Shell->force(qw(install $MODULE))"
# Test if module is properly installed
$PERL_HOME/bin/perl -e "use $MODULE"
done

# Add permissions
chmod -R ug+rwX,o+rX-w $PERL_HOME
