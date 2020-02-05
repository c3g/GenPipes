#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

# Perl module where to install libs
#PERL_MODULE=mugqic/perl/5.10.1
#PERL_MODULE=mugqic/perl/5.18.2
PERL_MODULE=mugqic/perl/5.28.2
module load $PERL_MODULE
echo "Perl module loaded :" $PERL_MODULE
# Install Perl modules from CPAN

# MakeMaker's prompt function will always return the default without waiting for user input.
export PERL_MM_USE_DEFAULT=1
# Install module prerequisites automatically
export PERL_EXTUTILS_AUTOINSTALL="--defaultdeps"
echo "starting for loop..."
echo for MODULE in \
Algorithm::Munkres \
Array::Compare \
Bio::ASN1::EntrezGene \
Bio::DB::Fasta \
Bio::DB::SeqFeature::Store \
Bio::Phylo \
Bio::PhyloNetwork \
Bit::Vector \
Carp \
Clone \
Config::General \
Config::Simple \
Convert::Binary::C \
Data::Dumper \
DBD::Gofer \
DBD::mysql \
DBD::SQLite \
DBI \
DBI::SQL::Nano \
Devel::Size \
Digest::MD5 \
Error \
File::Basename \
File::Slurp \
File::Which \
File::ShareDir \
File::Spec::Functions \
File::Spec::Link \
Filesys::Df \
FindBin \
Font::TTF::Font \
forks \
GD \
GD::Polyline \
Getopt::Long \
Graph \
Graph::Directed \
Graphics::ColorObject \
GraphViz \
HTML::TableExtract \
Inline \
IO::File \
IO::Scalar \
List::MoreUtils \
List::Util \
Math::Bezier \
Math::BigFloat \
Math::Round \
Math::VecStat \
Memoize \
Parallel::ForkManager \
Params::Validate \
Parse::Range \
PDF::API2 \
PDF::Create \
Params::Validate \
Perl::Unsafe::Signals \
PerlIO::gzip \
Pod::Usage \
POSIX \
PostScript::TextBlock \
Proc::ParallelLoop \
Readonly \
Regexp::Common \
Set::IntSpan \
SOAP::Lite \
Sort::Naturally \
Spreadsheet::ParseExcel \
Statistics::Descriptive \
Statistics::R \
String::Approx \
Set::IntSpan \
Statistics::Basic \
Storable \
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
Text::Format \
Text::NSP::Measures::2D::Fisher::twotailed \
Time::HiRes \
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
echo perl -MCPAN -e "CPAN::Shell->force(qw(install $MODULE))"
# Test if module is properly installed
$PERL_HOME/bin/perl -e "use $MODULE"
done
echo "for loop finished...."
# Add permissions
chmod -R ug+rwX,o+rX-w $PERL_HOME
