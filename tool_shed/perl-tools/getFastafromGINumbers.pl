#!/usr/bin/perl

use LWP::Simple;
#$gi_list = '24475906,224465210,50978625,9507198';
$gi_list=@ARGV[0];

#assemble the URL
$base = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
$url = $base . "efetch.fcgi?db=nucleotide&id=$gi_list&rettype=acc";

print STDERR $url ."\n";

#post the URL
$output = get($url);

@accessionNumber = split(/\n/,$output);
$acc_list= join(",",@accessionNumber);
@acc_array = split(/,/, $acc_list);

#append [accn] field to each accession
for ($i=0; $i < @acc_array; $i++) {
   $acc_array[$i] .= "[accn]";
}

#join the accessions with OR
$query = join('+OR+',@acc_array);

#assemble the esearch URL
$base = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
$url = $base . "esearch.fcgi?db=nucleotide&term=$query&usehistory=y";

print STDERR $url."\n" ;

#post the esearch URL
$output = get($url);


#parse WebEnv and QueryKey
$web = $1 if ($output =~ /<WebEnv>(\S+)<\/WebEnv>/);
$key = $1 if ($output =~ /<QueryKey>(\d+)<\/QueryKey>/);

#assemble the efetch URL
$url = $base . "efetch.fcgi?db=nucleotide&query_key=$key&WebEnv=$web";
$url .= "&rettype=fasta&retmode=text";

print STDERR $url."\n" ;

#post the efetch URL
$fasta = get($url);
print "$fasta";
