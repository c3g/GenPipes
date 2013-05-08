#!/usr/bin/perl  

use LWP::Simple;
$acc_list=@ARGV[0];
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

print STDERR $url . "\n";

#post the esearch URL
$output = get($url);


#parse WebEnv and QueryKey
$web = $1 if ($output =~ /<WebEnv>(\S+)<\/WebEnv>/);
$key = $1 if ($output =~ /<QueryKey>(\d+)<\/QueryKey>/);
$count = $1 if ($output =~ /<Count>(\d+)<\/Count>/);

#assemble the efetch URL
#$url = $base . "efetch.fcgi?db=nucleotide&query_key=$key&WebEnv=$web";
#$url .= "&rettype=fasta&retmode=text";
#print STDERR $url . "\n" ;
#post the efetch URL
#$fasta = get($url);
#print "$fasta";

#retrieve data in batches of 500
$retmax = 500;
for ($retstart = 0; $retstart < $count; $retstart += $retmax) {
        $efetch_url = $base ."efetch.fcgi?db=nucleotide&WebEnv=$web";
        $efetch_url .= "&query_key=$key&retstart=$retstart";
        $efetch_url .= "&retmax=$retmax&rettype=fasta&retmode=text";
        $efetch_out = get($efetch_url);
	if ($efetch_out =~ /<ERROR>(\S+)<\/ERROR>/) {
   		print STDERR "Error in fetching results for $efetch_url : $errorMessage";
	}else{
	        print "$efetch_out";
	}
	print STDERR $efetch_url . "\n" ;
}
