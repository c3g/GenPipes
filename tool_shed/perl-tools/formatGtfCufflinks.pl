# Formats a denovo cufflinks merged.gtf to get rid of redundant info
# Maxime Caron - Jan 2012
#Mathieu Bourgey - Jan 2013

use Switch;

print "Input: $ARGV[0]";
print "\n********************\n";
print "Output: $ARGV[1]\n";

open(INFO, $ARGV[0]);
open(OUT, ">".$ARGV[1]);

while(<INFO>) {

	chomp($_);
	if($_ =~ "exon_number \"1\"") {
	@splitA = split(/\t/, $_);
	@splitB = split(";", $splitA[8]);
	$sizeB = @splitB;
	my $splitOUT= $splitA[0] ."\t".$splitA[1] ."\t".$splitA[2] ."\t".$splitA[3] ."\t".$splitA[4] ."\t".$splitA[5] ."\t".$splitA[6];

	
	
	if($sizeB == 8) {
		@splitCode = split(/\s+/, $splitB[6]);
		$splitCode[2] =~ s/"//g;
		$code = &classify_code($splitCode[2]);
		print OUT $splitOUT ."\t". $splitB[1] ."\t".$splitB[4]."\t".$splitB[5]."\tinfo \"".$code." \"\n";	
	}
	elsif($sizeB == 7) {
		@splitCode = split(/\s+/, $splitB[5]);
		$splitCode[2] =~ s/"//g;
		$code = &classify_code($splitCode[2]);
		print OUT $splitOUT ."\t".$splitB[1] ."\t".$splitB[3]."\t".$splitB[4]."\tinfo \"".$code." \"\n";
	}
	elsif($sizeB == 6) {
		@splitCode = split(/\s+/, $splitB[4]);
		$splitCode[2] =~ s/"//g;
		$code = &classify_code($splitCode[2]);
                print OUT $splitOUT ."\t".$splitB[1] ."\t".$splitB[3]."\t - - \tinfo \"".$code." \"\n";
        }
	else {
		@splitCode = split(/\s+/, $splitB[7]);
		$splitCode[2] =~ s/"//g;
		$code = &classify_code($splitCode[2]);
		$splitB[7] = &classify_code($splitB[7]);
		print OUT $splitOUT ."\t".$splitB[1] ."\t".$splitB[4]."\t".$splitB[6]."\tinfo \"".$code." \"\n";
	}
}

}
close(OUT); 
close(INFO);  

sub classify_code
{
	switch ($_[0]) { 
			case '=' { return 'Complete match of intron chain' }
			case 'c' { return 'Contained' }
			case 'j' { return 'Potentially novel isoform (fragment): at least one splice junction is shared with a reference transcript' }
			case 'e' { return 'Single exon transfrag overlapping a reference exon and at least 10 bp of a reference intron, indicating a possible pre-mRNA fragment' }
			case 'i' { return 'A transfrag falling entirely within a reference intron' }
                        case 'o' { return 'Generic exonic overlap with a reference transcript' }
                        case 'p' { return 'Possible polymerase run-on fragment (within 2Kbases of a reference transcript)' }
                        case 'r' { return 'Repeat. Currently determined by looking at the soft-masked reference sequence and applied to transcripts where at least 50% of the bases are lower case' }
			case 'u' { return 'Unknown, intergenic transcript' }
                        case 'x' { return 'Exonic overlap with reference on the opposite strand' }
                        case 's' { return 'An intron of the transfrag overlaps a reference intron on the opposite strand (likely due to read mapping errors)' }
                        case '.' { return '(.tracking file only, indicates multiple classifications)' }
		}
}	
