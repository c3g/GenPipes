#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use PDF::Create;

my $usage=<<'ENDHERE';
NAME:
itaggerTxtToPdf.pl

PURPOSE:
Convert a text file to pdf file

INPUT:
--infile <string>   : text file
--deal_with_headers : Optional. Fix to deal with date prepended to each lines.

OUTPUT:
--outfile <string>  : pdf file.

NOTES:
TODO: Find a way to accurately display tab indentation
in the final pdf file. Find a way to calculate space
would be a good start.

BUGS/LIMITATIONS:

AUTHOR/SUPPORT:
Julien Tremblay - jtremblay@lbl.gov

ENDHERE

## OPTIONS
my ($help, $infile, $outfile, $deal_with_date_headers);
my $verbose = 0;

## SCRIPTS
GetOptions(
    'infile=s' 			=> \$infile,
    'outfile=s' 		=> \$outfile,
    'verbose' 			=> \$verbose,
	'deal_with_headers' => \$deal_with_date_headers,
    'help' 				=> \$help
);
if ($help) { print $usage; exit; }

## MAIN
my $pdf = new PDF::Create('filename' => $outfile,
                          'Version'  => 1.7,
                          'PageMode' => 'UseOutlines',
                          'Author'   => 'Fabien Tassin',
                          'Title'    => 'My title',
                     );
my $root = $pdf->new_page('MediaBox' => [ 0, 0, 612, 792 ]);


# Prepare 2 fonts
my $f1 = $pdf->font('Subtype'  => 'Type1',
                    'Encoding' => 'WinAnsiEncoding',
                    'BaseFont' => 'Helvetica');
my $f2 = $pdf->font('Subtype'  => 'Type1',
                    'Encoding' => 'WinAnsiEncoding',
                    'BaseFont' => 'Courier');

open(IN, $infile) or die "Can't open file ".$infile."\n";
my @new_strings;
my $WIDTH = 120; 
my $barcode_flag = 0;

my @barcode_prefix;
my @barcode_names;
my @barcode_seq;
my @barcode_abs;
my @barcode_rel;

if($deal_with_date_headers){
	while(<IN>){
		chomp($_);
		
		if($_ =~ m/<BARCODES DISTRIBUTION>/){
			$barcode_flag = 1;
			print "BARCODES DISTRIBUTION\n";
	
			@barcode_prefix = ();
			@barcode_names = ();
			@barcode_seq = ();
			@barcode_abs = ();
			@barcode_rel = ();
			
			my @barcodes;
			my @chars = split(//, $_);
			my $i=0;
			my $string = "";
			if(@chars > $WIDTH){
				foreach my $char (@chars){
					if($i % $WIDTH == 0 && $i != 0){
						push(@new_strings, $string);
						$string = $char;
					}else{
						$string .= $char;
					}
					$i++;	
				}
				push(@new_strings, $string) if(length($string) < $WIDTH);
				
			}else{
				my $new_string = join("", @chars);
				push (@new_strings, $new_string);
			}
			
			next;
		}
	
		if($_ =~ m/<\\BARCODES DISTRIBUTION>/){
			$barcode_flag = 0;
			# PRINT BARCODES COUNT
			# find max length for barcode name
			my $max_barcode_prefix_length = 0;
			my $max_barcode_name_length = 0;
			my $max_barcode_seq_length = 0;
			my $max_barcode_abs_length = 0;
			my $max_barcode_rel_length = 0;
	
			foreach my $len (@barcode_prefix){
				my @el = split(//, $len);
				my $el = @el;
				$max_barcode_prefix_length = $el if($el > $max_barcode_prefix_length);
			}
			foreach my $len (@barcode_names){
				my @el = split(//, $len);
				my $el = @el;
				$max_barcode_name_length = $el if($el > $max_barcode_name_length);
			}
			foreach my $len (@barcode_seq){
				my @el = split(//, $len);
				my $el = @el;
				$max_barcode_seq_length = $el if($el > $max_barcode_seq_length);
			}
			foreach my $len (@barcode_abs){
				my @el = split(//, $len);
				my $el = @el;
				$max_barcode_abs_length = $el if($el > $max_barcode_abs_length);
			}
			foreach my $len (@barcode_rel){
				my @el = split(//, $len);
				my $el = @el;
				$max_barcode_rel_length = $el if($el > $max_barcode_rel_length);
			}
	
			# PRINT VALUES RESPECTING TAB VALUES	
			for(my $j=0; $j<@barcode_prefix;$j++){
				my $k;
				my $diff;
				my $new_string = "";
				
				my @chars_p = split(//, $barcode_prefix[$j]);
				my @chars_n = split(//, $barcode_names[$j]);
				my @chars_s = split(//, $barcode_seq[$j]);
				my @chars_a = split(//, $barcode_abs[$j]);
				my @chars_r = split(//, $barcode_rel[$j]);
			
				#=========				
				$k=0;
				foreach my $char (@chars_p){
					$new_string .= $char;
					$k++;
				}
				$diff = $max_barcode_prefix_length - $k;
				for(my $n=0;$n<$diff;$n++){
					$new_string .= " ";
				}
				$new_string .= "      ";
	
				#=========				
				$k=0;
				foreach my $char (@chars_n){
					$new_string .= $char;
					$k++;
				}
				$diff = $max_barcode_name_length - $k;
				for(my $n=0;$n<$diff;$n++){
					$new_string .= " ";
				}
				$new_string .= "      ";
	
				#=========				
				$k=0;
				foreach my $char (@chars_s){
					$new_string .= $char;
					$k++;
				}
				$diff = $max_barcode_seq_length - $k;
				for(my $n=0;$n<$diff;$n++){
					$new_string .= " ";
				}
				$new_string .= "      ";
				
				#=========				
				$k=0;
				foreach my $char (@chars_a){
					$new_string .= $char;
					$k++;
				}
				$diff = $max_barcode_abs_length - $k;
				for(my $n=0;$n<$diff;$n++){
					$new_string .= " ";
				}
				$new_string .= "      ";
	
				#=========			
				$k=0;
				foreach my $char (@chars_r){
					$new_string .= $char;
					$k++;
				}
				$diff = $max_barcode_rel_length - $k;
				for(my $n=0;$n<$diff;$n++){
					$new_string .= "      ";
				}
				
				push (@new_strings, $new_string);	
				#print $new_string."\n";	
			}
	
		}	
	
		## BEHAVE DEPENDING ON FLAG VALUES
		if($barcode_flag == 0){
			my @barcodes;
			my @chars = split(//, $_);
			my $i=0;
			my $string = "";
	
			if(@chars > $WIDTH){
				foreach my $char (@chars){
					if($i % $WIDTH == 0 && $i != 0){
						push(@new_strings, $string);
						$string = $char;
					}else{
						$string .= $char;
					}
					$i++;	
				}
				push(@new_strings, $string) if(length($string) < $WIDTH);
				
			}else{
				my $new_string = join("", @chars);
				push (@new_strings, $new_string);
			}
	
		}elsif($barcode_flag == 1){
			my $i=0;
			my $string = "";
			my @barcodes_values = split(/\s/, $_); #[2][3][4][5]
			
			push(@barcode_prefix, $barcodes_values[0]." ".$barcodes_values[1]);
			push(@barcode_names, $barcodes_values[2]);
			push(@barcode_seq, $barcodes_values[3]);
			push(@barcode_abs, $barcodes_values[4]);
			push(@barcode_rel, $barcodes_values[5]);
		}
	}
	close(IN);
}else{
	while(<IN>){
		chomp($_);
		push(@new_strings, $_);
	}
}

my $x = 20;
my $y = 10;
my $line_counter = 0;
my @pages = ();
my $page_number = 0;
my $page;

foreach my $new_string (@new_strings){
	#print $new_string."\n";
	if($line_counter == 0){
		$page = $root->new_page;	
		push(@pages, $page);
	}elsif($line_counter % 77 == 0 && $line_counter != 0){
		$page = $root->new_page;		
		push(@pages, $page);
		$page_number++;
		$x = 20;
		$y = 10;
		$line_counter = 0;
	}	
	$new_string =~ s/\t/    /g;
	$pages[$page_number]->stringl($f2, 8, $x, 792 - $y, $new_string);
	#$x = $x + 20;
	$y = $y + 10;
	$line_counter++;
}

$pdf->close;
exit;
