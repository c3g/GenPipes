#!/usr/bin/perl

use strict;
use POSIX qw( strftime );
use Net::SMTP;
use XML::Simple;
use Data::Dumper;
use Number::Format;

use POSIX qw(strftime);

#use File::stat;
use Getopt::Long;
use File::Basename;

my $version = "0.1";

my $EMAIL_SENT_FLAG_FILE = "nanuqRunIndex-email.empty";    # file written when the "operation" is completed
my $MAX_UNEXPECTED_TO_DISPLAY = 15;
my $MIN_READ_COUNT = 100;#0000;


&main();

sub main {
   my @runDirectories;
   my @to;
   my $verbose;
   my $fromEmail;
   my $nanuqAuthFile;
   my $printFullSummary;
   my $result = GetOptions(
      "from=s"   => \$fromEmail,
      "to=s"     => \@to,
      "runDir=s" => \@runDirectories,
      "verbose"  => \$verbose,
      "fullSummary"  => \$printFullSummary,
      "nanuqAuthFilei=s" => \$nanuqAuthFile,
   );
   
   $printFullSummary = $printFullSummary || $verbose;
  
   if(!defined($nanuqAuthFile) || !-e $nanuqAuthFile) {
     die "Missing nanuqAuthFile\n";
   }   
   
   my @summaryWaitingForData;
   my @summaryEmailAlreadySent;
   my @summaryWaitingForJobs;
   my @summaryJobSubmited;
   my @summaryNoIndexed;
   my @summaryAnalysisError;
   my @summaryAnalysisSuccess;

   for my $rootDir (@runDirectories) {
      if ( !-d $rootDir ) {
         warn "'" . $rootDir . "' is not a directory. Ignored\n";
      }

      opendir( ROOT_DIR, $rootDir ) or die "Couldn't open directory " . $rootDir . "\n";
      my @rootFiles = grep { /^[^\.]/ && -d "$rootDir/$_" } readdir(ROOT_DIR);
      closedir(ROOT_DIR);

      for my $partialRunDir (@rootFiles) {
         my $runDir = $rootDir . '/' . $partialRunDir;
         if ( !-d $runDir ) {
            next;
         }
         
         print STDOUT "Processing run '$partialRunDir'\n" if ($verbose);
         
         if (!(( $partialRunDir =~ /.*_\d+HS\d\d[AB]/ ) || ( $partialRunDir =~ /.*\d+_[^_]+_\d+_.+/ ))) {
            print STDOUT "\tNot a valid run folder: Skipping\n" if ($verbose);
            next;
         }
         
         if ( !-e $runDir . "/RunInfo.xml" ) {
            print STDOUT "\tNo RunInfo.xml: Skipping\n" if ($verbose);
            next;
         }

         my $xml                = new XML::Simple( 'ForceArray' => ['Read'] );
         my $runInfo            = $xml->XMLin( $runDir . "/RunInfo.xml" );
         my $maxIndexReadNumber = 0;
         my $readNb             = 1;
         my $isIndexed;
         my $indexLength = 0;
         my @readsInfo;
         my $nbLanes = $runInfo->{'Run'}->{'FlowcellLayout'}->{'LaneCount'};
         for my $rH_readInfo ( @{ $runInfo->{'Run'}->{'Reads'}->{'Read'} } ) {
            my $nbCycles = $rH_readInfo->{'NumCycles'};
            my $isCurrentReadIndexed = $rH_readInfo->{'IsIndexedRead'};
            push(@readsInfo, {nbCycles=>$nbCycles , isIndexed=>$isCurrentReadIndexed});            
            if ($rH_readInfo->{'IsIndexedRead'} eq 'Y'){
               $maxIndexReadNumber = $readNb;
               $indexLength += $nbCycles;
               $isIndexed = "Y";
            }
            $readNb++;
         }

         my $fileToLookForIndex = "Basecalling_Netcopy_complete_Read" . $maxIndexReadNumber . ".txt";
         print STDOUT "\tLooking for file '$fileToLookForIndex'\n" if ($verbose);
         
         opendir( RUN_DIR, $runDir ) or die "Couldn't open possible run directory " . $runDir . "\n";
         my @runFiles = grep { /^[^\.]/ && -f "$runDir/$_" } readdir(RUN_DIR);
         closedir(RUN_DIR);

         my $indexReady       = 0;
         my $indexEmailSent   = 0;
         my $metricsJobCount  = 0;
         my $metricsFileCount = 0;
         for my $partialRunFile (@runFiles) {
            if ( ( $isIndexed eq "Y" ) && ( $partialRunFile eq $fileToLookForIndex ) ) {
               $indexReady = 1;
            } elsif ( $partialRunFile eq $EMAIL_SENT_FLAG_FILE ) {
               $indexEmailSent = 1;
            } elsif ( $partialRunFile =~ /^nanuqRunIndex\-.+_\d+\.metrics$/ ) {
               $metricsFileCount++;
            } elsif ( $partialRunFile =~ /^nanuqRunIndex-jobStarted_\d+\.empty$/ ) {
               $metricsJobCount++;
            }
         }
         

         if ( $isIndexed ne "Y" ) {
            print STDOUT "\tRun is not indexed\n" if ($verbose);
            push (@summaryNoIndexed, $partialRunDir);
         } elsif (!$indexReady) {
            print STDOUT "\tWaiting for data\n" if ($verbose);
            push (@summaryWaitingForData, $partialRunDir);
         } elsif ($indexEmailSent) {
            print STDOUT "\tEmail already sent\n" if ($verbose);
            push (@summaryEmailAlreadySent, $partialRunDir);            
         } elsif ( $metricsJobCount == 0 ) {
            @readsInfo = splice(@readsInfo,0,$maxIndexReadNumber);
            print STDOUT "\tLaunching jobs...\n" if ($verbose);
            launchJobs ($runDir, $indexLength, \@readsInfo, $nbLanes, $verbose, $fromEmail, $nanuqAuthFile);
            print STDOUT "\tJob Launching Completed\n" if ($verbose);
            push (@summaryJobSubmited, $partialRunDir);
         } elsif ( $metricsJobCount == $metricsFileCount ) {
            # validate index and send email
            print STDOUT "\tAnalysing results\n" if ($verbose);
            launchJobs ($runDir, $indexLength, \@readsInfo, $nbLanes, $verbose, $fromEmail, $nanuqAuthFile, 1);
            if (analyzeRun( $runDir, $indexLength, $fromEmail, \@to, $verbose )) {
               push (@summaryAnalysisError, $partialRunDir);
            } else {
               push (@summaryAnalysisSuccess, $partialRunDir);
            }
         } else {
            print STDOUT "\tWaiting for metrics files jobs to finish\n" if ($verbose);
            push (@summaryWaitingForJobs, $partialRunDir);
         }
      }
   }

   #### Summary
   if (1) {
      my $displaydate= strftime('%Y-%m-%d %H:%M:%S', localtime(time));
      print STDOUT "\n##################### SUMMARY #####################";
      print STDOUT "\nTime: $displaydate";
      print STDOUT "\nNon-Indexed:\n\t" . join("\n\t", sort @summaryNoIndexed) if (scalar(@summaryNoIndexed) && $printFullSummary);
      print STDOUT "\nWaiting data:\n\t" . join("\n\t", sort @summaryWaitingForData) if (scalar(@summaryWaitingForData) && $printFullSummary);
      print STDOUT "\nEmail already sent:\n\t" . join("\n\t", sort @summaryEmailAlreadySent) if (scalar(@summaryEmailAlreadySent) && $printFullSummary);
      print STDOUT "\nWaiting job completion:\n\t" . join("\n\t", sort @summaryWaitingForJobs) if (scalar(@summaryWaitingForJobs));
      print STDOUT "\nJob submitted:\n\t" . join("\n\t", sort @summaryJobSubmited) if (scalar(@summaryJobSubmited));
      print STDOUT "\nAnalysis Error:\n\t" . join("\n\t", sort @summaryAnalysisError) if (scalar(@summaryAnalysisError));
      print STDOUT "\nAnalysis Success:\n\t" . join("\n\t", sort @summaryAnalysisSuccess) if (scalar(@summaryAnalysisSuccess));
      print STDOUT "\n";
   }
   
}

##################
# Job Generation #
##################

sub launchJobs {
  my $runDirectory = shift;
  my $indexLength = shift;
  my $rAoH_readsInfo = shift;
  my $nbLanes = shift;
  my $verbose = shift;
  my $email = shift;
  my $nanuqAuthFile = shift;
  my $onlySampleSheet = shift;
  
  my ($runName, $runID) = getRunName($runDirectory);

  downloadSampleSheet($runDirectory, $nanuqAuthFile, $runID, $verbose);
  generateIndexCounts($email, $runDirectory, $runID, $runName, $nbLanes, $indexLength, $rAoH_readsInfo, $verbose) if (!$onlySampleSheet);
  
}

sub getRunName {
  my $runDirectory = shift;
  my $runID;
  my $runName;
  if ( $runDirectory =~ /.*_\d+HS\d\d[AB]/ ) {
    ($runName,$runID) = $runDirectory =~ /.*\/(\d+_[^_]+_\d+_[^_]+_(\d+)HS.+)/;
  }
  elsif ( $runDirectory =~ /.*\d+_[^_]+_\d+_.+/ ) {
    ($runName,$runID) = $runDirectory =~ /.*\/(\d+_([^_]+_\d+)_.*)/;
  }
  else {
    die "Not a valid run folder: " . $runDirectory . "\n";
  }
  return ($runName, $runID);
}

sub generateIndexCounts {
  my $email           = shift;
  my $runDirectory    = shift;
  my $runID           = shift;
  my $runName         = shift;
  my $nbLanes         = shift;
  my $indexLength     = shift;
  my $rAoH_readsInfo  = shift;
  my $verbose         = shift;

  my $baseCallDir     = $runDirectory.'/Data/Intensities/BaseCalls';
  
  my $mask = "";
  my $indexPrinted=0;
  for my $rH_readInfo (@{$rAoH_readsInfo}) {
    if($rH_readInfo->{'isIndexed'} eq 'Y') {
      if($indexPrinted == 0) {
        $mask .= $indexLength.'B';
        $indexPrinted=1;
      } 
    }
    else {
      $mask .= $rH_readInfo->{'nbCycles'}.'T';
    }
  }

  for(my $lane=1; $lane <= $nbLanes; $lane++) {
    my $command = "touch $runDirectory/nanuqRunIndex-jobStarted_$lane.empty";
    print STDOUT "\t\t$command\n" if ($verbose); 
    system($command);
    my $stdOutRedirect = ($verbose) ? "" : " > /dev/null";
    $command = 'echo "java -Xmx55G -jar ~/CountIlluminaBarcodes-0.2-alpha-jar-with-dependencies.jar BASECALLS_DIR='.$baseCallDir.' LANE='.$lane.' READ_STRUCTURE='.$mask.' METRICS_FILE='.$runDirectory.'/nanuqRunIndex-'.$runName.'_'.$lane.'.metrics MAX_MISMATCHES=1 NUM_PROCESSORS=12 BARCODE_FILE=~/barcodes.txt" | msub -d '.$runDirectory.' -V -l nodes=1:ppn=12 -l walltime=12:00:0 -q sw -j oe -N idx_'.$runID.'_'.$lane.' -m ae -M '.$email. $stdOutRedirect;
    print STDOUT "\t\t$command\n" if ($verbose); 
    system($command);
  }
}

sub downloadSampleSheet {
  my $runDirectory   = shift;
  my $nanuqAuthFile  = shift;
  my $runID          = shift;
  my $verbose        = shift;

  if ( !-e $runDirectory . "/SampleSheet.nanuq.csv" ) {
    my $instrument = 'HiSeq';
    if($runDirectory =~ /_M00/){
      $instrument = 'MiSeq';
    }
    my $stdOutRedirect = ($verbose) ? "" : " > /dev/null";
    my $command ='wget --no-cookies --directory-prefix '.$runDirectory.'/ --post-file '.$nanuqAuthFile.' https://genomequebec.mcgill.ca/nanuqMPS/sampleSheet/'.$instrument.'/'.$runID.'/SampleSheet.nanuq.csv'.$stdOutRedirect;
    print STDOUT "\t\t$command\n" if ($verbose);
    system($command);
    if ($? == -1) {
      print "failed to execute: $!\n";
      #exit(1);
    }
    elsif ($? & 127) {
      #printf "child died with signal %d, %s coredump\n", ($? & 127), ($? & 128) ? 'with' : 'without';
      #exit(1);
    }
    else {
      my $childValue = $? >> 8;
      if($childValue != 0) {
        #printf "child exited with value %d\n", $childValue;
        #exit(1);
      }
    }
  }
}



###############################################
# Result Analysis
###############################################

sub analyzeRun {
   my $runDirectory = shift;
   my $indexLength  = shift;
   my $fromEmail    = shift;
   my $rA_to        = shift;
   my $verbose      = shift;

   print STDOUT "\t\tParsing sample sheet: '$runDirectory/SampleSheet.nanuq.csv' ... " if ($verbose);
   my $rHoA_sampleSheetInfo = parseSampleSheet("$runDirectory/SampleSheet.nanuq.csv");
   if ($rHoA_sampleSheetInfo) {
      print STDOUT "Success\n" if ($verbose);
   } else {
      print STDOUT "Failed\n" if ($verbose);
      return 666;
   }
   
   my $text         = "Run: $runDirectory\n\n";
   my $nbFound      = 0;
   my $nbMissing    = 0;
   my $nbUnexpected = 0;

   my $formater = new Number::Format;
   my $format = "#,###,###,###";
   my $percentFormat = "###.#%";

   # Process each lane of the sample sheet ordered by number
   for my $laneNumber ( sort keys %$rHoA_sampleSheetInfo ) {
      my $rA_expectedBarcodes = $rHoA_sampleSheetInfo->{$laneNumber};

      my $indexRemaining      = scalar(@$rA_expectedBarcodes);

      my @foundList;         # barcode in sample sheet and metrics file
      my @unexpectedList;    # barcode in metrics file but not in the sample sheet
      my @rawUnexpectedList;
      my @missingList;       # barcode in sample sheet and metrics file (but with a low score)
      my @foundIndexes;      # list of barcode found, used to determine barcodes found in sample sheet but not in metrics file (not even in low score entries)
      my @foundNames;

      my $previousCount   = 0;
      my $thresholdPassed = 0;
      my $total;
      my $rAoH_metrics;
      my $currentIndex = "";
      
      if (scalar(@$rA_expectedBarcodes) == 1 && @$rA_expectedBarcodes[0] eq "NoIndex") {
         ($total, $rAoH_metrics) = parseMetricsFile( $runDirectory . "/nanuqRunIndex-" . basename($runDirectory) . "_" . $laneNumber . ".metrics", 6, 1);
         push (@foundList, { "index" => "NoIndex", "names" => "", "reads" => $total });
         push (@foundIndexes, "NoIndex");
      } else {
         my $laneIndexLength = 0;
         for my $barcode (@$rA_expectedBarcodes) {
            if (length($barcode) > $indexLength) {
               $barcode = substr($barcode, 0, $indexLength);
            }
            if (length($barcode) > $laneIndexLength) {
               $laneIndexLength = length($barcode);
            }
         }         
         
         ($total, $rAoH_metrics) = parseMetricsFile( $runDirectory . "/nanuqRunIndex-" . basename($runDirectory) . "_" . $laneNumber . ".metrics", $laneIndexLength );
         print STDOUT "\t\tParsing metrics file: '" . $runDirectory . "/nanuqRunIndex-" . basename($runDirectory) . "_" . $laneNumber . ".metrics' ... " if ($verbose);
         
         if ( scalar(@$rAoH_metrics) == 0 || $total == 0) {
            print STDOUT "Invalid metrics file for lane $laneNumber in run : '$runDirectory'\n" if ($verbose);
            return -1;
         } else {
            print STDOUT "Success (" . scalar(@$rAoH_metrics) . " found)\n" if ($verbose);
         }         
         for my $rH_currentMetric (@$rAoH_metrics) {
            # "valid" indexes are those that have at least 1/5 of the read count of the previous entry
            if (( $$rH_currentMetric{"reads"} < ( $previousCount / 5 ) ) || ($indexRemaining == 0 && ($$rH_currentMetric{"reads"} < ( $previousCount / 2 )) )) {
               $thresholdPassed = 1;
            }
   
            $currentIndex = $$rH_currentMetric{'index'};
            
            if (length($currentIndex) > $laneIndexLength) {
               $currentIndex = substr($currentIndex, 0, $laneIndexLength);
               $$rH_currentMetric{'index'} = $currentIndex;
            }
            
            my $currentNames = $$rH_currentMetric{'names'};
            if ( grep( /^\Q$currentIndex\E$/, @$rA_expectedBarcodes ) ) {
   
               # we have an exact match
               $indexRemaining--;
               push( @foundIndexes, $currentIndex );
               push( @foundNames, split(",", $currentNames));
               if ($thresholdPassed) {
                  # too bad, the current index doesn't have enough reads
                  push( @missingList, $rH_currentMetric );
               } else {
                  # yeah!
                  push( @foundList, $rH_currentMetric );
               }
            } else {
   
               my $found = 0;
               
               if ($currentNames ne "") {
                  my @currentNames = split(",", $currentNames);
                  for my $foundEntry (@foundList) {
                     my @candidateNames = split("," , $$foundEntry{'names'});
                     my %names;
                     @names{@currentNames};
                     @names{@candidateNames};
                     
                     if ((scalar(keys(%names)) < (scalar(@candidateNames) + scalar(@currentNames))) && scalar(keys(%names)) > 0) {
                        $$foundEntry{'reads'} += $$rH_currentMetric{"reads"};
                        $found = 1;
                        last;
                     }
                  }
                  
                  if (!$found) {
                     for my $unexpectedEntry (@unexpectedList) {
                        my @candidateNames = split("," , $$unexpectedEntry{'names'});
                        my %names;
                        @names{@currentNames};
                        @names{@candidateNames};
                        if ((scalar(keys(%names)) < (scalar(@candidateNames) + scalar(@currentNames))) && scalar(keys(%names)) > 0) {
                           $$unexpectedEntry{'reads'} += $$rH_currentMetric{"reads"};
                           $found = 1;
                           last;
                        }
                     }
                  }
   
                  if (!$found) {
                     for my $missingEntry (@missingList) {
                        my @candidateNames = split("," , $$missingEntry{'names'});
                        my %names;
                        @names{@currentNames};
                        @names{@candidateNames};
                        
                        if ((scalar(keys(%names)) < (scalar(@candidateNames) + scalar(@currentNames))) && scalar(keys(%names)) > 0) {
                           $$missingEntry{'reads'} += $$rH_currentMetric{"reads"};
                           $found = 1;
                           last;
                        }
                     }
                  }               
               }
               if ($found) {
                  next;
               }
               if (!$thresholdPassed) {
                  push( @unexpectedList, $rH_currentMetric );
                  push( @foundNames, split(",", $currentNames));
               }
            }
   
            $previousCount = $$rH_currentMetric{"reads"};
   
         }
      }
      

      print STDOUT "\t\tPreparing the email content for lane $laneNumber ... " if ($verbose);
      

      # Determine barcodes found in sample sheet but not in metrics file (not even in low score entries)
      my $dels = join '|', map quotemeta, @foundIndexes;
      my @missingIndexes = grep !/$dels/, @$rA_expectedBarcodes;

      $nbFound      += scalar(@foundList);
      $nbMissing    += scalar(@missingList) + scalar(@missingIndexes);
      $nbUnexpected += scalar(@unexpectedList);

      # Message construction...
      $text .= "Lane $laneNumber\n";
      $text .= "\tOK: \n";
      for my $H_currentMetrics (@foundList) {
         $text .= "\t\t" . $formater->format_picture($$H_currentMetrics{'reads'}, $format) . " ". $formater->format_picture($$H_currentMetrics{'reads'}/$total*100, $percentFormat) . "\t$$H_currentMetrics{'index'}\t$$H_currentMetrics{'names'}\n";
      }
      $text .= "\tMissing in the run data:\n" if ( scalar(@missingList) + scalar(@missingIndexes) > 0 );
      for my $H_currentMetrics (@missingList) {
         $text .= "\t\t" . $formater->format_picture($$H_currentMetrics{'reads'}, $format) . " ". $formater->format_picture($$H_currentMetrics{'reads'}/$total*100, $percentFormat) . "\t$$H_currentMetrics{'index'}\t$$H_currentMetrics{'names'}\n";
      }
      for my $rH_currentMetrics (@missingIndexes) {
         $text .= "\t\t" . $formater->format_picture(0, $format) . " ". $formater->format_picture(0, $percentFormat) . "\t" . $rH_currentMetrics . "\n";
      }
      $text .= "\tUnexpected (found in data, but not in nanuq):\n " if ( scalar(@unexpectedList) > 0 );
      my $count = 0;
      for my $H_currentMetrics (@unexpectedList) {
         if ( $count++ < $MAX_UNEXPECTED_TO_DISPLAY ) {
            $text .= "\t\t" . $formater->format_picture($$H_currentMetrics{'reads'}, $format) . " ". $formater->format_picture($$H_currentMetrics{'reads'}/$total*100, $percentFormat) . "\t$$H_currentMetrics{'index'}\t$$H_currentMetrics{'names'}\n";
         } else {
            $text .= "\t\t" . scalar( @unexpectedList - $MAX_UNEXPECTED_TO_DISPLAY ) . " more...\n";
            last;
         }
      }
      $text .= "\tTotal:\t". $formater->format_picture($total, $format) . "\n";      
      $text .= "\n";

      print STDOUT "\t\tSuccess\n" if ($verbose);
   }

   # Send email and create a flag file
   my $errorCode = 0;
   print STDOUT "\t\tSending email ... " if ($verbose);
   $errorCode = sendEmail( $fromEmail, $rA_to, $runDirectory, $text, $nbFound, $nbMissing, $nbUnexpected );
   print STDOUT "\t\tSent with code $errorCode \n" if ($verbose);
   if ( $errorCode == 0 ) {

      my $touchFile = "$runDirectory/$EMAIL_SENT_FLAG_FILE";
      open( TOUCH, ">", $touchFile ) or warn "\t\tCannot touch: " . $touchFile . "\n";
      close(TOUCH);
      if ($verbose) {
         print STDOUT "\t\tTouched: " . $touchFile . "\n";
      }
   }

   return $errorCode;
}

sub parseMetricsFile {
   my $file = shift;
   my $indexLength = shift;
   my $countAllReads = shift; # not only PF Reads
   my @AoH_sequences;
   my $total = 0;

   open my $info, $file or die "Could not open '$file': $!";
   while ( my $line = <$info> ) {

      # Sequence   Names   Reads   pfReads
      if ( $line =~ /(.*?)\t(.*?)\t(\d+)\t(\d+)\s*/g ) {
         my $sequence = $1;
         my $names    = $2;
         my $count    = ($countAllReads) ? $3 : $4;
         $names =~ s/\s*//g;

         $total += $count;
         
         if ( length($sequence) > $indexLength) {
            $sequence = substr($sequence, 0, $indexLength);
            
            my $found = 0;
            for my $rH_oldSequences (@AoH_sequences) {
               if ($$rH_oldSequences{'index'} eq $sequence) {
                  $$rH_oldSequences{'reads'} += $count;
                  $found = 1;
                  last;
               }
            }
            next if ($found);
         }
         
         if ( $sequence =~ /^\.+$/ ) {
            # The barcode sequence of sample in a non-indexed lane is "......"
            $sequence = "NoIndex";
         }
         
         if ($count >= $MIN_READ_COUNT) {
            push( @AoH_sequences, { "index" => $sequence, "names" => $names, "reads" => $count } );
         }
      }
   }
   close($info);

   #sort in descending order
   @AoH_sequences =  sort { $b->{reads} <=> $a->{reads} } @AoH_sequences;
   
   for my $rH_source (@AoH_sequences){
      my $currentIndex = $$rH_source{"index"};
      if ($currentIndex =~ /\./) {
         for my $rH_candidate (@AoH_sequences) {
            my $candidateIndex = $$rH_candidate{"index"};
            if ($candidateIndex =~ /\./) {
               next;
            }
            if ($candidateIndex =~ /^$currentIndex$/) {
               $$rH_candidate{'reads'} += $$rH_source{"reads"};
               $$rH_source{'reads'} = 0;
               last;
            }
         }
      }   
   }
   @AoH_sequences =  sort { $b->{reads} <=> $a->{reads} } @AoH_sequences;

   return ($total, \@AoH_sequences);
}

sub parseSampleSheet {
   my $file = shift;
   my %HoA_map;

   open my $info, $file or return;
   my $line = <$info>;
   while ( my $line = <$info> ) {
      my @tokens = split( ",", $line );
      if ( scalar(@tokens) > 4 ) {
         my $lane = $tokens[1];
         my $sequence = Trim( $tokens[4] );
         $sequence =~ s/\-//;    # remove nextera "-" between the 2 indexes
         if ( $sequence eq "" ) {

            # The barcode sequence of sample in a non-indexed lane
            $sequence = "NoIndex";
         }
         if ( defined( $HoA_map{ $lane } ) ) {
            my $rA_array = $HoA_map{ $lane };
            push( @$rA_array, $sequence );
         } else {
            my @array = ($sequence);
            $HoA_map{ $lane } = \@array;
         }

      }
   }
   return \%HoA_map;
}

sub sendEmail {
   my $fromEmail    = shift;
   my $rA_to        = shift;
   my $run          = shift;
   my $text         = shift;
   my $nbFound      = shift;
   my $nbMissing    = shift;
   my $nbUnexpected = shift;

   my $smtp = Net::SMTP->new('192.168.32.1');
   if ( !defined($smtp) || !($smtp) ) {
      print STDOUT "SMTP ERROR: Unable to open smtp session.\n";
      return 1;
   }

   if ( !( $smtp->mail($fromEmail) ) ) {
      print STDOUT "SMTP ERROR: Cannot set from: " . $fromEmail . ".\n";
      return 2;
   }

   if ( !( $smtp->to(@$rA_to) ) ) {
      print STDOUT "SMTP ERROR: Cannot set to: " . join( ',', @$rA_to ) . ".\n";
      return 3;
   }

   my $tos = join( ',', @$rA_to );
   
   my ($runName, $runID) = getRunName($run);

   my $subject = "Subject: Index Validation - $runID - (OK:$nbFound - Missing:$nbMissing - Unexpected:$nbUnexpected)\n";

   $smtp->data();
   $smtp->datasend( "To: " . $tos . "\n" );
   $smtp->datasend($subject);
   $smtp->datasend("\n");
   $smtp->datasend($text);
   $smtp->dataend();

   $smtp->quit;

   return 0;
}

sub Trim {
   my $text = shift;
   $text =~ s/^\s+//;
   $text =~ s/\s+$//;
   return $text;
}

