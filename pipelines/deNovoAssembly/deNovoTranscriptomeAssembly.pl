#!/usr/bin/perl

=head1 NAME

I<deNovoTranscriptomeAssemly>

=head1 SYNOPSIS

deNovoTranscriptomeAssemly.pl B<args> [-f -c -n -s -e]

=head1 DESCRIPTION

B<deNovoTranscriptomeAssemly> Is the main de novo RNA assembly pipeline.

=head1 AUTHORS

B<David Morais> - I<dmorais@cs.bris.ac.uk>

B<Mathieu Bourgey> - I<mbourgey@genomequebec.com>

=head1 DEPENDENCY

B<Pod::Usage> Usage and help output.

B<Data::Dumper> Used to debbug

B<Config::Simple> Used to parse config file

B<File::Basename> path parsing

B<Cwd> path parsing



B<The pipeline specific libs should be in a dir called lib/ placed in the same dir this as this script.>

I<Pipeline specific Libs:>

B<GetFastaAlias>

B<MergeFastq>

B<LoadConfig>

B<SampleSheet>

B<BAMtools>

B<SplitFile>

B<Trinity>

B<BLAST>

B<SampleSheet>

B<SequenceDictionaryParser>

B<SubmitToCluster>

B<Trimmomatic>

B<BWA>

B<ReadStats>

B<HtseqCount>

B<DiffExpression>


=head1 Scripts

This pipeline uses a set of scripts. 
You should create a dir called B<script>
and place it in the same dir as the 
pipeline.

B<This is the list of scripts>

AliasFastqMerge.py

create_gtf.sh

fastalength

fastasplit

filter-duplicates-1.0-jar-with-dependencies.jar

generate_BLAST_HQ.sh

getStat.sh

HtSeq_full_matrix.sh

HtSeq_temp_matrix.sh

ParallelBlast.pl

Parallelize

=cut

# Strict Pragmas
#---------------------
use strict;
use warnings;

#---------------------

BEGIN {
    #Makesure we can find the GetConfig::LoadModules module relative to this script install
    use File::Basename;
    use Cwd 'abs_path';

    my ( undef, $mod_path, undef ) = fileparse( abs_path(__FILE__) );
    unshift @INC, $mod_path . "lib";

}

# Dependencies
#-----------------
use Getopt::Std;
use File::Basename;
use GetFastaAlias;
use MergeFastq;
use LoadConfig;
use Data::Dumper;
use SampleSheet;
use BAMtools;
use SplitFile;
use Trinity;
use BLAST;
use SampleSheet;
use SequenceDictionaryParser;
use SubmitToCluster;
use Trimmomatic;
use BWA;
use ReadStats;
use HtseqCount;
use DiffExpression;

# Globals
#---------------------------------
my @steps;
push( @steps, { 'name' => 'trimming' } );
push( @steps, { 'name' => 'removeDuplicateReads' } );
push( @steps, { 'name' => 'deNovoAssembly' } );
push( @steps, { 'name' => 'blastContig' } );
push( @steps, { 'name' => 'abundance' } );
push( @steps, { 'name' => 'mergeCounts' } );
push( @steps, { 'name' => 'diffExpression' } );

my $merJobId = undef;
my %groupCounter;    # counts the number of sample in each group

# These variables allows control over group loop. In other words
# a group job will be runned only one time
#-----------------------------------------------------
my %groupDone;
my %blastGroupDone;
my %indexGroupDone;
my %diffExpressGroupDone;
my %htSeqGroupDone;
my %readStatsGroupDone;
my @database;
my $step1;

&main();

# SUB
#-----------------
sub printUsage {
    print "\nUsage: perl " . $0 . " args [-f -c -s -e -n]\n";
    print "\t-h  help and usage\n";
    print "\t-f  file Alias\n";
    print "\t-c  config file\n";
    print "\t-s  start step, inclusive\n";
    print "\t-e  end step, inclusive\n";
    print "\t-n  nanuq sample sheet\n";
    print "\n";
    print "Steps:\n";

    for ( my $idx = 0 ; $idx < @steps ; $idx++ ) {
        print "" . ( $idx + 1 ) . '- ' . $steps[$idx]->{'name'} . "\n";
    }
    print "\n";
}

sub main {
    my %opts;
    getopts( 'f:c:s:e:n:', \%opts );

    if ( defined( $opts{'h'} ) || !defined( $opts{'c'} ) || !defined( $opts{'s'} ) || !defined( $opts{'e'} ) || !defined( $opts{'f'} ) || !defined( $opts{'n'} ) ) {
        printUsage();
        exit(1);
    }

    my %cfg               = LoadConfig->readConfigFile( $opts{'c'} );
    my $rHoAoH_sampleInfo = SampleSheet::parseSampleSheetAsHash( $opts{'n'} );
    my ($rH_aliasSampleInfo, $rHoH_groupInfo) = GetFastaAlias::sampleInfo( $opts{'f'}, \%cfg, $rHoAoH_sampleInfo );

    # If there are more than one database
    @database = split /\s+/, $cfg{'blast.db'};

    # Get RunType
    #------------
    my $runType;
    for my $sampleName ( keys %{$rHoAoH_sampleInfo} ) {
        my $rAoH_sampleLanes = $rHoAoH_sampleInfo->{$sampleName};
        for my $rH_laneInfo (@$rAoH_sampleLanes) {
            $rH_laneInfo->{'group_name'} = $rH_aliasSampleInfo->{$sampleName}{'group_name'};
            $runType = $rH_laneInfo->{'runType'};
            SubmitToCluster::initSubmit( \%cfg, $sampleName );
        }
    }
    
    # Merge Step
    #------------
    my $rH_mergeDetails = MergeFastq::mergeFiles( \%cfg, $runType, $opts{'f'} );
    if ( $rH_mergeDetails->{'command'} ne "No merge" ) {
        print $rH_mergeDetails->{'command'};
        print "MERGE_JOB_ID=\"\"\n";
        $merJobId = SubmitToCluster::printSubmitCmd( \%cfg, "merge", "", 'MERGE', undef, "step_0", $rH_mergeDetails->{'command'} );
        $merJobId = '$' . $merJobId;
    }

    # add to %groupCounter the number of
    # samples per group and initialize %diffExpressGroupDone
    #--------------------------------------------------------
    %groupCounter = countGroups($rHoH_groupInfo);
    foreach my $key ( keys %groupCounter ) {
        $diffExpressGroupDone{$key} = 1;
        $rHoH_groupInfo->{$key}->{'mergeSamplesDone'} = 0;
    }

    print "GROUP_JOB_IDs=\"\"\n";
    my $counter = $step1 = $opts{'s'};
    my $cnt = 0;

    # Loop through the steps
    #-----------------------
  START:
    my $n = scalar( keys %{$rHoAoH_sampleInfo} );
    foreach my $sampleName ( keys %{$rHoAoH_sampleInfo} ) {
        my $rAoH_sampleLanes = $rHoAoH_sampleInfo->{$sampleName};

        for ( my $current = $counter - 1 ; $current <= ( $opts{'e'} - 1 ) ; $current++ ) {
            my $fname  = $steps[$current]->{'name'};
            my $subref = \&$fname;

            if ( ( $steps[$current]->{'name'} ne 'trimming' && $steps[$current]->{'name'} ne 'removeDuplicateReads' )
                && $n > 0 && $opts{'s'} <= 2 && $cnt == 0 ) {

                next;
            }

            if ( $n == 1 && $steps[$current]->{'name'} eq 'removeDuplicateReads' ) {

                &$subref( $current != ( $opts{'s'} - 1 ), \%cfg, $sampleName, $rAoH_sampleLanes, $rHoH_groupInfo,  $rH_aliasSampleInfo);

                #                print  " ", $sampleName, " ", $fname, "\n";                             # dry run

                $counter = $counter + 2;
                $cnt     = 1;
                $n--;
                goto START;
            }

            else {
                # Tests for the first step in the list. Used for dependencies.
                #--------------------------------------------------------------

                &$subref( $current != ( $opts{'s'} - 1 ), \%cfg, $sampleName, $rAoH_sampleLanes, $rHoH_groupInfo,  $rH_aliasSampleInfo );

                #                 print  " ", $sampleName, " ", $fname, "\n";                               # dry run

            }

        }
        $n--;

    }

}

sub trimming {
    my $depends          = shift;
    my $rH_cfg           = shift;
    my $sampleName       = shift;
    my $rAoH_sampleLanes = shift;
    my $rHoH_groupInfo   = shift;
    my $rH_aliasSampleInfo = shift;
    my $jobDependency    = ( defined $merJobId ) ? $merJobId : undef;

    print "TRIM_JOB_IDS=\"\"\n" unless $step1 > 1;
    for my $rH_laneInfo (@$rAoH_sampleLanes) {
        my $outputDir = 'reads/'.$sampleName;
        print "mkdir -p ".$outputDir."\n";
        my $rH_trimDetails = Trimmomatic::trim( $rH_cfg, $sampleName, $rH_laneInfo, $outputDir);

        my $trimJobId = undef;
        if ( length( $rH_trimDetails->{'command'} ) > 0 ) {

            $trimJobId = SubmitToCluster::printSubmitCmd( $rH_cfg, "trim", $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}, 'TRIM', undef, $sampleName, $rH_trimDetails->{'command'}, LoadConfig::getParam( $rH_cfg, "default", 'sampleOutputRoot' ).'/'.$sampleName);
            $trimJobId = '$' . $trimJobId;
            print 'TRIM_JOB_IDS=${TRIM_JOB_IDS}' . LoadConfig::getParam( $rH_cfg, 'default', 'clusterDependencySep' ) . $trimJobId . "\n\n";

        }
    }
}

sub removeDuplicateReads {
    my $depends          = shift;
    my $rH_cfg           = shift;
    my $sampleName       = shift;
    my $rAoH_sampleLanes = shift;
    my $rHoH_groupInfo   = shift;
    my $rH_aliasSampleInfo = shift;
    my $jobDependency    = undef;

    if ( $depends > 0 ) {
        $jobDependency = '$TRIM_JOB_IDS';
    }

    my $minQuality = $rH_cfg->{'trim.minQuality'};
    my $minLength = $rH_cfg->{'trim.minLength'};

    print "DUP_JOB_IDS=\"\"\n";
    for my $rH_laneInfo (@$rAoH_sampleLanes) {
        my $rH_dupDetails;

        my $workDir = 'reads/'.$sampleName.'/';

        my $inputFastqPair1Name = $workDir . $sampleName . '.t'.$minQuality.'l'.$minLength.'.pair1.fastq.gz';
        my $inputFastqPair2Name = $workDir . $sampleName . '.t'.$minQuality.'l'.$minLength.'.pair2.fastq.gz';
        my $inputFastqSingleName = $workDir . $sampleName . '.t'.$minQuality.'l'.$minLength.'.single1.fastq.gz';
        my $outputPrefix = $workDir . $sampleName;

        if($rH_laneInfo->{'runType'} eq "PAIRED_END") {
          $rH_dupDetails = BAMtools::deleteDuplicates( $rH_cfg, $sampleName, $inputFastqPair1Name, $inputFastqPair2Name, undef, $outputPrefix);
          $rH_aliasSampleInfo->{$sampleName}{'bwa_pair1'}   = $rH_dupDetails->{'pair1'};
          $rH_aliasSampleInfo->{$sampleName}{'bwa_pair2'}   = $rH_dupDetails->{'pair2'};
        } else {
          $rH_dupDetails = BAMtools::deleteDuplicates( $rH_cfg, $sampleName, undef, undef, $inputFastqSingleName, $outputPrefix);
        }
        $rH_aliasSampleInfo->{$sampleName}{'bwa_single1'} = $rH_dupDetails->{'single1'};
        $rH_aliasSampleInfo->{$sampleName}{'bwa_single2'} = $rH_dupDetails->{'single2'};

        # add the returning value to the hash ref $rH_aliasSampleInfo. To be used by BWA aln

        # return 0 if this step was called from BWA step (we are only interested on the $rH_aliasSampleInfo values)
        if ( $step1 > 2 ) {return 0;}

        my $dupJobId = undef;
        if ( length( $rH_dupDetails->{'command'} ) > 0 ) {
            $dupJobId = SubmitToCluster::printSubmitCmd( $rH_cfg, "duplicate", $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}, 'DUP', $jobDependency, $sampleName, $rH_dupDetails->{'command'}, LoadConfig::getParam( $rH_cfg, "default", 'sampleOutputRoot' ).'/'.$sampleName );
            $dupJobId = '$' . $dupJobId;
            print 'DUP_JOB_IDS=${DUP_JOB_IDS}' . LoadConfig::getParam( $rH_cfg, 'default', 'clusterDependencySep' ) . $dupJobId . "\n";
            print 'GROUP_JOB_IDS=${GROUP_JOB_IDS}' . LoadConfig::getParam( $rH_cfg, 'default', 'clusterDependencySep' ) . $dupJobId . "\n\n";
        }
    }

}

sub deNovoAssembly {
    my $depends          = shift;
    my $rH_cfg           = shift;
    my $sampleName       = shift;
    my $rAoH_sampleLanes = shift;
    my $rHoH_groupInfo   = shift;
    my $rH_aliasSampleInfo = shift;
    my $jobDependency    = undef;
    my $group            = $rH_aliasSampleInfo->{$sampleName}{'group_name'};
    if ( $depends > 0 ) {
        $jobDependency = '$GROUP_JOB_IDS';
    }
    
     
    # Chrysalis
    #---------------
    for my $rH_laneInfo (@$rAoH_sampleLanes) {

        #next if the one sample of the group is already done
        next if ( exists $groupDone{$group} );
        print "if [ ! -d " . $group . "/output_jobs ] ; then mkdir -p " . $group . "/output_jobs ; fi\n";
        print "if [ ! -d assembly/" . $group . "/output_jobs ]; then  mkdir -p assembly/" . $group . " ; fi\n";
        
        print "CHRYSALIS_JOB_IDS=\"\"\n";

        my $rH_chrysalisDetails = Trinity::chrysalis( $rH_cfg, $group, $rH_laneInfo, $rHoH_groupInfo->{$group}->{'left'} , $rHoH_groupInfo->{$group}->{'right'}); 
        $groupDone{$group} = 1;
        my $chrysalisJobId = undef;

        if ( length( $rH_chrysalisDetails->{'command'} ) > 0 ) {
            $chrysalisJobId = SubmitToCluster::printSubmitCmd( $rH_cfg, "trinity", $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}, 'CHRYSALIS', $jobDependency, $group, $rH_chrysalisDetails->{'command'} );
            $chrysalisJobId = '$' . $chrysalisJobId;
            print 'CHRYSALIS_JOB_IDS=${CHRYSALIS_JOB_IDS}' . LoadConfig::getParam( $rH_cfg, 'default', 'clusterDependencySep' ) . $chrysalisJobId . "\n\n";
            print "GROUP_JOB_IDS=\"\"\n";
        }

        # Split Chrysalis file
        #---------------------------
        print "SPLITCHRYS_JOB_IDS=\"\"\n";
        my $rH_splitChrysalisDetails = SplitFile::splitButterfly( $rH_cfg, $group, $rH_laneInfo );
        my $splitChrysalisJobId = undef;
        if ( length( $rH_splitChrysalisDetails->{'command'} ) > 0 ) {
            $jobDependency = '$CHRYSALIS_JOB_IDS';
            $splitChrysalisJobId = SubmitToCluster::printSubmitCmd( $rH_cfg, "splitchrysalis", $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}, 'SPLITCHRYSALIS', $jobDependency, $group, $rH_splitChrysalisDetails->{'command'} );
            $splitChrysalisJobId = '$' . $splitChrysalisJobId;
            print 'SPLITCHRYS_JOB_IDS=${SPLITCHRYS_JOB_IDS}' . LoadConfig::getParam( $rH_cfg, 'default', 'clusterDependencySep' ) . $splitChrysalisJobId . "\n\n";
        }

        #ButterFly
        #----------
        # This puts the files names in the @Files array (with a 4 digits padding)
        my @files;
        my $nbChunks = LoadConfig::getParam( $rH_cfg, 'butterfly', 'chunks');
        for ( my $i = 0 ; $i < $nbChunks ; $i++ ) {
            push( @files, 'chunk.' . sprintf( "%04d", $i ) . '.txt') ;
        }

        print "BUTTERFLY_JOB_IDS=\"\"\n";

        foreach my $file (@files) {
            my $rH_butterflyDetails = Trinity::butterfly( $rH_cfg, $group, $rH_laneInfo, $file );
            my $butterflyJobId = undef;
            if ( length( $rH_splitChrysalisDetails->{'command'} ) > 0 ) {
                $jobDependency = '$SPLITCHRYS_JOB_IDS';
                $butterflyJobId = SubmitToCluster::printSubmitCmd( $rH_cfg, "butterfly", $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} . '_'. basename($file), 'BUTTERFLY', $jobDependency, $group, $rH_butterflyDetails->{'command'} );
                $butterflyJobId = '$' . $butterflyJobId;
                print 'BUTTERFLY_JOB_IDS=${BUTTERFLY_JOB_IDS}' . LoadConfig::getParam( $rH_cfg, 'default', 'clusterDependencySep' ) . $butterflyJobId . "\n\n";

            }

        }

        # Concatenate files and create gft
        #-----------------------------------
        print "CONCAT_JOB_IDS=\"\"\n";
        my $rH_concatDetails = Trinity::concatFastaCreateGtf( $rH_cfg, $group, $rH_laneInfo );
        my $concatJobId = undef;
        if ( length( $rH_splitChrysalisDetails->{'command'} ) > 0 ) {
            $jobDependency = '$BUTTERFLY_JOB_IDS';
            $concatJobId = SubmitToCluster::printSubmitCmd( $rH_cfg, "concat", $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}, 'CONCAT', $jobDependency, $group, $rH_concatDetails->{'command'} );
            $concatJobId = '$' . $concatJobId;
            print 'CONCAT_JOB_IDS=${CONCAT_JOB_IDS}' . LoadConfig::getParam( $rH_cfg, 'default', 'clusterDependencySep' ) . $concatJobId . "\n\n";

        }

    }

}

sub blastContig {
    my $depends          = shift;
    my $rH_cfg           = shift;
    my $sampleName       = shift;
    my $rAoH_sampleLanes = shift;
    my $rHoH_groupInfo   = shift;
    my $rH_aliasSampleInfo = shift;
    my $jobDependency    = undef;
    my $group            = $rH_aliasSampleInfo->{$sampleName}{'group_name'};

    my $laneDirectory  = 'assembly/' . $group . "/";
    
    if ( $depends > 0 ) {
        $jobDependency = '$CONCAT_JOB_IDS';
    }

    my $fileName = 'Trinity.2.fasta';

    # Split Multifasta file
    #---------------------------
    for my $rH_laneInfo (@$rAoH_sampleLanes) {

        #next if one sample of the group is already done
        next if ( exists $blastGroupDone{$group} );

        print "SPLITFASTA_JOB_IDS=\"\"\n";
        my $rH_splitFastaDetails = SplitFile::splitFasta( $fileName, $rH_cfg, $group, $rH_laneInfo, $laneDirectory );
        $blastGroupDone{$group} = 1;
        my $splitFastaJobId = undef;
        if ( length( $rH_splitFastaDetails->{'command'} ) > 0 ) {
            $splitFastaJobId = SubmitToCluster::printSubmitCmd( $rH_cfg, "splitfasta", $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}, 'SPLITFASTA', $jobDependency, $group, $rH_splitFastaDetails->{'command'} );
            $splitFastaJobId = '$' . $splitFastaJobId;
            print 'SPLITFASTA_JOB_IDS=${SPLITFASTA_JOB_IDS}' . LoadConfig::getParam( $rH_cfg, 'default', 'clusterDependencySep' ) . $splitFastaJobId . "\n\n";
        }

        #BLAST
        #-------
        print "BLAST_JOB_IDS=\"\"\n";
        $jobDependency = '$SPLITFASTA_JOB_IDS';
        my @files;
        for ( my $i = 1 ; $i <= $rH_cfg->{'blast.chunks'} ; $i++ ) {
            push( @files, $fileName.'_chunk_' . sprintf( "%07d", $i ) );
        }

        foreach my $db (@database) {
            $jobDependency = '$SPLITFASTA_JOB_IDS';
            foreach my $file (@files) {
                my $rO_blastDetailsJob = BLAST::alignParallel( $rH_cfg, $group, $rH_laneInfo, $file, $db, $laneDirectory );
                if ( $rO_blastDetailsJob->isUp2Date() == 0) {
                    SubmitToCluster::printSubmitCmd( $rH_cfg, "blast", $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}, 'BLAST', $jobDependency, $group, $rO_blastDetailsJob);
                    print 'BLAST_JOB_IDS=${BLAST_JOB_IDS}' . LoadConfig::getParam( $rH_cfg, 'default', 'clusterDependencySep' ) . $rO_blastDetailsJob->getCommandJobId(0) . "\n\n";

                }

            }

            # BLAST BEST HIT
            #----------------
            $jobDependency = '$BLAST_JOB_IDS';
            print "BLASTBESTHIT_JOB_IDS=\"\"\n";
            my $rO_blastBestHitDetailsJob = BLAST::bestHit( $rH_cfg, $group, $rH_laneInfo, $db, $laneDirectory );
            my $blastbesthitJobId = undef;
            if ( $rO_blastBestHitDetailsJob->isUp2Date() == 0) {
                SubmitToCluster::printSubmitCmd( $rH_cfg, "blastbesthit", $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}, 'BLASTBESTHIT', $jobDependency, $group, $rO_blastBestHitDetailsJob);
                print 'BLASTBESTHIT_JOB_ID=${BLASTBESTHIT_JOB_IDS}' . LoadConfig::getParam( $rH_cfg, 'default', 'clusterDependencySep' ) . $rO_blastBestHitDetailsJob->getCommandJobId(0) . "\n\n";
            }
        }
    }
}

sub abundance {
    my $depends          = shift;
    my $rH_cfg           = shift;
    my $sampleName       = shift;
    my $rAoH_sampleLanes = shift;
    my $rHoH_groupInfo   = shift;
    my $rH_aliasSampleInfo = shift;
    my $jobDependency    = undef;
    my $group            = $rH_aliasSampleInfo->{$sampleName}{'group_name'};

    if ( $depends > 0 ) {
        $jobDependency = '$CONCAT_JOB_IDS';
    }

    # Step needed to get the pair and single names. It only executes this if
    # the script did not start from the removeDuplicateReads step
    if ( $step1 != 1 ) {
        for my $rH_laneInfo (@$rAoH_sampleLanes) {
            removeDuplicateReads( 0, $rH_cfg, $sampleName, $rAoH_sampleLanes, $rHoH_groupInfo, $rH_aliasSampleInfo);
        }
    }

    print 'mkdir -p ' . "alignment/" . $group . "\n";
    my $assembly = 'assembly/' . $group . '/Trinity.fasta';
    for my $rH_laneInfo (@$rAoH_sampleLanes) {
        my $outputPrefix = "alignment/" . $group.'/'.$sampleName.'.rsem';

        my $command = Trinity::abundance($rH_cfg, $assembly, $outputPrefix, $rH_aliasSampleInfo->{$sampleName}{'bwa_pair1'}, $rH_aliasSampleInfo->{$sampleName}{'bwa_pair2'});
        # One needs to run before others to setup the reference.
        if(defined($rHoH_groupInfo->{$group}->{'abundancePrepRef'})){
          $jobDependency = $rHoH_groupInfo->{$group}->{'abundancePrepRef'};
        }

        my $jobId = SubmitToCluster::printSubmitCmd( $rH_cfg, "abundance", $sampleName, 'ABUNDANCE', $jobDependency, $group, $command);
        if(!defined($rHoH_groupInfo->{$group}->{'abundancePrepRef'})) {
          my $firstJobId = "FIRST_".$jobId;
          print $firstJobId.'='.'$'.$jobId."\n";
          $rHoH_groupInfo->{$group}->{'abundancePrepRef'} = '$'.$firstJobId;
          print "ABUNDANCE_JOB_IDS=\"\"\n";
        }
        print 'ABUNDANCE_JOB_IDS=${ABUNDANCE_JOB_IDS}'.LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep').'$'.$jobId."\n";
    }

    return '$ABUNDANCE_JOB_IDS'
}

sub mergeCounts {
    my $depends          = shift;
    my $rH_cfg           = shift;
    my $sampleName       = shift;
    my $rAoH_sampleLanes = shift;
    my $rHoH_groupInfo   = shift;
    my $rH_aliasSampleInfo = shift;

    my $jobDependency    = undef;
    my $group            = $rH_aliasSampleInfo->{$sampleName}{'group_name'};

    # if this is not the last sample in the group next
    if ( $rHoH_groupInfo->{$group}->{'mergeSamplesDone'} < $groupCounter{$group} ) {
      $rHoH_groupInfo->{$group}->{'mergeSamplesDone'}++;
      return;
    }

    if ( $depends > 0 ) {
        $jobDependency = '$ABUNDANCE_JOB_IDS';
    }

    my %inputs;
    for my $rH_laneInfo (@$rAoH_sampleLanes) {
      $inputs{$rH_laneInfo->{'name'}} = 1;
    }

    my @inputList = keys(%inputs);
    my $output = "bad";
    my $command = Trinity::mergeCounts($rH_cfg, \@inputList, $output);
    my $jobId = undef;
    if(defined($command)) {
      $jobId = SubmitToCluster::printSubmitCmd( $rH_cfg, "mergeCounts", $group, 'MERGECNT', $jobDependency, $group, $command);
      $jobId = '$'. $jobId;
    }

    return $jobId;
}

sub readStats {
    my $depends          = shift;
    my $rH_cfg           = shift;
    my $sampleName       = shift;
    my $rAoH_sampleLanes = shift;
    my $rHoH_groupInfo   = shift;
    my $rH_aliasSampleInfo = shift;
    my $jobDependency    = undef;
    my $group            = $rH_aliasSampleInfo->{$sampleName}{'group_name'};

    if ( $depends > 0 ) {
        $jobDependency = '$BWA_JOB_IDS';
    }
    print "READSTATS_JOB_IDS=\"\"\n";
    for my $rH_laneInfo (@$rAoH_sampleLanes) {

        my $read = 'reads/' . $sampleName . '.t' . $rH_cfg->{'trim.minQuality'} . 'l' . $rH_cfg->{'trim.minLength'} . '.pair1.fastq.gz';
        my $sortedBam = 'alignment/' . $group . '/' . $sampleName . '.sorted.bam';
        my $rH_readStatDetails = ReadStats::stats( $rH_cfg, $sampleName, $rH_laneInfo, $read, $sortedBam, $group );
        my $readStatJobId = undef;

        my $rH_readStatConcatDetails = ReadStats::concatStats( $rH_cfg, $sampleName, $rH_laneInfo, $group );
        my $readStatConcatJobId = undef;

        if ( length( $rH_readStatDetails->{'command'} ) > 0 ) {
            $readStatJobId = SubmitToCluster::printSubmitCmd( $rH_cfg, "readstats", $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}, 'READSTATS', $jobDependency, $group, $rH_readStatDetails->{'command'} );
            $readStatJobId = '$' . $readStatJobId;
            print 'READSTATS_JOB_IDS=${READSTATS_JOB_IDS}' . LoadConfig::getParam( $rH_cfg, 'default', 'clusterDependencySep' ) . $readStatJobId . "\n\n";

        }

        # Concatenate
        #------------
        $jobDependency = '$READSTATS_JOB_IDS';
        print "READSTATSCONCAT_JOB_IDS=\"\"\n";
        if ( length( $rH_readStatConcatDetails->{'command'} ) > 0 ) {
            $readStatConcatJobId = SubmitToCluster::printSubmitCmd( $rH_cfg, "readstatsconcat", $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}, 'READSTATSCONCAT', $jobDependency, $group, $rH_readStatConcatDetails->{'command'} );
            $readStatConcatJobId = '$' . $readStatConcatJobId;
            print 'READSTATSCONCAT_JOB_IDS=${READSTATSCONCAT_JOB_IDS}' . LoadConfig::getParam( $rH_cfg, 'default', 'clusterDependencySep' ) . $readStatConcatJobId . "\n\n";
        }

        # Contig Stats
        #------------

        $jobDependency = '$CONCAT_JOB_IDS';

        #next if one sample of the group is already done
        next if ( exists $readStatsGroupDone{$group} );

        my $rH_readStatContigDetails = ReadStats::contigStats( $rH_cfg, $sampleName, $rH_laneInfo, $group );
        $readStatsGroupDone{$group} = 1;
        my $readStatContigJobID = undef;

        print "READSTATCONTIG_JOB_IDS=\"\"\n";
        if ( length( $rH_readStatContigDetails->{'command'} ) > 0 ) {
            $readStatContigJobID = SubmitToCluster::printSubmitCmd( $rH_cfg, "readstatcontig", $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}, 'READSTATCONTIG', $jobDependency, $group, $rH_readStatContigDetails->{'command'} );
            $readStatContigJobID = '$' . $readStatContigJobID;
            print 'READSTATCONTIG_JOB_IDS=${READSTATCONTIG_JOB_IDS}' . LoadConfig::getParam( $rH_cfg, 'default', 'clusterDependencySep' ) . $readStatContigJobID . "\n\n";
        }

    }

}

sub htseqCount {

    my $depends          = shift;
    my $rH_cfg           = shift;
    my $sampleName       = shift;
    my $rAoH_sampleLanes = shift;
    my $rHoH_groupInfo   = shift;
    my $rH_aliasSampleInfo = shift;
    my $jobDependency    = undef;
    my $group            = $rH_aliasSampleInfo->{$sampleName}{'group_name'};

    if ( $depends > 0 ) {
        $jobDependency = '$READSTATS_JOB_IDS';
    }

    no warnings;    # disable the warning created by the next
                    #next if one sample of the group is already done
    next if ( exists $htSeqGroupDone{$group} );

    # Picard sorting
    #----------------
    print "HTSEQSORT_JOB_IDS=\"\"\n";
    for my $rH_laneInfo (@$rAoH_sampleLanes) {

        my $rH_htseqSortDetails = HtseqCount::sortRead( $rH_cfg, $sampleName, $rH_laneInfo, $group );
        $htSeqGroupDone{$group} = 1;
        my $htseqSortJobId = undef;

        if ( length( $rH_htseqSortDetails->{'command'} ) > 0 ) {
            $htseqSortJobId = SubmitToCluster::printSubmitCmd( $rH_cfg, "htseqsort", $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}, 'HTSEQSORT', $jobDependency, $group, $rH_htseqSortDetails->{'command'} );
            $htseqSortJobId = '$' . $htseqSortJobId;
            print 'HTSEQSORT_JOB_IDS=${HTSEQSORT_JOB_IDS}' . LoadConfig::getParam( $rH_cfg, 'default', 'clusterDependencySep' ) . $htseqSortJobId . "\n\n";
        }

    }

    # Htseq-count
    #-------------
    $jobDependency = '$HTSEQSORT_JOB_IDS';
    print "HTSEQCOUNT_JOB_IDS=\"\"\n";
    for my $rH_laneInfo (@$rAoH_sampleLanes) {

        my $rH_htseqCountDetails = HtseqCount::readCount( $rH_cfg, $sampleName, $rH_laneInfo, $group );
        my $htseqCountJobId = undef;

        if ( length( $rH_htseqCountDetails->{'command'} ) > 0 ) {
            $htseqCountJobId = SubmitToCluster::printSubmitCmd( $rH_cfg, "htseqcount", $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}, 'HTSEQCOUNT', $jobDependency, $sampleName, $rH_htseqCountDetails->{'command'} );
            $htseqCountJobId = '$' . $htseqCountJobId;
            print 'HTSEQCOUNT_JOB_IDS=${HTSEQCOUNT_JOB_IDS}' . LoadConfig::getParam( $rH_cfg, 'default', 'clusterDependencySep' ) . $htseqCountJobId . "\n\n";

        }

    }

    # Htseq-matrixMake
    #------------------
    $jobDependency = '$HTSEQCOUNT_JOB_IDS';
    print "HTSEQMATRIX_JOB_IDS=\"\"\n";
    foreach my $db (@database) {
        for my $rH_laneInfo (@$rAoH_sampleLanes) {

            my $rH_htseqMatrixDetails = HtseqCount::matrixMake( $rH_cfg, $sampleName, $rH_laneInfo, $db, $group );
            my $htseqMatrixJobId = undef;

            if ( length( $rH_htseqMatrixDetails->{'command'} ) > 0 ) {
                $htseqMatrixJobId = SubmitToCluster::printSubmitCmd( $rH_cfg, "htseqmatrix", $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}, 'HTSEQMATRIX', $jobDependency, $sampleName, $rH_htseqMatrixDetails->{'command'} );
                $htseqMatrixJobId = '$' . $htseqMatrixJobId;
                print 'HTSEQMATRIX_JOB_IDS=${HTSEQMATRIX_JOB_ID}' . LoadConfig::getParam( $rH_cfg, 'default', 'clusterDependencySep' ) . $htseqMatrixJobId . "\n\n";
            }

        }

    }

}

sub diffExpression {
    my $depends          = shift;
    my $rH_cfg           = shift;
    my $sampleName       = shift;
    my $rAoH_sampleLanes = shift;
    my $rHoH_groupInfo   = shift;
    my $rH_aliasSampleInfo = shift;
    my $jobDependency    = undef;
    my $group            = $rH_aliasSampleInfo->{$sampleName}{'group_name'};
    
    my $outputDir = 'DGE/' . $group . "/";
    my $matrix = $outputDir . 'matrix.csv';
    my $desingFile = $rH_cfg->{'diffExpress.designFile'};
    if ( $depends > 0 ) {
        $jobDependency = '$HTSEQMATRIX_JOB_IDS';
    }

    # if this is not the last sample in the group next
    if ( $diffExpressGroupDone{$group} < $groupCounter{$group} ) {

        $diffExpressGroupDone{$group}++;
        return;
    }

    print "DIFFEXPRESS_JOB_IDS=\"\"\n";
    for my $rH_laneInfo (@$rAoH_sampleLanes) {

        my $rH_diffExpresstDetails = DiffExpression::edger( $rH_cfg, $desingFile, $matrix, $outputDir );
        my $diffExpressJobId = undef;

        if ( length( $rH_diffExpresstDetails->{'command'} ) > 0 ) {
            $diffExpressJobId = SubmitToCluster::printSubmitCmd( $rH_cfg, "diffexpress", $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}, 'DIFFEXPRESS', $jobDependency, $sampleName, $rH_diffExpresstDetails->{'command'} );
            $diffExpressJobId = '$' . $diffExpressJobId;
            print 'DIFFEXPRESS_JOB_IDS=${DIFFEXPRESS_JOB_IDS}' . LoadConfig::getParam( $rH_cfg, 'default', 'clusterDependencySep' ) . $diffExpressJobId . "\n\n";

        }

    }
}

sub countGroups {
    # Counts the number of samples per group
    my $rHoH_group = shift;
    my %group_count;

    foreach my $group ( keys(%{ $rHoH_group }) ) {
        my @aux = split(/\s+/, $rHoH_group->{$group}->{'left'});
        $group_count{$group} = scalar(@aux) - 1;
    }

    return %group_count;
}

