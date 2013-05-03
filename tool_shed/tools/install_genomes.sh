#/bash/src	
######## Installing a genome from a source path, accesion number or GI numbers list

## The purpose of this script is to provide a template for standardized installation of new genomes.
## Specially for new organisms / small genome size / small lists of NCBI GI/accession numbers
## TO DO : Output to bash script, containing submissions to the moab scheduler

# USAGE : . install_genomes.sh genomes_list.tsv
# genomes_list.tsv IS A TAB SEPARATED FILE CONTINING 4 COLUMNS
# species,  group,  source , type, alias
# if type is zip or tar.gz, source must be the url link to the compressed fasta file(s)
# if type is gi or acc, source must contain a list (comma separated) of GI or accession numbers respectively
# if type is fasta, source must be the path to a fasta file(s) (file will be copied to the $MUGQIC_INSTALL_HOME/genome directory)

# TODO: create a subdirectory for ech version of samtools installed in the 

# Examples of genomes_list.tsv file lines:
#potato	NUCLEAR http://solanaceae.plantbiology.msu.edu/data/PGSC_DM_v3_2.1.11_pseudomolecules.zip	zip     S_tuberosum_Group_Phureja_DM1-3_Assembly_Version_3_DM_PGSC_Version_2.1.11_Pseudomolecule_Sequences
#maize	CHLOROPLAST	11994090	gi	Zea_mays_chloroplast_complete_genome
#carrot	SELECTED_GENES	DQ192201,DQ192196,DQ192186,DQ192193,DQ192191,DQ222430,AY368232,DQ192189,DQ192197,DQ192188,DQ192187,DQ192192,EU664593	acc	Daucus_carota_subsp_Sativus_selected_genes
#maize	NUCLEAR_MAIZE_SEQ_ORG	http://ftp.maizesequence.org/current/assembly/ZmB73_RefGen_v2.tar.gz	tar.gz	Maize_Genome_Sequencing_Project_RefGen_v2
#potato	NUCLEAR	/sb/programs/analyste/genomes/potato/tmp	fasta	complete_PGSC_DM_genome_v2.1.11_plus_superscaffolds_v3_2.1.9


function main {
  InstallGenomes $1
}

function epilogue {
	###################
	################### Epilogue
	###################
	# Set the MUGQIC_INSTALL_HOME environment variable
	# Load MUGQIC module files
	umask 0002
	## mugic module
	HOST=`hostname`;
	DNSDOMAIN=`dnsdomainname`;
	if [[ $HOST == abacus* ]]; then
	export MUGQIC_INSTALL_HOME=/sb/programs/analyste
	elif [[ $HOST == lg-* ]]; then
	export MUGQIC_INSTALL_HOME=/software/areas/genomics
	elif [[ $HOST == ip03 ]]; then
 	export MUGQIC_INSTALL_HOME=/mnt/lustre03/bourque/bourque_group/opt
	elif [[ $DNSDOMAIN == ferrier.genome.mcgill.ca ]]; then
	export MUGQIC_INSTALL_HOME=/sb/programs/analyste
	elif [[ $DNSDOMAIN == guillimin.clumeq.ca ]]; then
	export MUGQIC_INSTALL_HOME=/software/areas/genomics
	elif [[ $DNSDOMAIN == m ]]; then
	export MUGQIC_INSTALL_HOME=/mnt/lustre03/bourque/bourque_group/opt
	fi
	module use $MUGQIC_INSTALL_HOME/modulefiles
}

function isAvailable {
	###################
	################### isAvailable
	###################
	# Check if a module is available, matching by the module name, 
	# return the output of module avail command, extracting the latest version if many versions are available 
	moduleName=$1;
	available=`module avail -l 2>&1 | grep $moduleName | sort -t "/" -k 3 -gr | sed -e 's/^.*\/*\('$moduleName'[/].*\) .$/\1/g' | head -1 | awk '{print $1}'`;
	echo $available;
}


function loadModules {
	###################
	################### load Modules needed to manage fasta files
	###################
	# Check if a module is available, matching by the module name, 
	# return the output of module avail command, extracting the latest version if many versions are available 
	module use $MUGQIC_INSTALL_HOME/modulefiles
	for mymodule in bwa samtools picard "mugqic\/tools" "mugqic\/R";
	do
		moduleLoad=$( isAvailable $mymodule )
		if  [[ ! -z "$moduleLoad" ]]; then
				echo "Loading module " $moduleLoad
				module load $moduleLoad;
		else 
				echo "ERROR: Module " $mymodules " is not available";
				exit 1;				
		fi;
	done;
}


function InstallGenome () {
	##
	######### Installing a genome from a source path, accesion number of GI numbers list (tab separated)
	# POTATO	NUCLEAR	http://solanaceae.plantbiology.msu.edu/data/PGSC_DM_v3_2.1.11_pseudomolecules.zip	zip	S_tuberosum_Group_Phureja_DM1-3_Assembly_Version_3_DM_PGSC_Version_2.1.11_Pseudomolecule_Sequences
	genomeName=$1; 
	genomeGroup=$2;
	genomeSource=$3;
	gType=$4;
	genomeAlias=$5;
	fastaFiles="";

	## Fetch from source
	INSTALL_PATH=$MUGQIC_INSTALL_HOME/genomes/$genomeName/$genomeGroup
	mkdir -p $INSTALL_PATH
	cd $INSTALL_PATH
	if [[ -z "$genomeAlias" ]]; then
    genomeAlias=`echo $genomeSource | awk -F"/" '{print $NF}'`
	fi
	case $gType in
      "zip")  
        #if [ ! -f $genomeAlias.zip ]; then
					wget $genomeSource -O $genomeAlias.zip
				#fi;
				fastaFiles=`unzip -l $genomeAlias.zip | grep -v zip | grep "\.fa[s]*[t]*[a]*" | awk '{print $NF}'`
				unzip -o $genomeAlias.zip
				cat $fastaFiles > $genomeAlias.fasta
				fastaFiles=`echo $fastaFiles $genomeAlias.fasta`
        ;;
     "tar.gz") 
				#if [ ! -f $genomeAlias.tar.gz ]; then
					wget $genomeSource -O $genomeAlias.tar.gz
				#fi;
 				fastaFiles=`tar -ztvf $genomeAlias.tar.gz | grep "\.fa[s]*[t]*[a]*" | awk '{print $NF}'`
				tar -xzvf $genomeAlias.tar.gz
				# merge all fasta files in a general fasta and add it too the list files to index/faidx,etc
				cat $fastaFiles > $genomeAlias.fasta
				fastaFiles=`echo $fastaFiles $genomeAlias.fasta`
				;;
     "uniq.gz") 
				#if [ ! -f $genomeAlias.tar.gz ]; then
					wget $genomeSource -O $genomeAlias.uniq.gz
				#fi;
 				fastaFiles=`gunzip -l $genomeAlias.uniq.gz | grep "\.uniq" | awk '{print $NF}'`
				gunzip $genomeAlias.tar.gz
				# merge all fasta files in a general fasta and add it too the list files to index/faidx,etc
				cat $fastaFiles > $genomeAlias.fasta
				fastaFiles=`echo $fastaFiles $genomeAlias.fasta`
				;;
				"fasta") 
			  cp $genomeSource $genomeAlias.fasta
			  fastaFiles=$genomeAlias.fasta
        ;;

        *)  
			  echo  "Type of genome file is not yet treated " $gType
			  return 1
        ;;
  esac
  
  mkdir -p fasta;
  for fa in $fastaFiles ;
  do
		rm fasta/$fa;
		mv $fa fasta/$fa;	
		ln -fs fasta/$fa $fa;
    dictName=`echo $fa | sed -e 's/\.fasta/\.dict/g' | sed -e 's/\.fa/\.dict/g'`
    chrsize=`echo $fa | sed -e 's/\.fasta/\.dict/g' | sed -e 's/\.fa/\.chromsize\.txt/g'`
    # bowtie-build Index reference (already done for iGenomes)
		#    bowtie-build /lb/project/mugqic/epigenome/genome_files/hg19/fasta/hg19.fasta /lb/project/mugqic/epigenome/genome_files/hg19/fasta/hg19.fasta # example Maxime
		# Index reference with samtools faidx reference (already done for iGenomes)
		samtools faidx fasta/$fa 	
		# reference dictionary (already done for iGenomes)
		java -jar $PICARD_HOME/CreateSequenceDictionary.jar REFERENCE=fasta/$fa OUTPUT=fasta/$dictName 
		# Generate tab delimited chromosome size file
		sed 1d fasta/$dictName | awk '{print $2"\t"$3}' | sed -e 's/SN://g' | sed -e 's/LN://g' > $chrsize
		# bwa index  reference (already done for iGenomes)
		bwa index -a bwtsw fasta/$fa 
	done;
}

function InstallGenomes () {
	curDir=`pwd`;
  epilogue
  loadModules  
	sourceFile=$1
	while read line
	do
		echo "Installing" $line;
		InstallGenome $line
	  cd $curDir
	done < $sourceFile;
}

main $1;
