#!/bin/bash

# TODO: not sure how GeneMark integrate.. plus anticipate perl version problems
# 2. GeneMark-ES. Download from http://exon.biology.gatech.edu
# 3. FGENESH 2.4 or higher. Purchase from http://www.softberry.com
# 4. GeneMarkS. Download from http://exon.biology.gatech.edu



MODULES="mugqic/blast/2.3.0+ mugqic/exonerate/2.2.0  postgresql/9.5.0 bioinformatics/RepeatMasker/4-0-6 mugqic_dev/augustus/2.7 gcc/4.7.0 openmpi_gcc64/1.6.4" # mugqic/perl/5.18.2  
module load  $MODULES
export PATH="$PATH:/cvmfs/soft.mugqic/root/software/snap/snap-2013-11-29/"
ld_preload_path="/opt/mpi/gcc/openmpi-1.6.4/lib64/libmpi.so"
export LD_PRELOAD="$ld_preload_path"


SOFTWARE="maker"  
VERSION="3.00.0-beta"
ARCHIVE=$SOFTWARE-$VERSION.tgz
ARCHIVE_URL="http://yandell.topaz.genetics.utah.edu/maker_downloads/54D1/F1B6/A363/16BBFFD0E46178E85A1EE9B25BF8/$ARCHIVE" 
SOFTWARE_DIR=$SOFTWARE-$VERSION  ## TO BE MODIFIED WITH SPECIFIC SOFTWARE DIRECTORY IF NECESSARY

INSTALL_DOWNLOAD="$MUGQIC_INSTALL_HOME_DEV/software/maker/maker-$VERSION"
mkdir -p $INSTALL_DOWNLOAD
cd $INSTALL_DOWNLOAD
wget $ARCHIVE_URL

# rm -rf maker
tar zxvf $ARCHIVE 
cd $SOFTWARE
cd src
perl Build.PL
./Build status          
./Build installdeps      
./Build installexes     
./Build install         
./Build status     




cd $INSTALL_DOWNLOAD
MODULEFILE="$MUGQIC_INSTALL_HOME_DEV/modulefiles/mugqic_dev/$SOFTWARE/$VERSION"
mkdir -p `dirname $MODULEFILE`
echo "\
#%Module1.0
proc ModulesHelp { } {
  puts stderr \"\tMUGQIC - $SOFTWARE \"
}
module-whatis \"$SOFTWARE\"

module load $MODULES
prepend-path    PATH                $MUGQIC_INSTALL_HOME/software/snap/snap-2013-11-29/ ; 
set             root                $INSTALL_DOWNLOAD/$SOFTWARE
setenv          LD_PRELOAD          $ld_preload_path
setenv          OMPI_MCA_mpi_warn_on_fork 0
prepend-path    PATH                \$root/bin ;  
" > $MODULEFILE



# echo " for m in gcc/4.7.0 openmpi_gcc64/1.6.4 mugqic_dev/ray/2.3.1; do module load \$m; done && export ENABLE_PARALLEL_LOCAL_SCRATCH=1 && cp lighter_k21/all.cor.fq \$PARALLEL_LOCAL_SCRATCH/ && mpiexec -n 1024 Ray -read-write-checkpoints /mnt/parallel_scratch_mp2_wipe_on_december_2015/bourque/bourque_group/jfliu_swine_pb_illumina_assembly_PRJBFX-1092/ray_checkpoints/41 -k 41 -write-kmers -route-messages -connection-type polytope -routing-graph-degree 16 -s \$PARALLEL_LOCAL_SCRATCH/all.cor.fq -o ray/41 &> ray/41.log " | qsub -N ray_k41 -d /mnt/parallel_scratch_mp2_wipe_on_december_2015/bourque/bourque_group/jfliu_swine_pb_illumina_assembly_PRJBFX-1092 -V -m ae -M francois.lefebvre3@mail.mcgill.ca -j oe -o jobs_output -l nodes=43 -l walltime=120:0:0
# 
# 
# 9:10:34 AM Gary Leveque: $ module load mugqic_dev/maker/3.00.0-beta && mpiexec -n 4 maker  maker_bopts.ctl maker_exe.ctl maker_opts.ctl
# 9:10:52 AM Gary Leveque: mpiexec: /opt/pgi/linux86-64/11.10/libso/libnuma.so.1: no version information available (required by /opt/mpi/gcc/openmpi-1.6.4/lib64/libopen-rte.so.4)
#
# mpiexec: /opt/pgi/linux86-64/11.10/libso/libnuma.so.1: no version information available (required by /opt/mpi/gcc/openmpi-1.6.4/lib64/libopen-rte.so.4)
#
# [ip03:27586] mca: base: component_find: unable to open /opt/mpi/gcc/openmpi-1.6.4/lib64/openmpi/mca_paffinity_hwloc: /opt/mpi/gcc/openmpi-1.6.4/lib64/openmpi/mca_paffinity_hwloc.so: undefined symbol: opal_hwloc_topology (ignored)
#
# [ip03:27589] mca: base: component_find: unable to open /opt/mpi/gcc/openmpi-1.6.4/lib64/openmpi/mca_paffinity_hwloc: /opt/mpi/gcc/openmpi-1.6.4/lib64/openmpi/mca_paffinity_hwloc.so: undefined symbol: opal_hwloc_topology (ignored)    blah
#
# 
# 
# export LD_PRELOAD=Š/openmpi_location/lib/libmpi.so
# mpiexec -mca btl ^openib -n 40 maker


# cd /mnt/parallel_scratch_mp2_wipe_on_april_2017/bourque/bourque_group/srehan_allodapine_genomes_PRJBFX-1096/gapfiller/Exoneurella_tridentata/gapfiller/deconseq/maker
#   module load mugqic_dev/maker/3.00.0-beta && mpiexec -mca btl ^openib -n 4 maker  maker_bopts.ctl maker_exe.ctl maker_opts.ctl



***Installation Documentation***

How to Install Standard MAKER

!!IMPORTANT NOTE FOR MAC OS X USERS!!
You will need to install developer tools (i.e. Xcode) from the App Store or
your installation disk. Also install fink (http://www.finkproject.org/) and
then install glib2-dev via fink (i.e. 'fink install glib2-dev').


**EASY INSTALL

1.  Go to the .../maker/src/ directory and run 'perl Build.PL' to configure.

2.  Accept default configuration options by just pressing enter. See MPI INSTALL
    in next section if you decide to configure for MPI.

3.  type './Build install' to complete the installation.

4.  If anything fails, either use the ./Build file commands to retry the failed
    section (i.e. './Build installdeps' and './Build installexes') or follow the
    detailed install instructions in the next section to manually install missing
    modules or programs. Use ./Build status to see available commands.

       ./Build status           #Shows a status menu
       ./Build installdeps      #installs missing PERL dependencies
       ./Build installexes      #installs missing external programs
       ./Build install          #installs MAKER

    Note: You do not need to be root.  Just say 'yes' to 'local install' when
    running './Build installdeps' and dependencies will be installed under
    .../maker/perl/lib, also missing external tools will be installed under
    .../maker/exe when running './Build installexes'.

    Note: For failed auto-download of external tools, when using the command
    './Build installexes', the .../maker/src/locations file is used to identify
    download URLs. You can edit this file to point to any alternate locations.



**MPI INSTALL

!!IMPORTANT!!
MAKER is not compatible with MVAPICH2. Use OpenMPI or MPICH. If using MPICH, make
sure to enable shared libaries during installation (this is not the default). If
using OpenMPI, make sure to set LD_PRELOAD to the location of libmpi.so before
even trying to install MAKER. It must also be set before running MAKER (or any
program that uses OpenMPI's shared libraries), so it's best just to add it to
your ~/.bash_profile. (i.e. export LD_PRELOAD=/usr/local/openmpi/lib/libmpi.so).


1.  Say yes to the 'configure for MPI' question when running 'perl Build.PL' in
    step 1 of the EASY INSTALL.

2.  Give path to 'mpicc'. Note to make sure you do not give the path to 'mpicc'
    from another MPI flavor that might be installed on your system.

3.  Give path to the folder containing 'mpi,h'. Note to make sure you do not
    give the path to a folder from another MPI flavor that might be installed
    on your system. Mixing MPI flavors for 'mpicc' and 'mpi.h' will cause
    failures. Make sure to read and confirm the auto-detected paths.

4.  Finish installation according to steps 2-4 of the EASY INSTALL

    Note: For OpenMPI you may also want to set OMPI_MCA_mpi_warn_on_fork=0 in
    your ~/.bash_profile to turn off certain nonfatal warnings.

    Note: If jobs hang or freeze when using mpiexec under OpenMPI try adding
    the '-mca btl ^openib' flag to mpiexec command when running MAKER.

        Example: mpiexec -mca btl ^openib -n 20 maker



**DETAILED INSTALL (for installing prerequisites manually)

1.  You need to have perl 5.8.0 or higher installed.  You can verify this by
    typing perl -v on the command line in a terminal.

    You will also need to install the following perl modules from CPAN.
      *DBI
      *DBD::SQLite
      *forks
      *forks::shared
      *File::Which
      *Perl::Unsafe::Signals
      *Bit::Vector
      *Inline::C
      *IO::All
      *IO::Prompt

    a. Type 'perl -MCPAN -e shell' to access the CPAN shell.  You may
       have to answer some configuration questions if this is your first time
       starting CPAN. You can normally just hit enter to accept CPAN defaults.
       You may have to be logged in as 'root' or use sudo to install modules
       via CPAN. If you don't have root access, then install local::lib from
       http://www.cpan.org using the bootstrap method to setup a non-root CPAN
       install location.

    b. Type 'install DBI' in CPAN to install the first module, then type
      'install DBD::SQLite' to install the next one, and so on.

    c. Alternatively you can download moadules from http://www.cpan.org/.
       Just follow the instructions that come with each module to install.

2.  Install BioPerl 1.6 or higher.  Download the Core Package from
    http://www.bioperl.org

  -quick and dirty installation-
  (with this option, not all of bioperl will work, but what MAKER uses will)

  a.  Download and unpack the most recent BioPerl package to a directory of your
      choice, or use Git to access the most current version of BioPerl. See
      http://www.bioperl.org for details on how to download using Git.
      You will then need to set PERL5LIB in your .bash_profile to the location
      of bioperl (i.e. export PERL5LIB="/usr/local/bioperl-live:$PERL5LIB").

  -full BioPerl instalation via CPAN-
  (you will need sudo privileges, root access, or CPAN configured for local
   installation to continue with this option)

  a.  Type perl -MCPAN -e shell into the command line to set up CPAN on your
      computer before installing bioperl (CPAN helps install perl dependencies
      needed to run bioperl).  For the most part just accept any default options
      by hitting enter during setup.
  b.  Type install Bundle::CPAN on the cpan command line.  Once again just press
      enter to accept default installation options.
  c.  Type install Module::Build on the cpan command line.  Once again just
      press enter to accept default installation options.
  d.  Type install Bundle::BioPerl on the cpan command line.  Once again press
      enter to accept default installation options.

  -full BioPerl instalation from download-
  a.  Unpack the downloaded bioperl tar file to the directory of your choice or
      use Git to get the most up to date version.  Then use the terminal
      to find the directory and type perl Build.PL in the command line, then
      type ./Build test, then if all tests pass, type ./Build install.  This
      will install BioPerl where perl is installed.  Press enter to accept all
      defaults.

3.  Install either WuBlast or NCBI-BLAST using instruction in 3a and 3b

3a.  Install WuBlast 2.0 or higher (Alternatively install NCBI-BLAST in 3b)
    (WuBlast has become AB-Blast and is no longer freely available, so if you
    are not lucky enough to have an existing copy of WuBlast, you can use NCBI
    BLAST or BLAST+ instead)

  a.  Unpack the tar file into the directory of your choice (i.e. /usr/local).
  b.  Add the following in your .bash_profile (value depend on where you chose
      to install wublast):
		export WUBLASTFILTER="/usr/local/wublast/filter"
		export WUBLASTMAT="/usr/local/wublast/matrix"
  c.  Add the location where you installed WuBlast to your PATH variable in 
      .bash_profile (i.e. PATH="/usr/local/wublast:$PATH").

3b.  Install NCBI-BLAST 2.2.X or higher (Alternatively install WuBlast in 3a)

  a.  Unpack the tar file into the directory of your choice (i.e. /usr/local).
  b.  Add the location where you installed NCBI-BLAST to your PATH variable in
      .bash_profile (i.e. PATH="/usr/local/ncbi-blast:$PATH").

4.  Install SNAP.  Download from http://korflab.ucdavis.edu/software.html

  a.  Unpack the SNAP tar file into the directory of your choice (ie /usr/local)
  b.  Add the following to your .bash_profile file (value depends on where you
      choose to install snap):  export ZOE="/usr/local/snap/Zoe"
  c.  Navigate to the directory where snap was unpacked (i.e. /usr/local/snap)
      and type make
  d.  Add the location where you installed SNAP to your PATH variable in
      .bash_profile (i.e. export PATH="/usr/local/snap:$PATH").


5.  Install RepeatMasker. Download from http://www.repeatmasker.org

  a.  The most current version of RepeatMasker requires a program called TRF.
      This can be downloaded from http://tandem.bu.edu/trf/trf.html
  b.  The TRF download will contain a single executable file.  You will need to
      rename the file from whatever it is to just 'trf' (all lower case).
  c.  Make TRF executable by typing chmod a+x trf.  You can then move this file
      wherever you want.  I just put it in the .../RepeatMasker directory.
  d.  Unpack RepeatMasker to the directory of your choice (i.e. /usr/local).
  e.  If you do not have WuBlast installed, you will need to install RMBlast.
      We do not recomend using cross_match, as RepeatMasker performance will suffer.
  f.  Now in the RepeatMasker directory type perl ./configure in the command
      line. You will be asked to identify the location of perl, rmblast/wublast,
      and trf.  The script expects the paths to the folders containing the
      executables (because you are pointing to a folder the path must end in a 
      '/' character or the configuration script throws a fit).
  g.  Add the location where you installed RepeatMasker to your PATH variable in
      .bash_profile (i.e. export PATH="/usr/local/RepeatMasker:$PATH").
  h.  You must register at http://www.girinst.org and download the Repbase
      repeat database, Repeat Masker edition, for RepeatMasker to work.
  i.  Unpack the contents of the RepBase tarball into the RepeatMasker/Libraries
      directory.


6.  Install Exonerate 2.2.  Download from http://www.ebi.ac.uk/~guy/exonerate

  a.  Exonerate has pre-comiled binaries for many systems; however, for Mac OS-X
      you will have to download the source code and complile it yourself by
      following steps b though d.
  b.  You need to have Glib 2.0 installed.  The easiest way to do this on a Mac
      is to install fink and then type 'fink install glib2-dev' in the terminal.
  c.  Change to the directory containing the Exonerate package to be compiled.
  d.  To install exonerate in the directory /usr/local/exonerate, type:
      ./configure -prefix=/usr/local/exonerate -> then type make -> then type
      make install
  e.  Add the location where you installed exonerate to your PATH variable in
      .bash_profile (i.e. export PATH="/usr/local/exonerate/bin:$PATH").


7.  Install MAKER.  Download from http://www.yandell-lab.org

  a.  Unpack the MAKER tar file into the directory of your choice (i.e.
      /usr/local).
  b.  Go to the MAKER src/ directory.
  c.  Configure using --> perl Build.PL
  D.  Install using --> ./Build install
  b.  Remember to add the following to your .bash_profile if you haven't already:
	export ZOE="where_snap_is/Zoe"
	export AUGUSTUS_CONFIG_PATH="where_augustus_is/config
  c.  Add the location where you installed MAKER to your PATH variable in
      .bash_profile (i.e. export PATH=/usr/local/maker/bin:$PATH).
  d.  You can now run a test of MAKER by following the instructions in the MAKER
      README file.



**(OPTIONAL COMPONENTS)

1. Augustus 2.0 or higher. Download from http://bioinf.uni-greifswald.de/augustus/
2. GeneMark-ES. Download from http://exon.biology.gatech.edu
3. FGENESH 2.4 or higher. Purchase from http://www.softberry.com
4. GeneMarkS. Download from http://exon.biology.gatech.edu

!!Read their installation documentation.




