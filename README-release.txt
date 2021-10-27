# Here are the release instructions

# Version bump the value. Remove '-beta'
vim VERSION

# Update mugqic_pipelines module version
vim resources/modules/genpipes.sh
(VERSION=3.5.0)

# Recreate the pipelines/<pipeline>/README.md using --help with markdown output
for script in \
  pipelines/ampliconseq/ampliconseq.py \
  pipelines/chipseq/chipseq.py \
  pipelines/dnaseq/dnaseq.py \
  pipelines/dnaseq_high_coverage/dnaseq_high_coverage.py \
  pipelines/methylseq/methylseq.py \
  pipelines/pacbio_assembly/pacbio_assembly.py \
  pipelines/rnaseq/rnaseq.py \
  pipelines/rnaseq_denovo_assembly/rnaseq_denovo_assembly.py \
  pipelines/rnaseq_light/rnaseq_light.py \
  pipelines/tumor_pair/tumor_pair.py
do
  echo $script
  pipeline=$(cut -d '/' -f2 <<< "$script")
  diagram_folder_path=https://bitbucket.org/mugqic/genpipes/raw/master/resources/workflows/
  $script --help > `dirname $script`/README.md
  if grep -PA1 "^----$" `dirname $script`/README.md | grep -P "\w+:"
  then
    perl -pi -e 'undef $/; s#----\n(\w*):\n#----\n```\n!['$pipeline' $1 workflow diagram]('$diagram_folder_path'GenPipes_'$pipeline'_$1.resized.png)\n[download full-size diagram]('$diagram_folder_path'GenPipes_'$pipeline'_$1.png)\n```\n$1:\n#g' `dirname "$script"`/README.md
  else
    perl -pi -e 'undef $/; s#Steps:\n------\n#Steps:\n```\n!['$pipeline' workflow diagram]('$diagram_folder_path'GenPipes_'$pipeline'.resized.png)\n[download full-size diagram]('$diagram_folder_path'GenPipes_'$pipeline'.png)\n```\n------\n#g' `dirname "$script"`/README.md
  fi
done

# Special case for HiCSeq pipeline which needs an extra '-e' parameter to output the README properly
pipelines/hicseq/hicseq.py -e DpnII --help > pipelines/hicseq/README.md
perl -pi -e 'undef $/; s#----\n(\w*):\n#----\n```\n![hicseq $1 workflow diagram]('$diagram_folder_path'GenPipes_hicseq_$1.resized.png)\n[download full-size diagram]('$diagram_folder_path'GenPipes_hicseq_$1.png)\n```\n$1:\n#g' pipelines/hicseq/README.md

# For Tumor Pair pipeline, produce a regular README until Rob's branch is merge to master : will then contain 3 tumor pair protocols for which we alraedy have workflow diagrams ready.
pipelines/tumor_pair/tumor_pair.py --help > pipelines/tumor_pair/README.md

# For Illumina pipeline, just produce a regular README, without diagram
pipelines/illumina_run_processing/illumina_run_processing.py --help > pipelines/illumina_run_processing/README.md

# Tag the branch and push the tag. You'll need to have a gpg signature for this. Extra precaution
git tag -s 3.5.0 -m 'Release 3.5.0'
git push -u origin --tags

# Recreate the CHANGELOG.md
bash ~/work/repo/dump_ChangeLog.sh > CHANGELOG.md
git commit -a -m "Version bump to 3.5.0"

# Create a release tarball archive
git archive --format=tar --prefix=genpipes-3.5.0/ 3.5.0 | gzip > ~/genpipes-3.5.0.tar.gz

# Upload this archive in
https://bitbucket.org/mugqic/genpipes/downloads

# Version bump the value. Until the next release, add '-beta' e.g. 3.5.1-beta
vim VERSION
git commit -m "Version bump to 3.5.1-beta" VERSION
git push

# Deploy GenPipes-<VERSION> as a module on all clusters

# Send a message to the mailing list:
mugqic_pipelines@googlegroups.com

# In JIRA, add a release date to the 'Version' category of the administer project
