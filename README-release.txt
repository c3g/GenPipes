# Here are the release instructions

# Version bump the value. Remove '-beta'
vim VERSION

# Update mugqic_pipelines module version
vim resources/modules/mugqic_pipelines.sh
(VERSION=3.1.0)

# Recreate the pipelines/<pipeline>/README.md using --help with markdown output
for f in \
  pipelines/ampliconseq/ampliconseq.py \
  pipelines/chipseq/chipseq.py \
  pipelines/dnaseq/dnaseq.py \
  pipelines/dnaseq_high_coverage/dnaseq_high_coverage.py \
  pipelines/illumina_run_processing/illumina_run_processing.py \
  pipelines/methylseq/methylseq.py \
  pipelines/pacbio_assembly/pacbio_assembly.py \
  pipelines/rnaseq/rnaseq.py \
  pipelines/rnaseq_denovo_assembly/rnaseq_denovo_assembly.py \
  pipelines/rnaseq_light/rnaseq_light.py \
  pipelines/tumor_pair/tumor_pair.py \
; do echo $f; $f --help > `dirname $f`/README.md; done
pipelines/hicseq/hicseq.py -e DpnII --help > pipelines/hicseq/README.md

# Tag the branch and push the tag. You'll need to have a gpg signature for this. Extra precaution
git tag -s 3.1.0 -m 'Release 3.1.0'
git push -u origin --tags

# Recreate the CHANGELOG.md
bash ~/work/repo/dump_ChangeLog.sh > CHANGELOG.md
git commit -a -m "Version bump to 3.1.0"

# Create a release tarball archive
git archive --format=tar --prefix=genpipes-3.1.0/ 3.1.0 | gzip > genpipes-3.1.0.tar.gz

# Upload this archive in
https://bitbucket.org/mugqic/mugqic_pipelines/downloads

# Version bump the value. Until the next release, add '-beta' e.g. 3.1.1-beta
vim VERSION
git commit -m "Version bump to 3.1.1-beta" VERSION
git push

# Deploy GenAP_Pipes-<VERSION> as a module on all clusters

# Send a message to the mailing list:
mugqic_pipelines@googlegroups.com

# In JIRA, add a release date to the 'Version' category of the administer project
