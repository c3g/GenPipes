# Here are the release instructions

# Version bump the value. Remove '-beta'
vim VERSION

# Update mugqic_pipelines module version
vim resources/modules/mugqic_pipelines.sh
(VERSION=2.0.0)

# Recreate the pipelines/<pipeline>/README.md using --help with markdown output
for f in \
  pipelines/dnaseq/dnaseq.py \
  pipelines/rnaseq/rnaseq.py \
  pipelines/rnaseq_denovo_assembly/rnaseq_denovo_assembly.py \
  pipelines/pacbio_assembly/pacbio_assembly.py \
  pipelines/chipseq/chipseq.py \
  pipelines/illumina_run_processing/illumina_run_processing.py \
; do echo $f; $f --help > `dirname $f`/README.md; done

git commit -a -m "Version bump to 2.0.0"

# Tag the branch and push the tag. You'll need to have a gpg signature for this. Extra precaution
git tag -s 2.0.0 -m 'Release 2.0.0'
git push -u origin --tags

# Create a release tarball archive
git archive --format=tar --prefix=mugqic_pipelines-2.0.0/ 2.0.0 | gzip > mugqic_pipelines-2.0.0.tar.gz

# Upload this archive in
https://bitbucket.org/mugqic/mugqic_pipelines/downloads

# Version bump the value. Until the next release, add '-beta' e.g. 2.1.0-beta
vim VERSION
git commit -m "Version bump to 2.1.0-beta" VERSION
git push

# Deploy mugqic_pipelines-<VERSION> as a module on all clusters

# Send a message to the mailing list:
mugqic_pipelines@googlegroups.com

# In JIRA, add a release date to the 'Version' category of the administer project
