# Here are release instructions

# Version bump the value. Remove '-beta'
vim lib/Version.pm

# Tag the branch and push the tag. You'll need to have a gpg signature for this. Extra precaution
git tag -s 1.1 -m 'Release 1.1'
git push --tags

# Version bump the value. Until a realease add '-beta'. Like 1.2-beta
vim lib/Version.pm

# Create a zipped bundle and load it to bitbucket?
