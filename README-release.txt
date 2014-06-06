# Here are release instructions

# Version bump the value. Remove '-beta'
vim lib/Version.pm
git commit -m "Version bump to 1.3"

# Tag the branch and push the tag. You'll need to have a gpg signature for this. Extra precaution
git tag -s 1.3 -m 'Release 1.3'
git push --tags

# Version bump the value. Until a realease add '-beta'. Like 1.4-beta
vim lib/Version.pm
git commit -m "Version bump to 1.4-beta"
git push

# Add a blogpost here:
#  https://biowiki.atlassian.net/wiki/pages/viewrecentblogposts.action?key=PS
#  Press the CREATE button and choose "Blog Post"

# In JIRA Add a release date to the 'Version' category of the administer project
