# Here are release instructions

# Version bump the value. Remove '-beta'
vim lib/Version.pm

# Tag the branch and push the tag. You'll need to have a gpg signature for this. Extra precaution
git tag -s 1.1 -m 'Release 1.1'
git push --tags

# Version bump the value. Until a realease add '-beta'. Like 1.2-beta
vim lib/Version.pm

# Add a blogpost here:
#  https://biowiki.atlassian.net/wiki/pages/viewrecentblogposts.action?key=PS
#  Press the CREATE button and choose "Blog Post"

# In JIRA Add a release date to the 'Version' category of the administer project
