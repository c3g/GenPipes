# Here are the release instructions

# Version bump the value. Remove '-beta'
vim VERSION

# Tag the branch and push the tag. You'll need to have a gpg signature for this. Extra precaution
git tag -s mugqic_pipeline-2.0.0 -m 'mugqic_pipeline release 2.0.0'
git push --tags

# Version bump the value. Until the next release, add '-beta' e.g. 2.1.0-beta
vim VERSION

# Add a blogpost here:
#  https://biowiki.atlassian.net/wiki/pages/viewrecentblogposts.action?key=PS
#  Press the CREATE button and choose "Blog Post"

# In JIRA, add a release date to the 'Version' category of the administer project
