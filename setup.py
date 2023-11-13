import setuptools

setuptools.setup(
    name='genpipes-suite',
    version='0.0_alpha',
    packages=setuptools.find_packages(),
    url='https://bitbucket.org/mugqic/genpipes/',
    license='GNU Lesser General Public License',
    author='Canadian Center for Computational Genomics',
    author_email='pipelines@computationalgenomics.ca',
    maintainer_email='po.quirion@mcgill.com, edouard.henrion@mcgill.ca',
    description='Several bioinformatics pipelines developed at McGill University Genome Centre',
    entry_points = {
    "console_scripts": [
        "genpipes = genpipes.__main__:main"] },
    install_requires = ["packaging >=   20.9"],
    include_package_data=True
)

# TODO do that in CI/CD... when it is ripe...
# generate bash completeion
# shtab --shell={bash,tcsh,zsh} -u genpipes.__init__.make_parser
# build:
# python3 -m build --sdist
# upload to test pipy
# twine upload --repository testpypi dist/*

