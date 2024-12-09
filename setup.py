import setuptools

setuptools.setup(
    name='c3g-genpipes',
    version='5.0.3',
    packages=setuptools.find_packages(),
    url='https://bitbucket.org/mugqic/genpipes/',
    license='GNU Lesser General Public License',
    author='Canadian Center for Computational Genomics',
    author_email='pipelines@computationalgenomics.ca',
    maintainer_email='paul.stretenowich@mcgill.ca, mareike.janiak@computationalgenomics.ca, po.quirion@mcgill.com, edouard.henrion@computationalgenomics.ca',
    description='Several bioinformatics pipelines developed at McGill University Genome Centre',
    entry_points = {
    "console_scripts": [
        "genpipes = genpipes.__main__:main"] },
    scripts = [
        'genpipes/utils/cancel_running_pipeline.py',
        'genpipes/utils/ignstats.py',
        'genpipes/utils/__init__.py',
        'genpipes/utils/job2json_project_tracking.py',
        'genpipes/utils/job2json.py',
        'genpipes/utils/log_report.py',
        'genpipes/utils/mugqicValidator.py',
        'genpipes/utils/nanuq2mugqic_pipelines.py',
        'genpipes/utils/csvToreadset.R',
        'genpipes/utils/utils.py',
        'genpipes/utils/watch_portal_folder.py',
        'genpipes/utils/log_report.pl',
        'genpipes/utils/chunk_genpipes.sh',
        'genpipes/utils/dump_ChangeLog.sh',
        'genpipes/utils/submit_genpipes',
        'genpipes/utils/statistics/logToJson.py'],
    install_requires = [ "packaging >=   20.9" ],
    python_requires = '>=3.11',
    include_package_data=True
)

# TODO do that in CI/CD... when it is ripe...
# generate bash completeion
# shtab --shell={bash,tcsh,zsh} -u genpipes.__init__.make_parser
# build:
# python3 -m build --sdist
# upload to test pipy
# twine upload --repository testpypi dist/*

