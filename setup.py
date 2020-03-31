from setuptools import setup, find_packages

setup(
    name='Meta-Waffle',
    version='0.0.1',
    packages=find_packages(),
    include_package_data=True,
    package_dir={'test': 'test'},
    package_data={'test': ['*.pickle', 'data/*.bam*', 'data/*.tsv', 'data/*.bed', 'data/*.pickle']},
    scripts=['scripts/waffle-peaks.py',
             'scripts/waffle-plot.py',
             'scripts/waffle-predict.py',
             'scripts/waffle-bam2count.py'],
    license='GPLv3'
)
