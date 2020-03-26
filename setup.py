from setuptools import setup, find_packages

setup(
    name='Meta-Waffle',
    version='0.0.1',
    packages=find_packages(),
    scripts=['scripts/waffle-peaks.py',
             'scripts/waffle-plot.py',
             'scripts/waffle-predict.py',
             'scripts/waffle-bam2count.py'],
    license='GPLv3'
)
