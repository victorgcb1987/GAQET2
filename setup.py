from setuptools import setup, find_packages

setup(
        name='gaqet',
        version='1.11.13',
        packages=find_packages(),
        install_requires=[],
        entry_points={
                     'console_scripts': [
                                         'GAQET=GAQET2.GAQET:main',
                                         'GAQET_PLOT=GAQET2.GAQET_PLOT:main'
                                                                ],},
        author='Victor Garcia-Carpintero Burgos',
        description='Genome Annotation Quality Evaluation Tool (GAQET)',
        license='MIT',
        url='https://github.com/victorgcb1987/GAQET2',)