from setuptools import setup, find_packages

setup(
        name='gaqet',
        version='1.10.13',
        packages=find_packages(),
        install_requires=[],
        entry_points={
                     'console_scripts': [
                                         'GAQET=GAQET.__main__:main',
                                                                ],},
        author='Victor Garcia-Carpintero Burgos',
        description='Genome Annotation Quality Evaluation Tool (GAQET2)',
        license='MIT',
        url='https://github.com/victorgcb1987/GAQET2',)