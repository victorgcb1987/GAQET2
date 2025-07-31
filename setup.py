from setuptools import setup, find_packages

setup(
        name='gaqet',
        version='1.11.15',
        packages=find_packages(),
        include_package_data=True,
        install_requires=[],
        entry_points={
                     'console_scripts': [
                                         'GAQET=GAQET.gaqet:main',
                                         'GAQET_PLOT=GAQET.gaqet_plot:main'
                                                                ],},
        author='Victor Garcia-Carpintero Burgos',
        description='Genome Annotation Quality Evaluation Tool (GAQET)',
        license='MIT',
        url='https://github.com/victorgcb1987/GAQET2',)