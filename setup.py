import os
from setuptools import setup, find_packages
def calculate_version(inputfile):
    version_list =  [x.split('\'')[1] for x in open(inputfile, 'r')
                     if x.startswith('__version__')]
    if version_list:
        return version_list[0]
    else:
        return '1.0'

_version_f = os.path.join(os.path.split(os.path.abspath(__file__))[0],
                          'pmid_utils', '_version.py'
)
package_version = calculate_version(_version_f)

setup(
    name='pmid_utils',
    version=package_version,
    py_modules=['pmid_utils'],
    packages=find_packages(),
    install_requires=[
        'numpy',
        'pandas',
    ],

    description='PMID Utils',
    license='GPL-3.0',
    classifiers=[
        'Intended Audience :: Bioinformatics',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
    ],
    url='https://github.com/asalt/pmid_utils',
    author='Alexander Saltzman',
    author_email='saltzman@bcm.edu,'


)
