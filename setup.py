#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages

with open('README.md') as readme_file:
    readme = readme_file.read()

with open('CHANGELOG.md') as history_file:
    history = history_file.read()

#with open('requirements.txt', 'r') as f:
#    requirements = f.read().splitlines()

setup_requirements = ['setuptools_scm']

test_requirements = []

setup(
    author="LEA - Uni Paderborn",
    author_email='upblea@mail.upb.de',
    python_requires='>=3.7',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Natural Language :: English',
        'Topic :: Scientific/Engineering',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Operating System :: Microsoft :: Windows',
        'Operating System :: POSIX :: Linux',
        'Environment :: MacOS X'
    ],
    description="FEM Magnetics Toolbox",
    install_requires=['pandas>=1.2.1',
					  'numpy>=1.20.3',
					  'matplotlib>=3.3.4',
					  'gmsh>=4.8.4',
					  'onelab>=1.0',
					  'scipy>=1.7.2',
                      'pytest',
					  ],
    license="GNU General Public License v3",
    long_description=readme + '\n\n' + history,
    long_description_content_type="text/markdown",
    include_package_data=True,
    keywords='femmt',
    name='femmt',
    packages=find_packages(include=['femmt', 'femmt.*']),
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    extras_require={},
    url='https://github.com/upb-lea/FEM_Magnetics_Toolbox',
    project_urls={
        "Documentation": "https://github.com/upb-lea/FEM_Magnetics_Toolbox",
        "Source Code": "https://github.com/upb-lea/FEM_Magnetics_Toolbox",
    },
    version='0.1.2',
    zip_safe=False,
    data_files=[('', ['CHANGELOG.md'])]
)
