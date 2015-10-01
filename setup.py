#!/usr/bin/env python
# -*- coding: utf-8 -*-


try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup


with open('README.md') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read().replace('.. :changelog:', '')

requirements = [
    'scikit-bio>=0.2.2,<1.0',
    'tempdir',
]

test_requirements = [
    'nose'
]

exec(open('tree2tax/version.py').read()) # loads __version__

setup(
    name='tree2tax',
    version=__version__,
    description="Automatic taxonomy through consistent application of tree-based thresholding",
    long_description=readme + '\n\n' + history,
    author="Ben Woodcroft",
    author_email='b.woodcroft near uq.edu.au',
    url='https://github.com/wwood/tree2tax',
    packages=[
        'tree2tax',
    ],
    package_dir={'tree2tax':
                 'tree2tax'},
    include_package_data=True,
    install_requires=requirements,
    license="LGPL3",
    zip_safe=False,
    keywords='tree2tax',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: BSD License',
        'Natural Language :: English',
        "Programming Language :: Python :: 2",
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
    ],
    test_suite='tests',
    tests_require=test_requirements,
    scripts=[
        'bin/autotaxonomy',
        'bin/percent_identity_vs_tree_distance',
        'bin/threshold_estimator',
        'bin/tree2tax'
    ]
)
