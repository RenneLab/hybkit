#!/usr/bin/env python3
# Daniel Stribling
# Renne Lab, University of Florida
# Hybkit Project : https://www.github.com/RenneLab/hybkit

'''
Setup module for the hybkit project.
'''

import setuptools
import os

if __name__ == '__main__':
    
    # Set project directory
    proj_dir = os.path.abspath(os.path.dirname(__file__))
    
    # Get the long description from the README file
    with open(os.path.join(proj_dir, 'README.md'), encoding='utf-8') as f:
        long_description = f.read()
    
    # Get the remaining project details variables from the "__about__.py" file.
    about_vars = {}
    with open(os.path.join(proj_dir, 'hybkit', '__about__.py')) as f:
        exec(f.read(), about_vars)
    
    setuptools.setup(
        name='hybkit',
        version=about_vars['__version__'],
        description=about_vars['description'],
        long_description=long_description,
        long_description_content_type='text/markdown',
        url=about_vars['project_url'],
        author=about_vars['__author__'],
        author_email=about_vars['__contact__'],
        classifiers=about_vars['classifiers'],
        keywords='genetics genomics ribonomics bioinformatics CLASH qCLASH miRNA',
        packages=setuptools.find_packages(where='src'),
        python_requires='>=3.6',
        project_urls=about_vars['info_urls'],
    
        #install_requires=['biopython'],  # Optional
    
        # To provide executable scripts, use entry points in preference to the
        # "scripts" keyword. Entry points provide cross-platform support and allow
        # `pip` to create the appropriate form of executable for the target
        # platform.
        #
        # For example, the following would provide a command called `sample` which
        # executes the function `main` from this package when invoked:
        #entry_points={  # Optional
        #    'console_scripts': [
        #        'sample=sample:main',
        #    ],
        #},
    )
