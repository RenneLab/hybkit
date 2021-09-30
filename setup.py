#!/usr/bin/env python3
# Daniel Stribling
# Renne Lab, University of Florida
# Hybkit Project : https://www.github.com/RenneLab/hybkit

'''
Setup module for the hybkit project.
'''

import setuptools
import os
import hybkit
import glob
import fnmatch

# Set project directory
proj_dir = os.path.abspath(os.path.dirname(__file__))

# Get the long description from the README file
with open(os.path.join(proj_dir, 'README.rst'), encoding='utf-8') as f:
    long_description = ''
    for line in f:
        if not 'include::' in line:
            long_description += line

# Get the remaining project details variables from the "__about__.py" file.
about_vars = {}
with open(os.path.join(proj_dir, 'hybkit', '__about__.py')) as f:
    exec(f.read(), about_vars)

# Dynamically generate reference data file tuples:
data_files = []
data_file_dirs = ['', 'scripts', 'reference_data', 'hybkit']
sample_directory_dirs = glob.glob('sample_0*')
data_file_dirs += sample_directory_dirs
for item in glob.glob('docs/**', recursive=True):
    if os.path.isdir(item) and not item.startswith(os.path.join('docs','_')):
        data_file_dirs.append(item)

ignore_file_patterns = []
with open('.gitignore', 'r') as git_ignore:
    for line in git_ignore:
        line = line.strip()
        if line.startswith('#') or not line:
            continue
        ignore_file_patterns.append(line)

for dir_name in data_file_dirs:
    file_list = [f for f in glob.glob(os.path.join(dir_name, '*'))
                 if not (
                         os.path.isdir(f)
                         or any(fnmatch.fnmatch(f, ignore) for ignore in ignore_file_patterns)
                        )]
    target_dir_name = os.path.join(about_vars['name_and_version'], dir_name)
    data_files.append((target_dir_name, file_list))

setuptools.setup(
    name='hybkit',
    version=about_vars['__version__'],
    description=about_vars['description'],
    long_description=long_description,
    long_description_content_type='text/x-rst',
    url=about_vars['project_url'],
    author=about_vars['__author__'],
    author_email=about_vars['__contact__'],
    classifiers=about_vars['classifiers'],
    keywords=about_vars['keywords'],
    packages=['hybkit'],
    package_dir={'hybkit': 'hybkit'},
    package_data={
                  'hybkit': [os.path.basename(f) for f in glob.glob('hybkit/*')
                             if not (os.path.basename(f).endswith('.py')
                             or os.path.basename(f).endswith('__'))]
                 },
    scripts=glob.glob('scripts/*'),
    python_requires='>=3.6',
    project_urls=about_vars['info_urls'],
    data_files=data_files,
    install_requires=[
                      'matplotlib',
                      'importlib_resources',
                     ],
 
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
