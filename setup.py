#!/usr/bin/env python3
# Daniel Stribling  |  ORCID: 0000-0002-0649-9506
# Renne Lab, University of Florida
# Hybkit Project : https://www.github.com/RenneLab/hybkit

"""Setup module for the hybkit project."""

import fnmatch
import glob
import os
import pprint

import setuptools

# import hybkit

# Set project directory
proj_dir = os.path.abspath(os.path.dirname(__file__))

# Get the long description from the README file
with open(os.path.join(proj_dir, 'README.rst'), encoding='utf-8') as f:
    long_description = ''
    for line in f:
        if 'include::' not in line:
            long_description += line

# Get the remaining project details variables from the "__about__.py" file.
about_vars = {}
with open(os.path.join(proj_dir, 'hybkit', '__about__.py')) as f:
    exec(f.read(), about_vars)  #noqa: S102

# Dynamically generate reference data file tuples:
data_files = []
data_file_dirs = ['', 'scripts', 'ref_data', 'hybkit']
example_directory_dirs = glob.glob('example_0*')
data_file_dirs += example_directory_dirs
for item in glob.glob('docs/**', recursive=True):
    if os.path.isdir(item) and not item.startswith(os.path.join('docs', '_')):
        data_file_dirs.append(item)

ignore_file_patterns = []
with open('.gitignore') as git_ignore:
    for line in git_ignore:
        use_line = line.strip()
        if use_line.startswith('#') or not use_line:
            continue
        ignore_file_patterns.append(use_line)

for dir_name in data_file_dirs:
    file_list = [f for f in glob.glob(os.path.join(dir_name, '*'))
                 if not (os.path.isdir(f)
                         or any(fnmatch.fnmatch(f, ignore) for ignore in ignore_file_patterns)
                         )]
    target_dir_name = os.path.join(about_vars['name_and_version'], dir_name)
    if dir_name == '':
        file_list += ['.gitignore']
    data_files.append((target_dir_name, file_list))

print('\nData Files:')
pprint.pp(data_files)

# Get Package Data:
package_data = [
    os.path.basename(f) for f in glob.glob('hybkit/*')
    if not (os.path.basename(f).endswith('.py')
            or os.path.basename(f).endswith('__'))
]
print('\nPackage Data:')
pprint.pp(package_data)

# Get Scripts
scripts = [s for s in glob.glob('scripts/*') if not s.endswith('__')]
print('\nScripts:')
pprint.pp(scripts)

# Print About Vars
print('\nAbout Variables:')
pprint.pp(about_vars)

setuptools.setup(
    name=about_vars['project_name'],
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
    package_data={'hybkit': package_data},
    scripts=scripts,
    python_requires=about_vars['python_requires'],
    project_urls=about_vars['info_urls'],
    data_files=data_files,
    install_requires=[
        'matplotlib',
        'biopython',
    ],
)
