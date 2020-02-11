#!/usr/bin/env bash
# Daniel B. Stribling
# Renne Lab, University of Florida
# Hybkit Project : http://www.github.com/RenneLab/hybkit
# Run the sphinx-apidoc with custom settings to generate api documentation.

sphinx-apidoc --force --separate --module-first \
              --maxdepth 6 -o source_ex ../ ../setup.py
#grep -v inheritance source/hybkit.rst > source/hybkit.rst.modified
#mv source/hybkit.rst.modified source/hybkit.rst

