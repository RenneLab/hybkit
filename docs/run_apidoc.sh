#!/usr/bin/env bash
# Daniel B. Stribling
# Renne Lab, University of Florida
# Hybkit Project : http://www.github.com/RenneLab/hybkit
# Run the sphinx-apidoc with custom settings to generate api documentation.

sphinx-apidoc --force --maxdepth 6 -o source ../ ../setup.py

