#!/usr/bin/env bash
# Daniel Stribling  |  ORCID: 0000-0002-0649-9506
# Renne Lab, University of Florida
# Hybkit Project : https://www.github.com/RenneLab/hybkit

# Ensure code conforms to PEP8 style guide, with exceptions:
#    Maximum Line length = 100
#    E402 - Allow module imports after top of file. (Necessary due to code in __init__)
#    W503 - Put binary operators at beginning of continuation lines.
# 
# Note: Will raaise E127 errors due to lack of pretty wrapping of "with" statement continuation.

echo -e "\n\n\n\n\n\n"

if [ "$1" ]; then
  pycodestyle $1 --show-pep8 --max-line-length=99 --ignore=E402,W503
else
  for fn in $(ls *.py); do
    pycodestyle $fn --show-pep8 --max-line-length=99 --ignore=E402,W503
  done
fi

# Check code docstrings in Google style, with exceptions:
#    D107 - Allow __init__ not to have docstrings 
#    D212 - Allow(/require) mutli-line docstrings to start on second line

echo -e "\n\n\n\n\n\n"

pydocstyle -e --convention=google \
           --add-select=D213,D408,D409 \
           --add-ignore=D107,D212



