#!/bin/bash

set -e
set -u

# Transplant the kaori code.
if [ -e source-kaori ]
then
    (cd source-kaori && git pull)
else
    git clone https://github.com/LTLA/kaori source-kaori
fi
rm -rf kaori
mkdir -p kaori
cp source-kaori/include/kaori/* kaori

# Transplant the byteme code.
if [ -e source-byteme ]
then
    (cd source-byteme && git pull)
else
    git clone https://github.com/LTLA/byteme source-byteme
fi
rm -rf byteme
mkdir -p byteme
cp source-byteme/include/byteme/* byteme
