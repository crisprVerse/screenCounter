#!/bin/bash

set -e
set -u

# Transplant the kaori code.
if [ ! -e source-kaori ]
then
    git clone https://github.com/crisprVerse/kaori source-kaori
else
    (cd source-kaori && git checkout master && git pull)
fi
(cd source-kaori && git checkout v1.1.1)

rm -rf kaori
mkdir -p kaori
cp -r source-kaori/include/kaori/* kaori

# Transplant the byteme code.
if [ ! -e source-byteme ]
then
    git clone https://github.com/LTLA/byteme source-byteme
else
    (cd source-byteme && git checkout master && git pull)
fi
(cd source-byteme && git checkout v1.0.1)

rm -rf byteme
mkdir -p byteme
cp source-byteme/include/byteme/* byteme
