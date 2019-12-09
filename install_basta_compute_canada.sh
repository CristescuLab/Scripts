#!/usr/bin/env bash
module purge
module load nixpkgs/16.09  gcc/7.3.0 leveldb/1.20 scipy-stack/2018b python/2.7.14
pip install Plyvel --user
pip install git+https://github.com/timkahlke/BASTA.git --user
basta taxonomy
rm $HOME/.basta/taxonomy/*.gz*
basta download gb
echo "BEFORE YOU RUN BASTA, EXECUTE:"
echo "module purge; module load nixpkgs/16.09  gcc/7.3.0 leveldb/1.20 scipy-stack/2018b python/2.7.14"