#!/usr/bin/env bash

## USAGE: bash sync_software.sh
# You need to be IN the folder you want to sync

h=`hostname`
if [[ $h = *"cedar"* ]]; then
    current_address="cedar.computecanada.ca"
    other_address="graham.computecanada.ca"
    else
    other_address="cedar.computecanada.ca"
    current_address="graham.computecanada.ca"
fi

if [[ $PWD = *"$1"* ]]; then
   rsync -rltv ${current_address}:$PWD/* ./
   rsync -rltv ./* ${other_address}:$PWD/; else
   echo "YOU ARE NOT IN THE RIGHT FOLDER!!!"
fi