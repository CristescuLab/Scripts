#!/usr/bin/env bash

## USAGE: bash sync_software.sh <path2folder> <source_server> <target_server>
# You need to be IN the folder you want to sync
path2folder=$1
source_server=$2
target_server=$3

echo "Syncing ${source_server}:${path2folder} to ${target_server}:${path2folder}"
rsync -rltv ${path2folder}/* ${target_server}:${path2folder}/
rsync -rltv ${target_server}:${path2folder}/* ${path2folder}/
