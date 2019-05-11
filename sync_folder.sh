#!/usr/bin/env bash

## USAGE: bash sync_software.sh <path2folder> <target_server>
# You need to be IN the folder you want to sync
path2folder=$1
target_server=$2

echo "Syncing ${path2folder} to ${target_server}:${path2folder}"
rsync -rltv --ignore-existing --human-readable ${path2folder}/* ${target_server}:${path2folder}/
rsync -rltv --ignore-existing --human-readable ${target_server}:${path2folder}/* ${path2folder}/
