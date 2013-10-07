#! /bin/bash

shopt -s nullglob

for file in ../data/*.root
do
  noext=${file%SyncedSkim.root}
  hadd "${noext}SyncSkim.root" "$file"
done

wait
echo "Done"
exit 0;