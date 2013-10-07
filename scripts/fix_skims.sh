#! /bin/bash

shopt -s nullglob

for file in ../data/*.root
do
  noext=${file%.root}
  ./scripts/fix_skimmed_file.exe "$file" &> "logs/fix_skimmed_fFile_${noext##*/}.log" &
done

wait
echo "Done"
exit 0;