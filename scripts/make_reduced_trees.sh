#! /bin/bash

shopt -s nullglob

for file in ../data/*.root
do
  noext=${file%.root}
  ./scripts/make_reduced_tree.exe -i "$file" &> "logs/make_reduced_tree_${noext##*/}.log" &
done

wait
echo "Done"
exit 0;
