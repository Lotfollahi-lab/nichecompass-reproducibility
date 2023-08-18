#!/bin/bash

dataset="$1"

cd ../../../datasets/benchmark_data/deeplinc
cd "$dataset"
mkdir runs
cd runs

# Nested loop to execute run1 to run10
for i in {1..10}
do
    mkdir "run$i"
    cd "run$i"

    # Set the adj file based on the current run number
    if [[ "$i" -le 2 ]]; then
        adj_file="adj4.csv"
    elif [[ "$i" -le 4 ]]; then
        adj_file="adj8.csv"
    elif [[ "$i" -le 6 ]]; then
        adj_file="adj12.csv"
    elif [[ "$i" -le 8 ]]; then
        adj_file="adj16.csv"
    else
        adj_file="adj20.csv"
    fi

    # Build and execute the Python script command with line breaks
    cmd="python ../../../../../../scripts/single_sample_method_benchmarking/deeplinc/deeplinc.py \
        -e ../../counts.csv \
        -a ../../$adj_file \
        -c ../../coords.csv \
        -r ../../cell_types.csv \
        -n run$i \
        -i 40 \
        --seed $((i-1))"
    eval "$cmd"

    cd ..
done

cd ../..

