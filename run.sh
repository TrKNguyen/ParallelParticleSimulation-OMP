#!/bin/bash

if [ ! -f "./sim" ]; then
    echo "Executable 'sim' not found. Please compile your project first."
    exit 1
fi

for n in 10 20 30 40 50 60 70; do
    input_file="tests/small/random${n}.in"
    if [ ! -f "$input_file" ]; then
        echo "Input file '$input_file' not found. Skipping..."
        continue
    fi
    echo "Running simulation with $input_file using 8 threads..."
    ./sim "$input_file" 8
    echo "-------------------------------------------"
done
