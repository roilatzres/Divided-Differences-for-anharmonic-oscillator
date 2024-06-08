#!/bin/bash

# Compile the program
g++ main.cpp -o main

# Check if the compilation was successful
if [ $? -eq 0 ]; then
    echo "Compilation successful. Running the program:"
    
    # Array of ampof values
    ampof_values=(0.015 )

    # Array of q values
    q_values=(10)

    # Run main with each combination of ampof and q values
    for ampof in "${ampof_values[@]}"; do
        for q in "${q_values[@]}"; do
            echo "Running with ampof=$ampof and q=$q"
            ./main $ampof $q > "amp_${ampof}_q_${q}.txt"
        done
    done
else
    echo "Compilation failed."
fi