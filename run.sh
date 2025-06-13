#!/bin/bash
echo "Compiling the program..."

# Compile the program
g++ -o main_orig main_orig.cpp  

echo "Compilation finished."

# Check if the compilation was successful
if [ $? -eq 0 ]; then
    echo "Compilation successful. Running the program:"
    
    # Array of ampof values
    ampof_values=(0.08 0.1 0.2)

    # Array of q values
    q_values=(10 11 12 13 14)

    # Run main with each combination of ampof and q values
    for ampof in "${ampof_values[@]}"; do
        for q in "${q_values[@]}"; do
            echo "Running with ampof=$ampof and q=$q"
            ./main_orig.exe $ampof $q > "amp_${ampof}_q${q}.txt"
        done
    done
else
    echo "Compilation failed."
fi