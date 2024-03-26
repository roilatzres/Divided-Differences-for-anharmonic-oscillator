#!/bin/bash

# Compile the program
g++ main.cpp -o main

# Check if the compilation was successful
if [ $? -eq 0 ]; then
    echo "Compilation successful. Running the program:"
    ./main
else
    echo "Compilation failed."
fi