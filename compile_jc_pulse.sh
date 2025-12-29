!/bin/bash
echo "Compiling the program..."

# Compile the program

g++ main_jc_multi_piece.cpp include/cnpy/cnpy.cpp -Iinclude/cnpy -Idlib/dlib -I"%CONDA_PREFIX%\Library\include" -L"%CONDA_PREFIX%\Library\lib" -lz -std=c++17 -O2 -o main_jc_multi_piece.exe


g++ main_const_pulse_multi_piece_remake.cpp -std=c++17 -O2 -o main_const_pulse_multi_piece_remake_eem
# g++ main_jc_multi_piece.cpp -Idlib/dlib -std=c++17 -O2 -o main_jc_multi_piece
# g++ main_const_pulse_multi_piece_test.cpp -Idlib/dlib -std=c++17 -O2 -o main_const_pulse_multi_piece_test_eem
# g++ main_sin_pulse_multi_piece.cpp -Idlib/dlib -std=c++17 -O2 -o main_sin_pulse_multi_piece_eem

echo "Compilation finished."


# # # Check if the compilation was successful
# if [ $? -eq 0 ]; then
#     echo "Compilation successful. Running the program:"
    
#     # Array of ampof values
#     ampof_values=(0.035 0.04 0.045)

#     # Array of q values
#     q_values=(16 17 18 19 20)

#     # target values
#     target_values=(16 17 18 19 20)  

#     # Run main with each combination of ampof and q values
#     for ampof in "${ampof_values[@]}"; do
#         for q in "${q_values[@]}"; do
#             for max_target in "${target_values[@]}"; do
#                 echo "Running with ampof=$ampof, q=$q, and max_target=$max_target"
#                 ./main_const_pulse_multi_piece_eem0.exe $ampof $q $max_target  32 
#                 # > "test_multi_piece_amp_${ampof}_q${q}_n${max_target}_p64.txt"

#             done
#         done
#     done
# else
#     echo "Compilation failed."
# fi#!/bin/bash
# # echo "Compiling the program..."

# # Compile the program

# # g++ main_const_pulse_multi_piece.cpp -Idlib/dlib -std=c++17 -O2 -o main_const_pulse_multi_piece

# # echo "Compilation finished."


# # # Check if the compilation was successful
# if [ $? -eq 0 ]; then
#     echo "Compilation successful. Running the program:"
    
#     # Array of ampof values
#     ampof_values=(0.035)

#     # Array of q values
#     q_values=(12 13 14 15)

#     # target values
#     target_values=(19 20)  

#     # Run main with each combination of ampof and q values
#     for ampof in "${ampof_values[@]}"; do
#         for q in "${q_values[@]}"; do
#             for max_target in "${target_values[@]}"; do
#                 echo "Running with ampof=$ampof, q=$q, and max_target=$max_target"
#                 ./main_const_pulse_multi_piece2.exe $ampof $q $max_target  32 > "test_multi_piece_amp_${ampof}_q${q}_n${max_target}_p64.txt"
#             done
#         done
#     done
# else
#     echo "Compilation failed."
# fi

# echo "Press ENTER to close..."
# read -r