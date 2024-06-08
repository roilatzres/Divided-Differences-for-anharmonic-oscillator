#include <iostream>
#include <vector>
#include <tuple>
#include <complex>
#include <fstream>
#include <functional>
#include <sstream>
#include <cmath>
#include "permutation.h"

// Function to save permutations and coefficients to a file
void save_results(const std::vector<std::vector<int>>& permutations, const std::vector<std::tuple<std::vector<std::vector<DivdiffElement>>, std::complex<double>, int>>& coefficients, const std::string& filename) {
    std::ofstream file(filename, std::ios::out);
    if (file.is_open()) {
        // Save permutations
        file << "Permutations:\n";
        for (const auto& perm : permutations) {
            for (int step : perm) {
                file << step << " ";
            }
            file << "\n";
        }

        // Save coefficients
        file << "Coefficients:\n";
        for (const auto& coeff : coefficients) {
            const auto& divdiff = std::get<0>(coeff);
            for (const auto& elements : divdiff) {
                for (const auto& elem : elements) {
                    file << elem.energy << " " << elem.omega << " " << elem.chi << " " << elem.coefficient << " ";
                }
                file << "\n";
            }
            std::complex<double> energy_coefficient = std::get<1>(coeff);
            file << "Energy coefficient real:\n";
            file << energy_coefficient.real() << "\n"; 
            file << "Energy coefficient imag:\n";
            file << energy_coefficient.imag() << "\n";
            file << "Final state:\n";
            file << std::get<2>(coeff) << "\n";
        }

        file.close();
    } else {
        std::cerr << "Unable to open file for writing: " << filename << std::endl;
    }
}

// Combined function to generate permutations, calculate coefficients, and save results
void generate_and_save_results(int total_steps, int target_step, int moves, int start_state, const std::string& filename) {
    // Generate permutations
    std::vector<std::vector<int>> permutations = ladder_permutations(total_steps, target_step, moves);

    // Calculate coefficients for each permutation
    std::vector<std::tuple<std::vector<std::vector<DivdiffElement>>, std::complex<double>, int>> coefficients;
    for (const std::vector<int>& perm : permutations) {
        coefficients.push_back(cal_coefficient(perm, start_state));
    }

    // Save results to file
    save_results(permutations, coefficients, filename);
}

// Function to load results from a file
std::tuple<std::vector<std::vector<int>>, std::vector<std::tuple<std::vector<std::vector<DivdiffElement>>, std::complex<double>, int>>> load_results(const std::string& filename) {
    std::vector<std::vector<int>> permutations;
    std::vector<std::tuple<std::vector<std::vector<DivdiffElement>>, std::complex<double>, int>> coefficients;

    std::ifstream file(filename, std::ios::in);
    //print file name
    // std::cout << "FILE NAME IN SERIALIZER:" << filename << std::endl;
    if (file.is_open()) {
        std::string line;
        bool reading_permutations = true;

        while (std::getline(file, line)) {
            if (line == "Permutations:") {
                reading_permutations = true;
                continue;
            } else if (line == "Coefficients:") {
                reading_permutations = false;
                continue;
            }

            if (reading_permutations) {
                std::istringstream iss(line);
                std::vector<int> perm;
                int step;
                while (iss >> step) {
                    perm.push_back(step);
                }
                permutations.push_back(perm);
            } else {
                // Read divdiff elements
                std::vector<std::vector<DivdiffElement>> divdiff;
                
                //first divdiff element
                std::istringstream iss(line);
                std::vector<DivdiffElement> elements;
                DivdiffElement elem;
                while (iss >> elem.energy >> elem.omega >> elem.chi >> elem.coefficient) {
                    elements.push_back(elem);
                }
                divdiff.push_back(elements);

                //rest of divdiff elements
                while (std::getline(file, line) && line != "Energy coefficient real:") {
                    //print line
                    // std::cout << "LINE IN SERIALIZER:" << line << std::endl;

                    std::istringstream iss(line);
                    std::vector<DivdiffElement> elements;
                    DivdiffElement elem;
                    while (iss >> elem.energy >> elem.omega >> elem.chi >> elem.coefficient) {
                        elements.push_back(elem);
                    }
                    divdiff.push_back(elements);
                }

                // Read energy coefficient
                int real_part=0, imag_part=0;
                std::getline(file, line);
                // std::cout << "LINE IN SERIALIZER:" << line << std::endl;
                std::istringstream iss1(line);
                iss1 >> real_part;
                // if ("Energy coefficient real:" == line){
                //     std::getline(file, line);
                //     //print line
                    
                // }

                std::getline(file, line);
                if (line == "Energy coefficient imag:")
                {
                    std::getline(file, line);
                    // std::cout << "LINE IN SERIALIZER:" << line << std::endl;
                    std::istringstream iss1(line);
                    iss1 >> imag_part;
                }
                
                // //print real and imag part
                // std::cout << "real_part coefficient: " << real_part << std::endl;
                // std::cout << "imag_part coefficient: " << imag_part << std::endl;
                
                std::complex<double> energy_coefficient(real_part, imag_part);
                
                // //print energy_coefficient
                // std::cout << "COEF IN SERIALIZER:" << std::endl;
                // std::cout << "Energy coefficient: " << energy_coefficient.real() << " + " << energy_coefficient.imag() << "i" << std::endl;

                // Read final state
                int final_state = -1;
                getline(file, line);
                if (line == "Final state:"){
                    std::getline(file, line);
                    std::istringstream iss2(line);
                    iss2 >> final_state;
                }
                

                coefficients.push_back(std::make_tuple(divdiff, energy_coefficient, final_state));
            }
        }

        file.close();
    } else {
        std::cerr << "Unable to open file for reading: " << filename << std::endl;
    }

    return std::make_tuple(permutations, coefficients);
}

// int main() {
//     int total_steps = 10;
//     int start_state = 0;
//     std::vector<std::tuple<int, int, std::string>> filenames;

//     // Loop over target steps and moves
//     for (int target_step = 0; target_step <= 7; ++target_step) {
//         for (int moves = 1; moves <= 12; ++moves) {
//             std::string filename = "permutations\\dispersive_perms\\perm_" + std::to_string(target_step) + "_" + std::to_string(moves) + ".txt";
//             filenames.emplace_back(target_step, moves, filename);

//             // Check if file exists
//             if (fs::exists(filename)) {
//                 std::cout << "Loading results from " << filename << "\n";
                
//                 auto results = load_results(filename);
//                 // Optionally process the results as needed
//             } else {
//                 std::cout << "Generating results for target_step=" << target_step << " and moves=" << moves << "\n";
//                 generate_and_save_results(total_steps, target_step, moves, start_state, filename);
//             }
//         }
//     }

//     // Now filenames contains all the filenames that were checked
//     for (const auto& [target_step, moves, filename] : filenames) {
//         auto results = load_results(filename);
//         // Process the loaded results as needed
//         std::cout << "Processed results from " << filename << "\n";
//     }

//     return 0;
// }