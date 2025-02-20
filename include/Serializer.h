#include <iostream>
#include <vector>
#include <tuple>
#include <complex>
#include <fstream>
#include <functional>
#include <sstream>
#include <cmath>
#include <cstring>
#include "permutation.h"

using std::vector;
using std::complex;

void save_binary_perm(const std::vector<std::vector<int>>& data, std::string& filename) {
    std::ofstream file(filename, std::ios::binary);
    if (!file) {
        std::cerr << "Error opening file for writing!" << std::endl;
        return;
    }

    size_t rows = data.size();
    size_t cols = (rows > 0) ? data[0].size() : 0;

    file.write(reinterpret_cast<const char*>(&rows), sizeof(rows));
    file.write(reinterpret_cast<const char*>(&cols), sizeof(cols));

    for (const auto& row : data) {
        if (!row.empty()) {  // Ensure valid data pointer
            file.write(reinterpret_cast<const char*>(row.data()), row.size() * sizeof(int));
        }
    }

    file.close();
}

vector<vector<int>> load_binary_perm(std::string& filename) {
    std::ifstream file(filename, std::ios::binary);
    if (!file) {
        std::cerr << "Error opening file for reading!" << std::endl;
        return {};
    }

    size_t rows, cols;
    file.read(reinterpret_cast<char*>(&rows), sizeof(rows));
    file.read(reinterpret_cast<char*>(&cols), sizeof(cols));

    vector<vector<int>> data(rows, vector<int>(cols));

    for (auto& row : data) {
        if (cols > 0) {  // Ensure valid memory
            file.read(reinterpret_cast<char*>(row.data()), cols * sizeof(int));
        }
    }

    file.close();
    return data;
}


void save_binary_coef(const std::vector<std::tuple<double, int, std::vector<int>>>& data, const std::string& filename) {
    std::ofstream file(filename, std::ios::binary);
    if (!file) {
        std::cerr << "Error opening file for writing: " << filename << std::endl;
        return;
    }

    size_t num_entries = data.size();
    file.write(reinterpret_cast<const char*>(&num_entries), sizeof(num_entries));

    std::vector<char> buffer;
    buffer.reserve(num_entries * (sizeof(double) + sizeof(int) + sizeof(size_t))); // Pre-allocate memory

    for (const auto& entry : data) {
        double val;
        int num;
        std::vector<int> vec;
        std::tie(val, num, vec) = entry;

        buffer.insert(buffer.end(), reinterpret_cast<const char*>(&val), reinterpret_cast<const char*>(&val) + sizeof(val));
        buffer.insert(buffer.end(), reinterpret_cast<const char*>(&num), reinterpret_cast<const char*>(&num) + sizeof(num));

        size_t vec_size = vec.size();
        buffer.insert(buffer.end(), reinterpret_cast<const char*>(&vec_size), reinterpret_cast<const char*>(&vec_size) + sizeof(vec_size));

        if (!vec.empty()) {
            buffer.insert(buffer.end(), reinterpret_cast<const char*>(vec.data()), reinterpret_cast<const char*>(vec.data()) + vec_size * sizeof(int));
        }
    }

    file.write(buffer.data(), buffer.size()); // Write all at once
    file.close();
}

std::vector<std::tuple<double, int, std::vector<int>>> load_binary_coef(const std::string& filename) {
    std::ifstream file(filename, std::ios::binary);
    if (!file) {
        std::cerr << "Error opening file for reading: " << filename << std::endl;
        return {};
    }

    size_t num_entries;
    file.read(reinterpret_cast<char*>(&num_entries), sizeof(num_entries));

    std::vector<char> buffer(std::istreambuf_iterator<char>(file), {}); // Load file into memory
    file.close();

    std::vector<std::tuple<double, int, std::vector<int>>> data;
    data.reserve(num_entries);

    const char* ptr = buffer.data();
    const char* end = buffer.data() + buffer.size();

    for (size_t i = 0; i < num_entries; ++i) {
        if (ptr + sizeof(double) + sizeof(int) + sizeof(size_t) > end) break;

        double val = *reinterpret_cast<const double*>(ptr);
        ptr += sizeof(double);

        int num = *reinterpret_cast<const int*>(ptr);
        ptr += sizeof(int);

        size_t vec_size = *reinterpret_cast<const size_t*>(ptr);
        ptr += sizeof(size_t);

        if (ptr + vec_size * sizeof(int) > end) break;

        std::vector<int> vec(vec_size);
        std::memcpy(vec.data(), ptr, vec_size * sizeof(int));
        ptr += vec_size * sizeof(int);

        data.emplace_back(val, num, std::move(vec));
    }

    return data;
}

void save_binary_complexEx(complex_Ex& data, const std::string& filename) {
    std::ofstream file(filename, std::ios::binary);
    if (!file) {
        std::cerr << "Error opening file for writing: " << filename << std::endl;
        return;
    }

    double real_val = data.real.get_double();
    double imag_val = data.imag.get_double();

    file.write(reinterpret_cast<const char*>(&real_val), sizeof(double));
    file.write(reinterpret_cast<const char*>(&imag_val), sizeof(double));

    file.close();
}

complex_Ex load_binary_complexEx(const std::string& filename) {
    std::ifstream file(filename, std::ios::binary);
    if (!file) {
        std::cerr << "Error opening file for reading: " << filename << std::endl;
        return {};
    }

    double real_val, imag_val;
    file.read(reinterpret_cast<char*>(&real_val), sizeof(double));
    file.read(reinterpret_cast<char*>(&imag_val), sizeof(double));

    complex_Ex result;
    result.real = ExExFloat(real_val);
    result.imag = ExExFloat(imag_val);

    file.close();
    return result;
}




// // Function to save permutations and coefficients to a file
// void save_results(const vector<vector<int>>& permutations, const vector<std::tuple<vector<DivdiffElement>, double, int, vector<int>>>& coefficients, const std::string& filename) {
//     std::ofstream file(filename, std::ios::out);
//     if (file.is_open()) {
//         // Save permutations
//         file << "Permutations:\n";
//         for (const auto& perm : permutations) {
//             for (int step : perm) {
//                 file << step << " ";
//             }
//             file << "\n";
//         }

//         // Save coefficients
//         file << "Coefficients:\n";
//         for (const auto& coeff : coefficients) {
//             //print states from coeff
//             vector<int> states = std::get<3>(coeff);
//             std::cout << "States in coeff: ";
//             for (const auto& state : states){
//                 std::cout << state << " ";
//             }
//             std::cout << std::endl;

//             const auto& beta = std::get<0>(coeff);
//             for (const auto& elements : beta) {
                
//                 //print omega vector of each element
//                 std::cout << "Omega vector: ";
//                 for (const auto& omega : elements.omega){
//                     std::cout << omega << " ";
//                 }
//                 std::cout << std::endl;
//                 //print coefficient of each element
//                 std::cout << "Coefficient: " << elements.coefficient << std::endl;
//                 file << "\n";
//             }
//             double energy_coefficient = std::get<1>(coeff);
//             file << "Energy coefficient real:\n";
//             file << energy_coefficient << "\n"; 
//             file << "Final state:\n";
//             file << std::get<2>(coeff) << "\n";
//         }

//         file.close();
//     } else {
//         std::cerr << "Unable to open file for writing: " << filename << std::endl;
//     }
// }

// // Combined function to generate permutations, calculate coefficients, and save results
// void generate_and_save_results(int total_steps, int target_step, int moves, int start_state, const std::string& filename) {
//     // Generate permutations
//     vector<vector<int>> permutations = ladder_permutations(total_steps, target_step, moves);

//     // Calculate coefficients for each permutation
//     vector<std::tuple<vector<DivdiffElement>, double, int, vector<int>>> coefficients;
//     for (const vector<int>& perm : permutations) {
//         coefficients.push_back(cal_coefficient(perm, start_state));
//     }

//     // Save results to file
//     save_results(permutations, coefficients, filename);
// }

// // Function to load results from a file
// std::tuple<vector<vector<int>>, vector<std::tuple<vector<vector<DivdiffElement>>, complex<double>, int>>> load_results(const std::string& filename) {
//     vector<vector<int>> permutations;
//     vector<std::tuple<vector<vector<DivdiffElement>>, complex<double>, int>> coefficients;

//     std::ifstream file(filename, std::ios::in);
//     //print file name
//     // std::cout << "FILE NAME IN SERIALIZER:" << filename << std::endl;
//     if (file.is_open()) {
//         std::string line;
//         bool reading_permutations = true;

//         while (std::getline(file, line)) {
//             if (line == "Permutations:") {
//                 reading_permutations = true;
//                 continue;
//             } else if (line == "Coefficients:") {
//                 reading_permutations = false;
//                 continue;
//             }

//             if (reading_permutations) {
//                 std::istringstream iss(line);
//                 vector<int> perm;
//                 int step;
//                 while (iss >> step) {
//                     perm.push_back(step);
//                 }
//                 permutations.push_back(perm);
//             } else {
//                 // Read divdiff elements
//                 vector<vector<DivdiffElement>> divdiff;
                
//                 //first divdiff element
//                 std::istringstream iss(line);
//                 vector<DivdiffElement> elements;
//                 DivdiffElement elem;
//                 while (iss >> elem.energy >> elem.omega >> elem.chi >> elem.coefficient) {
//                     elements.push_back(elem);
//                 }
//                 divdiff.push_back(elements);

//                 //rest of divdiff elements
//                 while (std::getline(file, line) && line != "Energy coefficient real:") {
//                     //print line
//                     // std::cout << "LINE IN SERIALIZER:" << line << std::endl;

//                     std::istringstream iss(line);
//                     vector<DivdiffElement> elements;
//                     DivdiffElement elem;
//                     while (iss >> elem.energy >> elem.omega >> elem.chi >> elem.coefficient) {
//                         elements.push_back(elem);
//                     }
//                     divdiff.push_back(elements);
//                 }

//                 // Read energy coefficient
//                 int real_part=0, imag_part=0;
//                 std::getline(file, line);
//                 // std::cout << "LINE IN SERIALIZER:" << line << std::endl;
//                 std::istringstream iss1(line);
//                 iss1 >> real_part;
//                 // if ("Energy coefficient real:" == line){
//                 //     std::getline(file, line);
//                 //     //print line
                    
//                 // }

//                 std::getline(file, line);
//                 if (line == "Energy coefficient imag:")
//                 {
//                     std::getline(file, line);
//                     // std::cout << "LINE IN SERIALIZER:" << line << std::endl;
//                     std::istringstream iss1(line);
//                     iss1 >> imag_part;
//                 }
                
//                 // //print real and imag part
//                 // std::cout << "real_part coefficient: " << real_part << std::endl;
//                 // std::cout << "imag_part coefficient: " << imag_part << std::endl;
                
//                 complex<double> energy_coefficient(real_part, imag_part);
                
//                 // //print energy_coefficient
//                 // std::cout << "COEF IN SERIALIZER:" << std::endl;
//                 // std::cout << "Energy coefficient: " << energy_coefficient.real() << " + " << energy_coefficient.imag() << "i" << std::endl;

//                 // Read final state
//                 int final_state = -1;
//                 getline(file, line);
//                 if (line == "Final state:"){
//                     std::getline(file, line);
//                     std::istringstream iss2(line);
//                     iss2 >> final_state;
//                 }
                

//                 coefficients.push_back(std::make_tuple(divdiff, energy_coefficient, final_state));
//             }
//         }

//         file.close();
//     } else {
//         std::cerr << "Unable to open file for reading: " << filename << std::endl;
//     }

//     return std::make_tuple(permutations, coefficients);
// }

// int main() {
//     int total_steps = 10;
//     int start_state = 0;
//     vector<std::tuple<int, int, std::string>> filenames;

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



// ****************************************************************************************
// ***************** Save and Load Functions **********************************************
// ****************************************************************************************

void save_permutations(const vector<vector<int>>& permutations, const std::string& filename) {
    std::ofstream file(filename, std::ios::binary);
    if (!file) {
        throw std::runtime_error("Cannot open file for writing.");
    }

    size_t size = permutations.size();
    file.write(reinterpret_cast<const char*>(&size), sizeof(size));
    for (const auto& permutation : permutations) {
        size_t inner_size = permutation.size();
        file.write(reinterpret_cast<const char*>(&inner_size), sizeof(inner_size));
        file.write(reinterpret_cast<const char*>(permutation.data()), inner_size * sizeof(int));
    }
}

vector<vector<int>> load_permutations(const std::string& filename) {
    std::ifstream file(filename, std::ios::binary);
    if (!file) {
        throw std::runtime_error("Cannot open file for reading.");
    }

    size_t size;
    file.read(reinterpret_cast<char*>(&size), sizeof(size));
    vector<vector<int>> permutations(size);

    for (auto& permutation : permutations) {
        size_t inner_size;
        file.read(reinterpret_cast<char*>(&inner_size), sizeof(inner_size));
        permutation.resize(inner_size);
        file.read(reinterpret_cast<char*>(permutation.data()), inner_size * sizeof(int));
    }

    return permutations;
}


// void save_coefficients(const std::tuple<vector<DivdiffElement>, double, int, vector<int>>& coefficients, const std::string& filename) {
//     std::ofstream file(filename, std::ios::binary);
//     if (!file) {
//         throw std::runtime_error("Cannot open file for writing.");
//     }

//     const auto& [beta, energy_coefficient_real, final_state, states] = coefficients;

//     size_t beta_size = beta.size();
//     file.write(reinterpret_cast<const char*>(&beta_size), sizeof(beta_size));
//     for (const auto& elem : beta) {
//         size_t omega_size = elem.omega.size();
//         file.write(reinterpret_cast<const char*>(&omega_size), sizeof(omega_size));
//         file.write(reinterpret_cast<const char*>(elem.omega.data()), omega_size * sizeof(int));
//         file.write(reinterpret_cast<const char*>(&elem.coefficient), sizeof(elem.coefficient));
//     }

//     file.write(reinterpret_cast<const char*>(&energy_coefficient_real), sizeof(energy_coefficient_real));
//     file.write(reinterpret_cast<const char*>(&final_state), sizeof(final_state));

//     size_t states_size = states.size();
//     file.write(reinterpret_cast<const char*>(&states_size), sizeof(states_size));
//     file.write(reinterpret_cast<const char*>(states.data()), states_size * sizeof(int));
// }

// std::tuple<vector<DivdiffElement>, double, int, vector<int>> load_coefficients(const std::string& filename) {
//     std::ifstream file(filename, std::ios::binary);
//     if (!file) {
//         throw std::runtime_error("Cannot open file for reading.");
//     }

//     size_t beta_size;
//     file.read(reinterpret_cast<char*>(&beta_size), sizeof(beta_size));
//     vector<DivdiffElement> beta(beta_size);

//     for (auto& elem : beta) {
//         size_t omega_size;
//         file.read(reinterpret_cast<char*>(&omega_size), sizeof(omega_size));
//         elem.omega.resize(omega_size);
//         file.read(reinterpret_cast<char*>(elem.omega.data()), omega_size * sizeof(int));
//         file.read(reinterpret_cast<char*>(&elem.coefficient), sizeof(elem.coefficient));
//     }

//     double energy_coefficient_real;
//     file.read(reinterpret_cast<char*>(&energy_coefficient_real), sizeof(energy_coefficient_real));

//     int final_state;
//     file.read(reinterpret_cast<char*>(&final_state), sizeof(final_state));

//     size_t states_size;
//     file.read(reinterpret_cast<char*>(&states_size), sizeof(states_size));
//     vector<int> states(states_size);
//     file.read(reinterpret_cast<char*>(states.data()), states_size * sizeof(int));

//     return std::make_tuple(beta, energy_coefficient_real, final_state, states);
// }


//////////////////// coeff for vector

void save_coefficients(const vector<std::tuple<vector<DivdiffElement>, double, int, vector<int>>>& coefficients, const std::string& filename) {
    std::ofstream file(filename, std::ios::binary);
    if (!file) {
        throw std::runtime_error("Cannot open file for writing.");
    }
    size_t size = coefficients.size();
    file.write(reinterpret_cast<const char*>(&size), sizeof(size));
    for (const auto& coeff : coefficients) {
        // const auto& [beta, energy_coefficient_real, final_state, states] = coeff;
        const auto& beta = std::get<0>(coeff);
        const auto& energy_coefficient_real = std::get<1>(coeff);
        const auto& final_state = std::get<2>(coeff);
        const auto& states = std::get<3>(coeff);


        size_t beta_size = beta.size();
        file.write(reinterpret_cast<const char*>(&beta_size), sizeof(beta_size));
        for (const auto& elem : beta) {
            size_t omega_size = elem.omega.size();
            file.write(reinterpret_cast<const char*>(&omega_size), sizeof(omega_size));
            file.write(reinterpret_cast<const char*>(elem.omega.data()), omega_size * sizeof(int));
            file.write(reinterpret_cast<const char*>(&elem.coefficient), sizeof(elem.coefficient));
        }
        file.write(reinterpret_cast<const char*>(&energy_coefficient_real), sizeof(energy_coefficient_real));
        file.write(reinterpret_cast<const char*>(&final_state), sizeof(final_state));
        size_t states_size = states.size();
        file.write(reinterpret_cast<const char*>(&states_size), sizeof(states_size));
        file.write(reinterpret_cast<const char*>(states.data()), states_size * sizeof(int));
    }
}

vector<std::tuple<vector<DivdiffElement>, double, int, vector<int>>> load_coefficients(const std::string& filename) {
    std::ifstream file(filename, std::ios::binary);
    if (!file) {
        throw std::runtime_error("Cannot open file for reading.");
    }
    size_t size;
    file.read(reinterpret_cast<char*>(&size), sizeof(size));
    vector<std::tuple<vector<DivdiffElement>, double, int, vector<int>>> coefficients(size);
    for (auto& coeff : coefficients) {
        // auto& [beta, energy_coefficient_real, final_state, states] = coeff;
        auto& beta = std::get<0>(coeff);
        auto& energy_coefficient_real = std::get<1>(coeff);
        auto& final_state = std::get<2>(coeff);
        auto& states = std::get<3>(coeff);

        
        size_t beta_size;
        file.read(reinterpret_cast<char*>(&beta_size), sizeof(beta_size));
        beta.resize(beta_size);
        for (auto& elem : beta) {
            size_t omega_size;
            file.read(reinterpret_cast<char*>(&omega_size), sizeof(omega_size));
            elem.omega.resize(omega_size);
            file.read(reinterpret_cast<char*>(elem.omega.data()), omega_size * sizeof(int));
            file.read(reinterpret_cast<char*>(&elem.coefficient), sizeof(elem.coefficient));
        }
        file.read(reinterpret_cast<char*>(&energy_coefficient_real), sizeof(energy_coefficient_real));
        file.read(reinterpret_cast<char*>(&final_state), sizeof(final_state));
        size_t states_size;
        file.read(reinterpret_cast<char*>(&states_size), sizeof(states_size));
        states.resize(states_size);
        file.read(reinterpret_cast<char*>(states.data()), states_size * sizeof(int));
    }
    return coefficients;
}


// Function to generate all permutations and save all coefficients to a file
void generate_and_save_coefficients(int start_step, int total_steps, int target_step, int moves, const std::string& filename) {
    vector<vector<int>> permutations = ladder_permutations(start_step, total_steps, target_step, moves);
    int start_state = 0;
    vector<std::tuple<vector<DivdiffElement>, double, int, vector<int>>> all_coefficients;

    for (const auto& permutation : permutations) {
        auto coefficients = cal_coefficient(permutation, start_state);
        all_coefficients.push_back(coefficients);
    }

    save_coefficients(all_coefficients, filename);
}

// save complex number to file
void save_complex_to_file(const complex<double>& number, const std::string& filename) {
    std::ofstream file(filename, std::ios::binary);
    if (!file) {
        throw std::runtime_error("Cannot open file for writing.");
    }
    file.write(reinterpret_cast<const char*>(&number), sizeof(number));
}

// load complex number from file
complex<double> load_complex_from_file(const std::string& filename) {
    std::ifstream file(filename, std::ios::binary);
    if (!file) {
        throw std::runtime_error("Cannot open file for reading.");
    }
    complex<double> number;
    file.read(reinterpret_cast<char*>(&number), sizeof(number));
    return number;
}


void load_all_coefficients(std::vector<std::vector<std::vector<std::tuple<std::vector<DivdiffElement>, double, int, std::vector<int>>>>> &all_coefficients, int start_state, int total_steps, int max_target, int q_max)
{
    // Loop over target steps and moves
    for (int target_step = 0; target_step <= max_target; ++target_step)
    {
        all_coefficients.push_back(vector<vector<std::tuple<vector<DivdiffElement>, double, int, vector<int>>>>());

        for (int q = 0; q <= q_max; ++q)
        {
            std::string filename = "./include/parms/dispersive/coeff/dispersive_coeff_" + std::to_string(target_step) + "_" + std::to_string(q) + ".bin";

            // Check if file exists
            std::ifstream infile(filename);
            if (!infile)
            {
                std::cout << "File " << filename << " does not exist." << std::endl;
                // generate and save coefficients
                generate_and_save_coefficients(start_state, total_steps, target_step, q, filename);
            }
            else
            {

                // Add a new vector for this q
                all_coefficients[target_step].push_back(vector<std::tuple<vector<DivdiffElement>, double, int, vector<int>>>());

                // Load all_coefficients from file
                vector<std::tuple<vector<DivdiffElement>, double, int, vector<int>>> value = load_coefficients(filename);
                all_coefficients[target_step][q] = value;
            }
        }
    }
}
