#include <iostream>
#include <vector>
#include <functional>
#include <cmath>
#include <complex>
#include <fstream>
#include <string>
#include "include\divdiffcomplex.h"
#include "include\permutation.h"
#include "include\json-develop\single_include\nlohmann\json.hpp"
#include "include\Serializer.h"


#define  MAX_STEP 10
#define  MAX_TARGET 7


// int main() {
//     // ***********************************Example usage for coefficient calculation*********************************
//     std::vector<int> permutation = {1, 1};
//     int start_state = 1;
//     auto result = cal_coefficient(permutation, start_state);
//     auto beta = std::get<0>(result);
//     double coefficient = std::get<1>(result);
//     int final_state = std::get<2>(result);

//     // Output result
//     //print beta
//     std::cout << "Beta: " << std::endl;
//     for (const auto& combination : beta) {
//         for (const auto& divdiff : combination) {
//             std::cout << "(" << divdiff.energy << ", " << divdiff.omega << ") ";
//         }
//         std::cout << std::endl;
//     }
//     std::cout << "Coefficient: " << coefficient << std::endl;
//     std::cout << "Final state: " << final_state << std::endl;


//     double multiplier = 4;
//     double sigma = 48;
//     double timestep = 1;
//     double chi = -279e-6 * 2 * M_PI; // Assigning the value to chi
//    double detuning = chi / 100 ;

// //    ***********************************Example usage for divided differences*********************************
//    divdiff_init();

//    for (const auto& path : beta) {
//        divdiffcomplex d(MAX_STEP,500);

//        d.CurrentLength=0;

//        for (const auto& element : path) {
//            std::cout << "(" << element.energy << ", " << element.omega << ") ";
//            std::complex<double> n= -1i * (element.energy); // TODO: add time dependence, currently only for t=1
//            d.AddElement(n);
//        }
//        std::cout << std::endl << "finished path" << std::endl ;
//        //print divdiffs elements
//         d.PrintList(d.divdiffs, d.CurrentLength, "Divdiffs");
       
//         //  for (int i = 0; i < d.CurrentLength; i++) {
//         //         std::cout << std::endl;
//         //  }

//         std::cout << std::endl << "final element" << std::endl ;
//         std::cout << d.divdiffs[d.CurrentLength -1].real().print() << " + " << d.divdiffs[d.CurrentLength -1].imag().print() << std::endl;
//    }


//     divdiff_clear_up();

// // ******************************* Example usage for permutation generation*********************************
//     int total_steps = 10;
//     int target_step = 4;
//     int moves = 6;

//     auto result = ladder_permutations(total_steps, target_step, moves);

//     // Output result
//     std::cout << "Possible permutations:" << std::endl;
//     for (const auto& combination : result) {
//         for (int step : combination) {
//             std::cout << step << " ";
//         }
//         std::cout << std::endl;
//     }

// // ******************************* Example usage for OLD serializer*********************************
//     int total_steps = 100;
//     int start_state = 0;

//     // Create a 3D vector to store the results
//     std::vector<std::vector<std::vector<
//     std::tuple<std::vector<std::vector<int>>, std::vector<std::tuple<std::vector<std::vector<DivdiffElement>>, std::complex<double>, int>>>
//     >>> results;

//     // Loop over target steps and moves
//     for (int target_step = 0; target_step <= 11; ++target_step) {
//         // Add a new 2D vector for this target step
//         results.push_back(std::vector<std::vector<
//             std::tuple<std::vector<std::vector<int>>, std::vector<std::tuple<std::vector<std::vector<DivdiffElement>>, std::complex<double>, int>>>
//         >>());

//         for (int q = 0; q <= 12; ++q) {
//             std::string filename = "C:\\Users\\Latzres\\OneDrive - Technion\\RoeeL\\Divided-Differences-for-anharmonic-oscillator\\include\\permutations\\dispersive_perms\\dispersive_perms_" 
//                 + std::to_string(target_step) + "_" + std::to_string(q) + ".txt";

//             // Check if file exists
//             std::ifstream infile(filename);
//             if (!infile) {
//                 std::cout << "File " << filename << " does not exist." << std::endl;
//                 // Generate and save permutations 
//                 generate_and_save_results(total_steps, target_step, q, start_state, filename);
//             } else {

//                 // Add a new vector for this q
//                 results[target_step].push_back(std::vector<
//                     std::tuple<std::vector<std::vector<int>>, std::vector<std::tuple<std::vector<std::vector<DivdiffElement>>, std::complex<double>, int>>>
//                 >());

//                 // Load results from file
//                 std::tuple<std::vector<std::vector<int>>, std::vector<std::tuple<std::vector<std::vector<DivdiffElement>>, std::complex<double>, int>>> value = load_results(filename);
//                 results[target_step][q].push_back(value);
//             }
//         }
//     }

//     // // Output results
//     // for (int target_step = 0; target_step <= 7; ++target_step) {
//     //     for (int q = 0; q <= 8; ++q) {
//     //         std::cout << "Target step: " << target_step << ", q: " << q << std::endl;
//     //         auto value = results[target_step][q][0];
//     //         auto permutations = std::get<0>(value);
//     //         auto coefficients = std::get<1>(value);

//     //         std::cout << "Permutations:" << std::endl;
//     //         for (const auto& combination : permutations) {
                
//     //             for (int step : combination) {
//     //                 std::cout << step << " ";
//     //             }
//     //             std::cout << std::endl;
//     //         }
//     //         // print permutation size
//     //         std::cout << "Permutation size: " << permutations.size() << std::endl;

//     //         std::cout << "Coefficients:" << std::endl;
//     //         for (const auto& coefficient : coefficients) {
//     //             auto divdiff = std::get<0>(coefficient);
//     //             auto energy_coefficient = std::get<1>(coefficient);
//     //             auto final_state = std::get<2>(coefficient);

//     //             std::cout << "Energy coefficient: " << energy_coefficient.real() << " + " << energy_coefficient.imag() << "i" << std::endl;
//     //             std::cout << "Final state: " << final_state << std::endl;
//     //             std::cout << "Divdiff elements:" << std::endl;
//     //             for (const auto& elements : divdiff) {
//     //                 for (const auto& element : elements) {
//     //                     std::cout << "(" << element.energy << ", " << element.omega << ", " << element.chi << ", " << element.coefficient << ") ";
//     //                 }
//     //                 std::cout << std::endl;
//     //             }
//     //         }
//     //     }
//     // }


// // ******************************* Example usage for NEW serializer*********************************

#include <iostream>
#include <vector>
#include <tuple>
#include <string>

using namespace std;

void print_permutations(const vector<vector<int>>& permutations) {
    for (const auto& permutation : permutations) {
        for (int step : permutation) {
            cout << step << " ";
        }
        cout << endl;
    }
}

void print_coefficients(const vector<tuple<vector<DivdiffElement>, double, int, vector<int>>>& loaded_coefficients) {
     for (const auto& coefficients : loaded_coefficients) {
        const auto& [beta, energy_coefficient_real, final_state, states] = coefficients;

        cout << "States: ";
        for (int state : states) {
            cout << state << " ";
        }
        cout << endl;
       
        cout << "Divdiff Elements:" << endl;
        for (const auto& elem : beta) {
            cout << "Omega: ";
            for (int omega : elem.omega) {
                cout << omega << " ";
            }
            cout << endl;
            cout << "Coefficient: " << elem.coefficient << endl;
        }

        cout << "Energy Coefficient (Real): " << energy_coefficient_real << endl;
        cout << "Final State: " << final_state << endl;

     }
}

//original parameters
double multiplier = 4;
float sigma = 48;
float timestep = 1;
double chi = -279e-6 * 2 * M_PI; 
// double detuning = chi *1e-2;

double tfinal = multiplier * sigma;
double omega = (2 * M_PI) / tfinal;


int main() {
    //check amount of permutation
    int total_steps = 100;
    int target_first_step = 0;
    int target_second_step = 0;
    int moves = 16;

    
    // print num permutation for each target step and q
    for (int target_step = 0; target_step <= 8; ++target_step) {
        for (int q = 0; q <= 14; ++q) {
            cout << "Target step: " << target_step << ", q: " << q << endl;
            auto result = jc_ladder_permutations(target_step, target_second_step, q);
            cout << "Number of permutations: " << result.size() << endl << endl;
        }
    }
    

    // /* test permutation loading */

    // int total_steps = 100;
    // int start_state = 0;

    // // Create a 3D vector to store the results
    // vector< //target_step
    // vector< //q
    // vector< //permutations
    // std::tuple<vector<DivdiffElement>, double, int, vector<int>>> // coefficients
    // >> results;

    // cout << "Loading permutations from file..." << endl;

    // // Loop over target steps and moves
    // for (int target_step = 0; target_step <= 8; ++target_step) {//FIXME: maybe 8 instead of 9
    //     results.push_back(vector<vector<std::tuple<vector<DivdiffElement>, double, int, vector<int>>>>());

    //     for (int q = 0; q <= 14; ++q) {

    //         std::string filename = "C:\\Users\\Latzres\\OneDrive - Technion\\RoeeL\\Divided-Differences-for-anharmonic-oscillator\\include\\parms\\dispersive\\coeff\\dispersive_coeff_" 
    //             + std::to_string(target_step) + "_" + std::to_string(q) + ".bin";

    //         // Check if file exists
    //         std::ifstream infile(filename);
    //         if (!infile) {
    //             std::cout << "File " << filename << " does not exist." << std::endl;
    //             // generate and save coefficients
    //             generate_and_save_coefficients(start_state, total_steps, target_step, q, filename);

    //         } else {
                
    //             //print loading
    //             cout << "Loading target step: " << target_step << ", q: " << q << endl;
    //             // Add a new vector for this q
    //             results[target_step].push_back(vector<std::tuple<vector<DivdiffElement>, double, int, vector<int>>>());
                


    //             // Load results from file
    //             vector<std::tuple<vector<DivdiffElement>, double, int, vector<int>>> value = load_coefficients(filename);
    //             results[target_step][q] = value;
                                
    //         }
    //     }
    // }

    // cout << "Permutations loaded successfully." << endl;

    // // // Output results
    // // for (int target_step = 0; target_step <= 1; ++target_step) {
    // //     for (int q = 0; q <= 4; ++q) {
    // //         std::cout << "Target step: " << target_step << ", q: " << q << std::endl;
    // //         auto value = results[target_step][q];
    // //         print_coefficients(value);
    // //     }
    // // }

    return 0;
}

// return 0;

// }