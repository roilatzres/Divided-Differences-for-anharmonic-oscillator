#include <iostream>
#include <vector>
#include <functional>
#include <cmath>
#include <complex>
#include "include\divdiffcomplex.h"
#include "include\permutation.h"
#include <fstream>
#include "include\json-develop\single_include\nlohmann\json.hpp"


#define  MAX_STEP 10
#define  MAX_TARGET 7


int main() {
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

// ******************************* Example usage for permutation generation*********************************
    int steps = 5;
    std::pair<int, int> start_step = {0, 0};
    std::pair<int, int> mid_step = {3, 0};
    std::pair<int, int> target_step_1 = {0, 0};
    std::pair<int, int> target_step_2 = {1, 0};

    std::pair<int, int> max_step = {4, 3};

    auto result_1 = ladder_permutations(start_step, mid_step, target_step_1, max_step, steps);
    auto result_2 = ladder_permutations(start_step, mid_step, target_step_2, max_step, steps);

    auto all_combinations_up = std::get<0>(result_1);
    auto all_combinations_down = std::get<1>(result_1);

    // Output all_combinations_up
    std::cout << "all_combinations_up_1:" << std::endl;
    for (const auto& combination : all_combinations_up) {
        for (std::pair<int, int> step : combination) {
            std::cout << "{" << step.first << ", " << step.second << "} ";
        }
        std::cout << std::endl;
    }

    // Output all_combinations_down
    std::cout << "all_combinations_down_1:" << std::endl;
    for (const auto& combination : all_combinations_down) {
        for (std::pair<int, int> step : combination) {
            std::cout << "{" << step.first << ", " << step.second << "} ";
        }
        std::cout << std::endl;
    }


    auto all_combinations_up2 = std::get<0>(result_2);
    auto all_combinations_down2 = std::get<1>(result_2);

    // Output all_combinations_up2
    std::cout << "all_combinations_up_2:" << std::endl;
    for (const auto& combination : all_combinations_up2) {
        for (std::pair<int, int> step : combination) {
            std::cout << "{" << step.first << ", " << step.second << "} ";
        }
        std::cout << std::endl;
    }

    // Output all_combinations_down2
    std::cout << "all_combinations_down_2:" << std::endl;
    for (const auto& combination : all_combinations_down2) {
        for (std::pair<int, int> step : combination) {
            std::cout << "{" << step.first << ", " << step.second << "} ";
        }
        std::cout << std::endl;
    }
return 0;

}