#include <iostream>
#include <vector>
#include <functional>
#include <cmath>
#include "permutation.h"
#include<complex>

#define  MAX_STEP 10
#define  MAX_TARGET 6


int main() {
    // Example usage for coefficient calculation
    std::vector<int> permutation = {1};
    int start_state = 1;
    auto result = cal_coefficient(permutation, start_state);
    auto beta = std::get<0>(result);
    double coefficient = std::get<1>(result);
    int final_state = std::get<2>(result);

    // Output result
    //print beta
    std::cout << "Beta: " << std::endl;
    for (const auto& combination : beta) {
        for (const auto& divdiff : combination) {
            std::cout << "(" << divdiff.energy << ", " << divdiff.omega << ") ";
        }
        std::cout << std::endl;
    }
    std::cout << "Coefficient: " << coefficient << std::endl;
    std::cout << "Final state: " << final_state << std::endl;

    // double multiplier = 4;
    // double sigma = 48;
    // double timestep = 1;
    // double chi = -279e-6 * 2 * M_PI; // Assigning the value to chi
    double detuning = chi / 100;

    //calculate divided differences
    divdiff_init();

    divdiffcomplex d(MAX_STEP,500);

    d.CurrentLength=0;
    for (const auto& combination : beta) {
        for (const auto& divdiff : combination) {
            std::cout << "(" << divdiff.energy << ", " << divdiff.omega << ") ";
            std::complex<double> n= -i * (divdiff.energy + divdiff.omega * detuning); // TODO: add time dependence, currently only for t=1
            d.AddElement(n);
        }
    }

    //pull the last element from the divdiffs
    std::complex<double> last = d.divdiffs[d.CurrentLength-1][0];

    divdiff_clear_up();



//********************************************************************************************************************//



    //original parameters
    double multiplier = 4;
    double sigma = 48;
    double timestep = 1;
    double chi = -279e-6 * 2 * M_PI; 
    double detuning = chi / 100;

    double tfinal = multiplier * sigma;

    std::vector<double> times;
    for (double t = 0; t < tfinal - timestep; t += timestep) {
        times.push_back(t);
    }


    //calculate the time evolution
    double t;

    // vector of all transition amplitudes for all times: [time][state][amp]
    std::vector<std::vector<std::complex<double>>> all_transition_amplitudes;
    
    //define start state
    int start_state = 0;

    //loop over all times
    for (int i = 0; i < times.size(); i++) {
        t = times[i];
    
        //vector of all transition amplitude for a given time
        std::vector<std::complex<double>> transition_amplitudes;
        
        //loop over all target steps
        for (int target = 0; target <= MAX_TARGET; target++) {
            std::complex<double> amp = 0;
            
            //ladder_permutations(max_step, target, Q)
            auto permutations = ladder_permutations(MAX_STEP, i, 6); //Q = 6

            // loop over all permutations
            for (const auto& permutation : permutations) {
                //calculate the coefficient
                auto result = cal_coefficient(permutation, start_state);
                auto beta = std::get<0>(result);
                double coefficient = std::get<1>(result);
                int final_state = std::get<2>(result);

                //calculate divided differences
                divdiff_init();

                divdiffcomplex d(MAX_STEP,500);

                d.CurrentLength=0;
                for (const auto& combination : beta) {
                    for (const auto& divdiff : combination) {
                        std::cout << "(" << divdiff.energy << ", " << divdiff.omega << ") ";
                        std::complex<double> n= -i * t *(divdiff.energy + divdiff.omega * detuning); 
                        d.AddElement(n);
                    }
                }

                //pull the last element from the divdiffs
                std::complex<double> last = d.divdiffs[d.CurrentLength-1][0];
                amp += coefficient * last;

                divdiff_clear_up();
            }
            transition_amplitudes.push_back(amp);
        }
        all_transition_amplitudes.push_back(transition_amplitudes);

    }

    return 0;
}


// // Example usage for permutation generation
//     int total_steps = 10;
//     int target_step = 4;
//     int moves = 6;
//
//     auto result = ladder_permutations(total_steps, target_step, moves);
//
//     // Output result
//     std::cout << "Possible permutations:" << std::endl;
//     for (const auto& combination : result) {
//         for (int step : combination) {
//             std::cout << step << " ";
//         }
//         std::cout << std::endl;
//     }
