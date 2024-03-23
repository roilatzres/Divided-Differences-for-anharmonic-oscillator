#include <iostream>
#include <vector>
#include <functional>
#include <cmath>
#include<complex>
#include "divdiffcomplex.h"
#include "permutation.h"

#define  MAX_STEP 10
#define  MAX_TARGET 6


int main() {
    // ***********************************Example usage for coefficient calculation*********************************
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

   // ***********************************Example usage for divided differences*********************************
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

//********************************************************************************************************************//



    //original parameters
    double multiplier = 4;
    double sigma = 48;
    double timestep = 1;
    double chi = -279e-6 * 2 * M_PI; 
    double detuning = chi / 100;

    double tfinal = multiplier * sigma;
    double omega = (2 * M_PI) / tfinal;

    std::vector<double> times;
    for (double t = 1; t <= tfinal - timestep; t += timestep) {//TODO: check if t=1 is correct (t=0 gives 0 amplitude)
        times.push_back(t);
    }


    //calculate the time evolution
    double t;

    // vector of all transition amplitudes for all times: [time][target_state][amp]
    std::vector<std::vector<std::complex<double>>> all_transition_amplitudes;
    
    //define start state
    int start_state = 0;

    //loop over all times
    // for (int i = 0; i < times.size(); i++) { 
    for (int i = 0; i < 2; i++) { //*************************** revert to all times************************
        t = times[i];
        //print time
        std::cout << "Time: " << t << std::endl;

        //vector of all transition amplitude for a given time
        std::vector<std::complex<double>> transition_amplitudes;
        
        //loop over all target steps
        // for (int target = 0; target <= MAX_TARGET; target++) {
        for (int target = 0; target <= 2; target++) {//*************************** revert to all targets************************
            std::complex<double> total_amp = 0;

            //print target
            std::cout << "Target: " << target << std::endl;
            
            //ladder_permutations(max_step, target, Q = total moves)
            // auto permutations = ladder_permutations(MAX_STEP, i, 6); //*************************** revert to Q = 6************************
            auto permutations = ladder_permutations(MAX_STEP, target, 2); //Q = 2
            
            //there are no permutaions - skip
            if (permutations.size() == 0) {
                continue;
            }

            //print all permutations
            std::cout << "Possible permutations:" << std::endl;
            for (const auto& combination : permutations) {
                for (int step : combination) {
                    std::cout << step << " ";
                }
                std::cout << std::endl;
            }

            // loop over all permutations
            for (const auto& permutation : permutations) {
                std::complex<double> perm_amp = 0;

                //calculate the coefficient
                auto result = cal_coefficient(permutation, start_state);
                auto beta = std::get<0>(result);
                std::complex<double> energy_coefficient = std::get<1>(result);
                int final_state = std::get<2>(result);

                //print beta
                std::cout << "Beta: " << std::endl;
                for (const auto& combination : beta) {
                    for (const auto& divdiff : combination) {
                        std::cout << "(" << divdiff.energy << ", " << divdiff.omega << ", " << divdiff.chi << ", " << divdiff.coefficient << ") ";
                    }
                    std::cout << std::endl;
                }
                std::cout << "Coefficient: " << energy_coefficient << std::endl;
                std::cout << "Final state: " << final_state << std::endl;

                //calculate divided differences
                divdiff_init();

                divdiffcomplex d(MAX_STEP,500);

                //calculate the divided differences for each path
                 for (const auto& combination : beta) {
                    double beta_coefficient = combination[0].coefficient;
                    // std::cout << "Beta coefficient: " << beta_coefficient << std::endl;

                    d.CurrentLength=0;
                    for (const auto& divdiff : combination) {
                        std::complex<double> n = (-1i * t * (divdiff.energy + divdiff.omega * omega  + divdiff.chi * chi));
                        //print n
                        std::cout << "n: " << n << std::endl;
                        d.AddElement(n);
                    }
                    int q = d.CurrentLength - 1;
                    double real = d.divdiffs[q].real().get_double();
                    double imag = d.divdiffs[q].imag().get_double();
                    
                    //print real and imag
                    std::cout << "Divdiff = " << real << " + " << imag << "i" << std::endl;

                    //print divdiffs elements
                    d.PrintList(d.divdiffs, d.CurrentLength, "Divdiffs");

                    //normalize the divided differences by (-it)^q / q! * beta_coefficient
                    std::complex<double> norm_real = cal_neg_i_pow(q) * std::pow(t, q) * real / factorial(q) ;
                    std::complex<double> norm_imag = cal_neg_i_pow(q+1) * std::pow(t, q) * (-imag) / factorial(q);// (-it)^q / q! * Imag = (-it)^q / q! * (i) * Imag.real = (-i)^q+1 * t^q * -Imag.real / q!
                    perm_amp = perm_amp + (norm_real + norm_imag) * beta_coefficient;
                    
                    //print perm_amp
                    std::cout << "Perm amplitude: " << perm_amp << std::endl;

                }
                total_amp += perm_amp * energy_coefficient;
                //print the total amplitude
                std::cout << "Total amplitude: " << total_amp << std::endl;

            }
            divdiff_clear_up();
            transition_amplitudes.push_back(total_amp);

        }
        all_transition_amplitudes.push_back(transition_amplitudes);

        //print all transition amplitudes
        std::cout << "All transition amplitudes: " << std::endl;
        for (const auto& transition : transition_amplitudes) {
            std::cout << transition << std::endl;
        }

    }

    return 0;
}





