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
    float sigma = 48;
    float timestep = 1;
    double chi = -279e-6 * 2 * M_PI; 
    // double detuning = chi *1e-2;

    //print vars
    std::cout << "Multiplier: " << multiplier << std::endl;
    std::cout << "Sigma: " << sigma << std::endl;
    std::cout << "Timestep: " << timestep << std::endl;
    std::cout << "Chi: " << chi << std::endl;
    // std::cout << "Detuning: " << detuning << std::endl;


    double tfinal = multiplier * sigma;
    double omega = (2 * M_PI) / tfinal;

    std::cout << "tfinal: " << tfinal << std::endl;
    std::cout << "Omega: " << omega << std::endl;
    std::vector<double> times;
    for (double t = timestep; t <= tfinal; t += timestep) {//TODO: check if t=1 is correct (t=0 gives 0 amplitude)
        times.push_back(t);
        //print times
        // std::cout << "Time: " << t << std::endl;
    }


    //calculate the time evolution
    double t;

    // vector of all transition ampli tudes for all times: [time][target_state][amp]
    std::vector<std::vector<std::complex<double>>> all_transition_amplitudes;
    
    //define start state
    int start_state = 0;

    int qubit = 1;
    // for (int qubit = 0; qubit <=1; qubit++){
        //loop over all times
        for (int i = 0; i < times.size(); i++) { 
            t = times[i];
            //print time
            std::cout << std::endl << "Time: " << t << std::endl;

            //vector of all transition amplitude for a given time
            std::vector<std::complex<double>> transition_amplitudes;
            
            //loop over all target steps
            for (int target = 0; target <= MAX_TARGET; target++) {
                std::complex<double> total_amp = 0;

                //print target
                std::cout << "Target: " << target << std::endl;

                for(int q = 1; q < 8; q++){

                    std::complex<double> q_amp = 0;
                    
                    //ladder_permutations(max_step, target, Q = total moves)
                    auto permutations = ladder_permutations(MAX_STEP, target, q);
                    
                    //there are no permutaions - skip
                    if (permutations.size() == 0) {
                        continue;
                    }

                    // //print all permutations
                    // std::cout << "Possible permutations:" << std::endl;
                    // for (const auto& combination : permutations) {
                    //     for (int step : combination) {
                    //         std::cout << step << " ";
                    //     }
                    //     std::cout << std::endl;
                    // }

                    // loop over all permutations
                    for (const auto& permutation : permutations) {
                        std::complex<double> perm_amp = 0;

                        //calculate the coefficient
                        auto result = cal_coefficient(permutation, start_state);
                        auto beta = std::get<0>(result);
                        std::complex<double> energy_coefficient = std::get<1>(result);
                        int final_state = std::get<2>(result);

                        // //print beta
                        // std::cout << "Beta: " << std::endl;
                        // for (const auto& combination : beta) {
                        //     for (const auto& divdiff : combination) {
                        //         std::cout << "(" << divdiff.energy << ", " << divdiff.omega << ", " << divdiff.chi << ", " << divdiff.coefficient << ") ";
                        //     }
                        //     std::cout << std::endl;
                        // }
                        // std::cout << "Coefficient: " << energy_coefficient << std::endl;
                        // std::cout << "Final state: " << final_state << std::endl;

                        //calculate divided differences
                        divdiff_init();

                        divdiffcomplex d(MAX_STEP,500);

                        //calculate the divided differences for each path
                        for (const auto& combination : beta) {
                            double beta_coefficient = combination[0].coefficient;
                            // std::cout << "Beta coefficient: " << beta_coefficient << std::endl;

                            d.CurrentLength=0;
                            for (const auto& divdiff : combination) {
                                std::complex<double> n;
                                if(qubit == 0){
                                    n = (-1j * t * (divdiff.omega * omega  + divdiff.chi * chi));
                                }else{//qubit == 1
                                    n = (-1j * t * (-divdiff.energy * chi + divdiff.omega * omega  + divdiff.chi * chi));
                                }
                                // //print n
                                // std::cout << "n: " << n << std::endl;
                                d.AddElement(n);
                            }
                            int q = d.CurrentLength - 1;
                            double real = d.divdiffs[q].real().get_double();
                            double imag = d.divdiffs[q].imag().get_double();
                            
                            //print real and imag
                            // std::cout << "Divdiff = " << real << " + " << imag << "i" << std::endl;

                            // //print divdiffs elements
                            // d.PrintList(d.divdiffs, d.CurrentLength, "Divdiffs");

                            //normalize the divided differences by (-it)^q / q! * beta_coefficient
                            std::complex<double> norm_real = cal_neg_i_pow(q) * std::pow(t, q) * real / factorial(q) ;
                            std::complex<double> norm_imag = cal_neg_i_pow(q+1) * std::pow(t, q) * (-imag) / factorial(q);// (-it)^q / q! * Imag = (-it)^q / q! * (i) * Imag.real = (-i)^q+1 * t^q * -Imag.real / q!
                            perm_amp = perm_amp + (norm_real + norm_imag) * beta_coefficient;
                            
                            // //print perm_amp
                            // std::cout << "Perm amplitude: " << perm_amp << std::endl;

                        }
                        q_amp += perm_amp * energy_coefficient;

                    }
                    //print the q amplitude
                    std::cout << q << " amplitude: " << q_amp << std::endl;
                    
                    //clear divdiff
                    divdiff_clear_up();
                    
                    //add to total amp
                    total_amp += q_amp;
                }

                if(target == 0){
                    total_amp += 1 ;// TODO: make sure this addition is correct
                }

                // // add the diagonal energy of the hamiltonian for excited state
                // if (qubit == 1){
                //     total_amp += target * (-chi);
                // }

                //print the Total amplitude
                std::cout << "Total amplitude: " << total_amp << std::endl;
                transition_amplitudes.push_back(total_amp);


            }
            // //print transition amplitudes
            // std::cout << "Transition amplitudes: " << std::endl;
            // for (const auto& amp : transition_amplitudes) {
            //     std::cout << amp << std::endl;
            // }

            // // save the unnormailzed transition amplitudes
            // reg_trans_amp = transition_amplitudes;
            
            // normalize transition amplitude 
            double norm = 0;
            for (const auto& amp : transition_amplitudes) {
                norm += std::norm(amp);
            }
            norm = sqrt(norm);
            for (auto& amp : transition_amplitudes) {
                amp = amp / norm;
            }

            all_transition_amplitudes.push_back(transition_amplitudes);

            std::cout << std::endl << std::endl;
        }

        //print all transition amplitudes with location in the vector
        std::cout << "All transition amplitudes: " << std::endl;
        for (int i = 0; i < all_transition_amplitudes.size(); i++) {
            std::cout << "Time: " << times[i] << std::endl;
            for (int j = 0; j < all_transition_amplitudes[i].size(); j++) {
                std::cout << "Target: " << j << " Amplitude: " << all_transition_amplitudes[i][j] << std::endl;
            }
        }


        // Assuming all_transition_amplitudes is a std::vector<std::vector<std::complex<double>>>
        nlohmann::json j;
        for (const auto& vec : all_transition_amplitudes) {
            nlohmann::json sub_j;
            for (const auto& c : vec) {
                sub_j.push_back({{"real", c.real()}, {"imag", c.imag()}});
            }
            j.push_back(sub_j);
        }
        if (qubit == 0){
            std::ofstream o("amplitudes_g_chg_amp.json");
            o << j << std::endl;
        }else{
            std::ofstream o("amplitudes_e_chg_amp1.json");
            o << j << std::endl;
        }
    // }

    return 0;
}





