#include <iostream>
#include <vector>
#include <functional>
#include <cmath>
#include <complex>
#include "include\divdiffcomplex.h"
#include "include\permutation.h"
#include "include\Serializer.h"
#include <fstream>
#include "include\json-develop\single_include\nlohmann\json.hpp"
#include <string>



#define  MAX_STEP 14
#define  MAX_TARGET 11


int main(int argc, char* argv[]) {
    //original parameters
    double multiplier = 4;
    float sigma = 48;
    float timestep = 1;
    double chi = -279e-6 * 2 * M_PI; 
    // double detuning = chi *1e-2;

    // initial values
    double amplitude = 0.005;
    int q_max = 8;
    
    if(argc > 1){
        amplitude = std::stod(argv[1]);
        q_max = std::stoi(argv[2]);
    }

    //print args
    std::cout << "Amplitude: " << amplitude << std::endl;
    std::cout << "Q_max: " << q_max << std::endl;


    // //print vars
    // std::cout << "Multiplier: " << multiplier << std::endl;
    // std::cout << "Sigma: " << sigma << std::endl;
    // std::cout << "Timestep: " << timestep << std::endl;
    // std::cout << "Chi: " << chi << std::endl;
    // // std::cout << "Detuning: " << detuning << std::endl;


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
    std::vector<std::vector<std::complex<double>>> all_transition_amplitudes_g;
    std::vector<std::vector<std::complex<double>>> all_transition_amplitudes_e;
    
    //define start state
    int start_state = 0;
    int total_steps = 100;

    // Create a 3D vector to store the results
    std::vector<std::vector<std::vector<
    std::tuple<std::vector<std::vector<int>>, std::vector<std::tuple<std::vector<std::vector<DivdiffElement>>, std::complex<double>, int>>>
    >>> perms_and_coefs;

    
    // Loop over target steps and moves
    for (int target_step = 0; target_step <= 11; ++target_step) {
        // Add a new 2D vector for this target step
        perms_and_coefs.push_back(std::vector<std::vector<
            std::tuple<std::vector<std::vector<int>>, std::vector<std::tuple<std::vector<std::vector<DivdiffElement>>, std::complex<double>, int>>>
        >>());

        for (int q = 0; q <= q_max; ++q) {
            std::string filename = "C:\\Users\\Latzres\\OneDrive - Technion\\RoeeL\\Divided-Differences-for-anharmonic-oscillator\\include\\permutations\\dispersive_perms\\dispersive_perms_" 
                + std::to_string(target_step) + "_" + std::to_string(q) + ".txt";


            // Check if file exists
            std::ifstream infile(filename);

            if (!infile) {
                std::cout << "File " << filename << " does not exist." << std::endl;
                // Generate and save permutations 
                generate_and_save_results(total_steps, target_step, q, start_state, filename);
            } else {

                // Add a new vector for this q
                perms_and_coefs[target_step].push_back(std::vector<
                    std::tuple<std::vector<std::vector<int>>, std::vector<std::tuple<std::vector<std::vector<DivdiffElement>>, std::complex<double>, int>>>
                >());

                // Load results from file
                std::tuple<std::vector<std::vector<int>>, std::vector<std::tuple<std::vector<std::vector<DivdiffElement>>, std::complex<double>, int>>> value = load_results(filename);
                perms_and_coefs[target_step][q].push_back(value);
            }
        }
    }


    // int qubit = 1;
    for (int qubit = 0; qubit <=1; qubit++){

        //loop over all times
        for (int i = 0; i < times.size(); i++) { 
            t = times[i];
            // t = 192;
            //print time
            std::cout << std::endl << "Time: " << t << std::endl;

            //vector of all transition amplitude for a given time
            std::vector<std::complex<double>> transition_amplitudes;
            
            //loop over all target steps
            for (int target = 0; target <= MAX_TARGET; target++) {
                std::complex<double> total_amp = 0;

                //print target
                std::cout << "Target: " << target << std::endl;

                for(int q = 0; q < q_max; q++){
                    //print q
                    std::cout << "Q: " << q << std::endl;

                    std::complex<double> q_amp = 0;
                    std::complex<double> energy_coefficient_amp = cal_neg_i_pow(q) * std::pow(amplitude, q);

                    if (target == 0 && q == 0){
                        total_amp += 1;
                        continue;
                    }
                    
                    auto value = perms_and_coefs[target][q][0];
                    auto permutations = std::get<0>(value);
                    auto coefficients = std::get<1>(value);

                    // std::cout << "Permutations:" << std::endl;
                    // for (const auto& combination : permutations) {
                    //     for (int step : combination) {
                    //         std::cout << step << " ";
                    //     }
                    //     std::cout << std::endl;
                    // }
                    
                    //there are no permutaions - skip
                    if (permutations.size() == 0) {
                        continue;
                    }

                    // std::cout << "Coefficients:" << std::endl;
                    


                    for (const auto& coefficient : coefficients) {
                        auto divdiff = std::get<0>(coefficient);
                        auto energy_coefficient = std::get<1>(coefficient);
                        auto final_state = std::get<2>(coefficient);

                        // std::cout << "Energy coefficient: " << energy_coefficient.real() << " + " << energy_coefficient.imag() << "i" << std::endl;
                        // std::cout << "Final state: " << final_state << std::endl;
                        // std::cout << "Divdiff elements:" << std::endl;
                         
                         //calculate divided differences
                        divdiff_init();

                        divdiffcomplex d(MAX_STEP,500);

                        std::complex<double> sum_final_elements = 0;

                        for (const auto& elements : divdiff) {
                            double beta_coefficient = elements[0].coefficient;
                            // std::cout << "Beta coefficient: " << beta_coefficient << std::endl;
                            d.CurrentLength=0;
                            for (const auto& element : elements) {
                                // std::cout << "(" << element.energy << ", " << element.omega << ", " << element.chi << ", " << element.coefficient << ") ";
                                // std::cout << std::endl;
                                
                                std::complex<double> n;
                                n = (-1j * t * (element.energy * chi * qubit + element.omega * omega  + element.chi * chi));
                                // // print n
                                // std::cout << "n: " << n << std::endl;
                                d.AddElement(n);
                            }
                            std::complex<ExExFloat> z_n = d.divdiffs[d.CurrentLength - 1];
                            double real = z_n.real().get_double();
                            double imag = z_n.imag().get_double();
                            
                            // //print real and imag
                            // std::cout << "Divdiff = " << real << " + " << imag << "i" << std::endl;

                            // //print divdiffs elements
                            // d.PrintList(d.divdiffs, d.CurrentLength, "Divdiffs");

                            //normalize the divided differences by (-it)^q / q! * beta_coefficient
                            std::complex<double> norm_real = cal_neg_i_pow(q) * std::pow(t, q) * real / factorial(q) ;
                            std::complex<double> norm_imag = cal_neg_i_pow(q+1) * std::pow(t, q) * (-imag) / factorial(q);// (-it)^q / q! * Imag = (-it)^q / q! * (i) * Imag.real = (-i)^q+1 * t^q * -Imag.real / q!
                            std::complex<double> final_element = (norm_real + norm_imag) * beta_coefficient;

                            //add the divdiff element to the total sum
                            sum_final_elements += final_element; 
                        }
                        q_amp = q_amp + sum_final_elements * energy_coefficient 
                                        * energy_coefficient_amp; 
                    }
                    /*
                    //ladder_permutations(max_step, target, Q = total moves)
                    auto permutations = ladder_permutations(MAX_STEP, target, q);
                    
                    //there are no permutaions - skip
                    if (permutations.size() == 0) {
                        if (target == 0 && q == 0){
                            q_amp += 1;
                        }
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
                                n = (-1j * t * (divdiff.energy * chi * qubit + divdiff.omega * omega  + divdiff.chi * chi));
                                // // print n
                                // std::cout << "n: " << n << std::endl;
                                d.AddElement(n);
                            }
                            int q = d.CurrentLength - 1;
                            std::complex<ExExFloat> z_n = d.divdiffs[q];

                            double real = d.divdiffs[q].real().get_double();
                            double imag = d.divdiffs[q].imag().get_double();
                            
                            // //print real and imag
                            // std::cout << "Divdiff = " << real << " + " << imag << "i" << std::endl;

                            // //print divdiffs elements
                            // d.PrintList(d.divdiffs, d.CurrentLength, "Divdiffs");

                            //normalize the divided differences by (-it)^q / q! * beta_coefficient
                            std::complex<double> norm_real = cal_neg_i_pow(q) * std::pow(t, q) * real / factorial(q) ;
                            std::complex<double> norm_imag = cal_neg_i_pow(q+1) * std::pow(t, q) * (-imag) / factorial(q);// (-it)^q / q! * Imag = (-it)^q / q! * (i) * Imag.real = (-i)^q+1 * t^q * -Imag.real / q!
                            perm_amp = perm_amp + (norm_real + norm_imag) * beta_coefficient;
                            //TODO: separate calculation for t^q / q! - change type to ExExfloat


                            // //print perm_amp
                            // std::cout << "Perm amplitude: " << perm_amp << std::endl;

                        }
                        // // print energy coefficient
                        // std::cout << "energy coefficient: " << energy_coefficient << std::endl;
                        // //print perm_amp
                        // std::cout << "Perm amplitude: " << perm_amp << std::endl;

                        q_amp += perm_amp * energy_coefficient;

                        //print q_amp
                        std::cout << "q_amp: " << q_amp << std::endl << std::endl;


                    }
                    */

                    //print the q amplitude
                    std::cout << q << " q_amp: " << q_amp << std::endl;
                    
                    //clear divdiff
                    divdiff_clear_up();
                    
                    //add to total amp
                    total_amp += q_amp;
                }

                //print the Total amplitude
                std::cout << "Total amplitude: " << total_amp << std::endl;
                transition_amplitudes.push_back(total_amp);
            }

            //print transition amplitudes
            std::cout << "Transition amplitudes: " << std::endl;
            for (const auto& amp : transition_amplitudes) {
                std::cout << amp << std::endl;
            }

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

            if (qubit == 0){
                all_transition_amplitudes_g.push_back(transition_amplitudes);
            }else{
                all_transition_amplitudes_e.push_back(transition_amplitudes);
            }
            // all_transition_amplitudes.push_back(transition_amplitudes);

            std::cout << std::endl << std::endl;
        }


        if (qubit == 0){
            all_transition_amplitudes = all_transition_amplitudes_g;
        }else{
            all_transition_amplitudes = all_transition_amplitudes_e;
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
            std::string filename = "amplitudes_g_chg_amp" + std::to_string(amplitude) + "_q" + std::to_string(q_max) + ".json";
            std::ofstream o(filename);
            o << j << std::endl;
        }else{
            std::string filename = "amplitudes_e_chg_amp" + std::to_string(amplitude) + "_q" + std::to_string(q_max) + ".json";
            std::ofstream o(filename);
            o << j << std::endl;
        }
    }

    return 0;
}





