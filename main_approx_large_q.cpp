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

using std::vector;
using std::complex;

// #define  MAX_STEP 14
#define  MAX_TARGET 6
#define MAX_Q 13



int main(int argc, char* argv[]) {
    complex<double> i_num(0, 1);
    //print i_num
    std::cout << "i_num: " << i_num << std::endl;
    //original parameters
    double multiplier = 4;
    float sigma = 48;
    float timestep = 1;
    double chi = -279e-6 * 2 * M_PI; 
    // double detuning = chi *1e-2;

    // initial values
    double amplitude = 0.001;
    int q_max = 10; // TODO: try 11
    
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
    vector<double> times;
    for (double t = timestep; t <= tfinal; t += timestep) {//TODO: check if t=1 is correct (t=0 gives 0 amplitude)
        times.push_back(t);
        //print times
        // std::cout << "Time: " << t << std::endl;
    }


    //calculate the time evolution
    double t;

    // vector of all transition ampli tudes for all times: [time][target_state][amp]
    vector<vector<complex<double>>> all_transition_amplitudes;
    vector<vector<complex<double>>> all_transition_amplitudes_g;
    vector<vector<complex<double>>> all_transition_amplitudes_e;
    
    //define start state
    int start_state = 0;
    int total_steps = 100;

    // Create a 3D vector to store the results
    vector< //target_step
    vector< //q
    vector< //permutations
    std::tuple<vector<DivdiffElement>, double, int, vector<int>>> // coefficients
    >> all_coefficients;

    
    // Loop over target steps and moves
    for (int target_step = 0; target_step <= MAX_TARGET; ++target_step) {
        all_coefficients.push_back(vector<vector<std::tuple<vector<DivdiffElement>, double, int, vector<int>>>>());

        for (int q = 0; q <= MAX_Q; ++q) {
            std::string filename = "./include/parms/dispersive/coeff/dispersive_coeff_" 
                + std::to_string(target_step) + "_" + std::to_string(q) + ".bin";

            // Check if file exists
            std::ifstream infile(filename);
            if (!infile) {
                std::cout << "File " << filename << " does not exist." << std::endl;
                // generate and save coefficients
                generate_and_save_coefficients(start_state, total_steps, target_step, q, filename);

            } else {

                // Add a new vector for this q
                all_coefficients[target_step].push_back(vector<std::tuple<vector<DivdiffElement>, double, int, vector<int>>>());
                
                // Load all_coefficients from file
                vector<std::tuple<vector<DivdiffElement>, double, int, vector<int>>> value = load_coefficients(filename);
                all_coefficients[target_step][q] = value;
                                
            }
        }
    }


    // int qubit = 1;
    for (int qubit = 0; qubit <=1; qubit++){

        //loop over all times
        // for (int i = 0; i < times.size(); i++) { 
            // t = times[i];
            t = 192;
            
            //print time
            std::cout << std::endl << "Main-Time: " << t << std::endl;

            //vector of all transition amplitude for a given time
            vector<complex<double>> transition_amplitudes;
            
            //loop over all target steps
            for (int target = 0; target <= MAX_TARGET; target++) {
                complex<double> total_amp = 0;

                //print target
                std::cout << "Target: " << target << std::endl;

                for(int q = 0; q <= q_max; q++){
                    //print q
                    // std::cout << "Q: " << q << std::endl;

                    complex<double> q_amp = 0; //TODO: change to EXEXFLOAT
                    complex<double> energy_coefficient_amp = cal_neg_i_pow(q) * std::pow(amplitude, q);

                    if (target == 0 && q == 0){
                        total_amp += 1;
                        continue;
                    }

                    if(q > 11){
                        if ((q + target) % 2 != 0){ // q + target is odd - skip TODO: check for asymetric case
                            continue;
                        }
                        
                        int mid_q = (q + target)/2;
                        int extra_steps = 1; //amount of extra steps we allow for the permutaions
                        vector<vector<int>> large_permutations;
                        
                        for(int i = 0; i <= extra_steps; i++){
                            vector<vector<int>> permutations_up = ladder_permutations(0, mid_q, mid_q, i + mid_q);
                            if(permutations_up.size() == 0){
                                continue;
                            }
                            // //print i 
                            // std::cout << "i: " << i << std::endl;
                            for(int j = 0; j <= extra_steps; j++){
                                vector<vector<int>> permutations_down = ladder_permutations(mid_q, mid_q, target, j + q - mid_q);
                                if(permutations_down.size() == 0){
                                    continue;
                                }
                                // //print j
                                // std::cout << "j: " << j << std::endl;

                                //combine permutations
                                for (const auto& up : permutations_up) {
                                    for (const auto& down : permutations_down) {
                                        vector<int> combined = up;
                                        combined.insert(combined.end(), down.begin(), down.end());
                                        large_permutations.push_back(combined);
                                    }
                                }
                                
                            }

                        }
                        std::cout << "num perm " << large_permutations.size() << std::endl;

                        // for (const auto& combination : large_permutations) {
                        for (int j = 0; j < large_permutations.size(); j++) {
                            const auto& combination = large_permutations[j];

                            double total_energy = 0;
                            int current_state = start_state;

                            // initialize vector for all states in the permutation
                            vector<int> states;
                            states.push_back(start_state);

                            // Iterate through each step in the permutation
                            for (int i = 0; i < combination.size(); ++i) {
                                // Update current state
                                int new_state = current_state + combination[i];
                                states.push_back(new_state);
                                
                                // Update energy_coefficient_real
                                if (combination[i] == 1) {
                                     total_energy += sqrt(current_state + 1); 
                                } else {//combination[i] == -1
                                     total_energy += sqrt(current_state);
                                }
                                current_state = new_state;
                            }
                            // //print size combination
                            // std::cout << "combinatin length: " << combination.size() << std::endl;
                            // // print q 
                            // std::cout << "q: " << q << std::endl;
                            // //prit states
                            // std::cout << "States: " << std::endl;
                            // for (int state : states) {
                            //     std::cout << state << " ";
                            // }
                            // std::cout << std::endl;

                            // //print combination
                            // std::cout << "Combination: " << std::endl;
                            // for (int step : combination) {
                            //     std::cout << step << " ";
                            // }
                            // std::cout << std::endl;
                            // std::cout << std::endl;

                            // // print total energy
                            // std::cout << "Total energy: " << total_energy << std::endl;

                            int effective_q = combination.size();  
                            
                            //total contribution to q_amp : gamma_q *
                            //                              total_energy *
                            //                              (-2t)^q / q! * 
                            //                              exp(-it <E>) * 
                            //                              (cos(omega * k* t / q+1) [where k goes from 1 to q])
                            ExExFloat effective_energy_coefficient_amp_real =  ExPow(amplitude, effective_q);
                            // std::cout << "Effective energy coefficient real: " << std::endl;
                            // effective_energy_coefficient_amp_real.print();
                            // std::cout << std::endl;
                            // // effective_energy_coefficient_amp_real.imag().print();
                            // // std::cout << std::endl;
                            complex<ExExFloat> effective_energy_coefficient_amp_i = cal_neg_i_pow(effective_q);
                            // std::cout << "Effective energy coefficient imag: " << std::endl;
                            // effective_energy_coefficient_amp_i.real().print();
                            // std::cout << std::endl;
                            // effective_energy_coefficient_amp_i.imag().print();
                            // std::cout << std::endl;

                            complex<ExExFloat> effective_energy_coefficient_amp = effective_energy_coefficient_amp_real * effective_energy_coefficient_amp_i;

                            //print effective energy coefficient
                            std::cout << "Effective energy coefficient: " << std::endl;
                            effective_energy_coefficient_amp.real().print();
                            std::cout << std::endl;
                            effective_energy_coefficient_amp.imag().print();
                            std::cout << std::endl;

                            //cal cos part
                            ExExFloat cos_part = 1;
                            for (int k = 1; k <= effective_q; k++) {
                                cos_part = cos_part * cos(omega * k * t / (effective_q + 1));
                            }
                            //print cos part
                            std::cout << "Cos part: " << std::endl;
                            cos_part.print();
                            std::cout << std::endl;

                            //cal total energy
                            complex<double> exp_power = total_energy * t / (effective_q + 1);
                            // std::cout << "Exp power: " << exp_power << std::endl;
                            complex<double> total_energy_exp = exp(exp_power * i_num);
                            //print total energy exp
                            std::cout << "Total energy exp: " << total_energy_exp << std::endl;
                            

                            divdiff_init(); // divdiff init is required for ExFactorial
                            ExExFloat e_q = effective_q;

                            ExExFloat time_coef_real = ExPow(2*t, effective_q) / ExFactorial(e_q, effective_q);                            
                            complex<ExExFloat> time_coef_i = cal_neg_i_pow(effective_q);
                            complex<ExExFloat> time_coef = time_coef_real * time_coef_i;
                            // //print time coef
                            // std::cout << "Time Coef: " << std::endl;
                            // time_coef.real().print();
                            // std::cout << std::endl;
                            // time_coef.imag().print();
                            // std::cout << std::endl;

                            std::cout << "Time Coef: " << time_coef.real().get_double() << 
                            "+ i" << time_coef.imag().get_double() << std::endl;

                            //print approx part of the sum of divdiff
                            complex<ExExFloat> approx_part = time_coef * total_energy_exp * cos_part ; // approx part of the sum of divdiff
                            // std::cout << "Approx part: " << std::endl;
                            // approx_part.real().print();
                            // std::cout << std::endl;
                            // approx_part.imag().print();
                            // std::cout << std::endl;
                            std::cout << "Approx part: " << approx_part.real().get_double() <<
                            "+ i" << approx_part.imag().get_double() << std::endl;

                            complex<ExExFloat> total_effective_amp = effective_energy_coefficient_amp * total_energy * //original part
                                                                     time_coef * total_energy_exp * cos_part ; // approx part of the sum of divdiff
                            // //print total effective amp
                            // std::cout << "Total effective amp: " << std::endl;
                            // total_effective_amp.real().print();
                            // std::cout << std::endl;
                            // total_effective_amp.imag().print();
                            // std::cout << std::endl;
                            // std::cout << "Total effective amp: " << total_effective_amp.real().get_double() << 
                            // "+ i" << total_effective_amp.imag().get_double() << std::endl;

                            //add to total amp
                            q_amp += total_effective_amp.real().get_double() + total_effective_amp.imag().get_double() * i_num;

                            divdiff_clear_up();
                            
                            
                        }
                        //print the q amplitude
                        std::cout << q << " new  q_amp: " << q_amp << std::endl;

                        //load original q_amp
                        std::string q_amp_filename = "./include/parms/dispersive/q_amp_recursive/dispersive_q_amp_qubit_" 
                            + std::to_string(qubit) + "_t_" + std::to_string((int)t) + "_target_" + std::to_string(target) + "_q" + std::to_string(q) + ".bin";
                        std::ifstream infile(q_amp_filename);
                        if (infile) {
                            // std::cout << "File " << q_amp_filename << " already exists." << std::endl;
                            
                            // load q_amp from file
                            complex<double> q_amp_orig = load_complex_from_file(q_amp_filename);

                            q_amp_orig = q_amp_orig * energy_coefficient_amp;
                            
                            // print the q amplitude
                            std::cout << q << " Original q_amp_orig: " << q_amp_orig << std::endl;
                        }
                        
                    // }
                    }else{//q<=11 
                        auto coefficients = all_coefficients[target][q];

                        // std::cout << "Permutations:" << std::endl;
                        // for (const auto& combination : permutations) {
                        //     for (int step : combination) {
                        //         std::cout << step << " ";
                        //     }
                        //     std::cout << std::endl;
                        // }
                        
                        //there are no permutaions - skip
                        if (coefficients.size() == 0) {
                            continue;
                        }

                        
                        // check if file for q_amp already exists
                        std::string q_amp_filename = "./include/parms/dispersive/q_amp_recursive/dispersive_q_amp_qubit_" 
                            + std::to_string(qubit) + "_t_" + std::to_string((int)t) + "_target_" + std::to_string(target) + "_q" + std::to_string(q) + ".bin";
                        std::ifstream infile(q_amp_filename);
                        if (infile) {//TODO: remove q != 12
                            // std::cout << "File " << q_amp_filename << " already exists." << std::endl;
                            
                            // load q_amp from file
                            q_amp = load_complex_from_file(q_amp_filename);

                            q_amp = q_amp * energy_coefficient_amp;
                            
                            // if(q > 12){ *************************************
                                // print the q amplitude
                                std::cout << q << " Original q_amp: " << q_amp << std::endl;
                            // }***********************************************
                            //add to total amp
                            total_amp += q_amp; 
                            continue;
                        }


                        // Need to calculate q_amp
                        divdiff_init();
                        for (const auto& coefficient : coefficients) {
                            
                            q_amp += cal_divdiff(coefficient, t, q, qubit, chi, omega, target);
                        }
                        //clear divdiff
                        divdiff_clear_up();


                        // save q_amp to a file
                        save_complex_to_file(q_amp, q_amp_filename);

                        q_amp = q_amp * energy_coefficient_amp;
                    // //else restore location********************
                    }

                    // // print the q amplitude
                    // std::cout << q << " q_amp: " << q_amp << std::endl;
                    if(q > 11){
                        // print the q amplitude
                        std::cout << q << " Original q_amp: " << q_amp << std::endl;
                    }
                    
                    //add to total amp
                    total_amp += q_amp;
                }

                // //print the Total amplitude
                // std::cout << "Total amplitude: " << total_amp << std::endl << std::endl;
                transition_amplitudes.push_back(total_amp);
            }

            // //print transition amplitudes
            std::cout << "Transition amplitudes: " << std::endl;
            int target_nums = 0;
            for (const auto& amp : transition_amplitudes) {
                std::cout << target_nums << ": " << amp << std::endl;
                target_nums++;
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
        
        // //end of time loop
        // }


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

        // Assuming all_transition_amplitudes is a vector<vector<complex<double>>>
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





