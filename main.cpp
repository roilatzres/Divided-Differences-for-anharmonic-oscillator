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

#define  MAX_STEP 14
#define  MAX_TARGET 8



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
    for (int target_step = 0; target_step <= 8; ++target_step) {
        all_coefficients.push_back(vector<vector<std::tuple<vector<DivdiffElement>, double, int, vector<int>>>>());

        for (int q = 0; q <= 8; ++q) {
            std::string filename = "./include/parms/dispersive/coeff/dispersive_coeff_" 
                + std::to_string(target_step) + "_" + std::to_string(q) + ".bin";

            // Check if file exists
            std::ifstream infile(filename);
            if (!infile) {
                std::cout << "File " << filename << " does not exist." << std::endl;
                // generate and save coefficients
                generate_and_save_coefficients(total_steps, target_step, q, filename);

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
        for (int i = 0; i < times.size(); i++) { 
            // t = 192;
            t = times[i];
            
            //print time
            std::cout << std::endl << "Time: " << t << std::endl;

            //vector of all transition amplitude for a given time
            vector<complex<double>> transition_amplitudes;
            
            //loop over all target steps
            for (int target = 0; target <= MAX_TARGET; target++) {
                complex<double> total_amp = 0;

                //print target
                std::cout << "Target: " << target << std::endl;

                for(int q = 0; q < q_max; q++){
                    //print q
                    std::cout << "Q: " << q << std::endl;

                    complex<double> q_amp = 0;
                    complex<double> energy_coefficient_amp = cal_neg_i_pow(q) * std::pow(amplitude, q);

                    if (target == 0 && q == 0){
                        total_amp += 1;
                        continue;
                    }
                    
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
                    if (infile) {
                        std::cout << "File " << q_amp_filename << " already exists." << std::endl;
                        
                        // load q_amp from file
                        q_amp = load_complex_from_file(q_amp_filename);

                         q_amp = q_amp * energy_coefficient_amp;
                        
                        //print the q amplitude
                        std::cout << q << " q_amp: " << q_amp << std::endl;
                        
                        //add to total amp
                        total_amp += q_amp;
                        continue;
                    }

                    //print starting divdiff
                    std::cout << "Divdiff:" << std::endl;
                    
                    
                    // Need to calculate q_amp
                    divdiff_init();
                    for (const auto& coefficient : coefficients) {
                        
                        q_amp += cal_divdiff(coefficient, t, q, qubit, chi, omega);
                    }
                    
                    //print  finishing divdiff
                    std::cout << "Finishing Divdiff:" << std::endl;


                    //clear divdiff
                    divdiff_clear_up();

                    // save q_amp to a file
                   save_complex_to_file(q_amp, q_amp_filename);

                    q_amp = q_amp * energy_coefficient_amp;

                    //print the q amplitude
                    std::cout << q << " q_amp: " << q_amp << std::endl;
                    
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





