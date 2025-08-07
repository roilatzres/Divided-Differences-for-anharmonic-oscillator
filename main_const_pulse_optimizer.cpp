#include <iostream>
#include <vector>
#include <functional>
#include <cmath>
#include <complex>
#include "include/divdiffcomplex.h"
#include "include/permutation.h"
// #include "include/optimizer_helper.h"
#include "include/Serializer.h"
#include <fstream>
#include "include/json-develop/single_include/nlohmann/json.hpp"
#include <string>
// #include "main_orig.h"
// #include <dlib/matrix.h>
#include "include/dlib/dlib/optimization.h"
#include "include/dlib/dlib/matrix.h"


using std::vector;
using std::complex;

#define  MAX_STEP 300
// #define  MAX_TARGET 15
#define MAX_Q 30
#define M_PI 3.14159265358979323846


//function to print complex_ex
void print_complex_Ex(complex_Ex a){
    a.real.print();
    std::cout << "+ ";
    a.imag.print();
    std::cout << "i" << std::endl; 
    std::cout << std::endl;
}


vector<vector<vector<vector<complex_Ex>>>> coef_matrices(double time, int q_max, int max_target, int qubit){
    complex<double> i_num(0, 1);
    
    //original parameters
    double multiplier = 4;
    float sigma = 48;
    float timestep = 1;
    double chi = -279e-6 * 2 * M_PI; 
    // double detuning = chi *1e-2;


    //TODO: notice that amp is double of the value to put in antiSymSim for sin pulse
    
    //print args
    std::cout << "q_max: " << q_max << std::endl;
    std::cout << "max target: " << max_target << std::endl;

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
    for (double t = timestep; t <= tfinal; t += timestep) {
        times.push_back(t);
        //print times
        // std::cout << "Time: " << t << std::endl;
    }

    // //calculate the time evolution
    double t;

    // all coef matrix - [time][power of q][start_state][target_state] 
    vector<vector<vector<vector<complex_Ex>>>> all_coefs = vector<vector<vector<vector<complex_Ex>>>>(time, vector<vector<vector<complex_Ex>>>(q_max + 1, vector<vector<complex_Ex>>(max_target, vector<complex_Ex>(max_target)))); 


    for (int i = 0; i < time; i++) {
            t = times[i];
        //     // t = 96;
            
        //print time
        std::cout << std::endl << "Main-Time: " << t << std::endl;

        
        //loop over all target steps
        for(int q = 0; q <= q_max; q++){
            cout << "Q: " << q << std::endl;
            
            for(int start_state = 0; start_state < max_target; start_state++){
                cout << "Start state: " << start_state << std::endl;
                
                auto transition_amplitudes_start_state = vector<complex_Ex> {};
                for (int target = 0; target < max_target; target++) {
                    //print target
                    std::cout << "Target: " << target << std::endl;
                    
                    complex_Ex q_coef;
                    
                    // easist case - if target is the same as start state and q is 0
                    if (target == start_state && q == 0){
                        // print case 1
                        std::cout << "Case 1: target is the same as start state and q is 0" << std::endl;
                        complex<double> exponent = std::exp(-i_num * t * (start_state * chi * qubit));
                        q_coef.real = exponent.real();
                        q_coef.imag = exponent.imag();
                        all_coefs[i][q][start_state][target] = q_coef;
                        cout << "q_coef" << std::endl;
                        print_complex_Ex(q_coef);
                        continue;
                    }

                    // second case - we already have the coefficient calculated
                    std::string coeff_filename = "./include/q_amp/optimizer/dispersive_coeff_start_state_" + std::to_string(start_state) 
                        + "_qubit_" + std::to_string(qubit) + "_t_" + std::to_string((int)t) + "_target_" + std::to_string(target) + 
                        "_q" + std::to_string(q) + "_time" + std::to_string(t) + ".bin";
                    std::ifstream infile_q_amp(coeff_filename);

                    if (infile_q_amp) {
                        std::cout << "Case 2: coefficient already exists" << std::endl;
                        q_coef = load_binary_complexEx(coeff_filename);
                        all_coefs[i][q][start_state][target] = q_coef;
                        print_complex_Ex(q_coef);
                        continue;
                    }

                    // third case - we need to calculate the coefficient
                    std::cout << "Case 3: coefficient needs to be calculated" << std::endl;
                    auto permutations = vector<vector<int>> {};
                    std::string perm_filename = "./include/permutations/dispersive_const_perms/start_state" +
                    std::to_string(start_state) + "_target_" + std::to_string(target) + "_q" + std::to_string(q) + ".txt";
                    
                    std::ifstream infile_perm(perm_filename);
                    if (infile_perm) {
                        permutations = load_binary_perm(perm_filename);
                    } else {
                        permutations = ladder_permutations(start_state, MAX_STEP, target, q);
                        if (permutations.empty()) {
                            std::cout << "No permutations found for start state " << start_state << " and target " << target << " with q = " << q << std::endl;
                            q_coef.real = 0;
                            q_coef.imag = 0;
                            all_coefs[i][q][start_state][target] = q_coef;
                            continue;
                        } else {
                            save_binary_perm(permutations, perm_filename);
                        }
                    }
                    
                    vector<std::tuple<double, int, vector<int>>> all_coefficients_const_pulse;
                    std::string all_coef_filename = "./include/parms/dispersive_const/start_state" +
                    std::to_string(start_state) + "_target_" + std::to_string(target) + "_q" + std::to_string(q) + ".txt";
                
                    std::ifstream infile_coef(all_coef_filename);
                    if (infile_coef) {
                        all_coefficients_const_pulse = load_binary_coef(all_coef_filename);
                    } else {
                        for (const auto& permutation : permutations) {
                            auto coefficients = cal_coefficient_const_pulse(permutation, start_state);
                            all_coefficients_const_pulse.push_back(coefficients);
                        }
                        save_binary_coef(all_coefficients_const_pulse, all_coef_filename);
                    }

                        
                    divdiff_init();
                    for (const auto& coefficient : all_coefficients_const_pulse) {
                        complex_Ex res = cal_divdiff_const_amp(coefficient, t, q, qubit, chi);
                        q_coef.real += res.real;
                        q_coef.imag += res.imag;
                    }
                    divdiff_clear_up();
                    
                    save_binary_complexEx(q_coef, coeff_filename);
                    
                    all_coefs[i][q][start_state][target] = q_coef;
                    cout << "q_coef" << std::endl;
                    print_complex_Ex(q_coef);

                }
            }

        }
        
    //end of time loop
    }

    return all_coefs;
}

// // int main(int argc, char* argv[]) {
// vector<vector<vector<complex_Ex>>> solve_for_t(double time, int q_max, double amplitude, int start_state, int max_target, int qubit){
//     complex<double> i_num(0, 1);
    
//     //original parameters
//     double multiplier = 4;
//     float sigma = 48;
//     float timestep = 1;
//     double chi = -279e-6 * 2 * M_PI; 
//     // double detuning = chi *1e-2;


//     //TODO: notice that amp is double of the value to put in antiSymSim for sin pulse
    
//     //print args
//     std::cout << "Amplitude: " << amplitude << std::endl;
//     std::cout << "q_max: " << q_max << std::endl;

//     // //print vars
//     // std::cout << "Multiplier: " << multiplier << std::endl;
//     // std::cout << "Sigma: " << sigma << std::endl;
//     // std::cout << "Timestep: " << timestep << std::endl;
//     // std::cout << "Chi: " << chi << std::endl;
//     // // std::cout << "Detuning: " << detuning << std::endl;


//     double tfinal = multiplier * sigma;
//     double omega = (2 * M_PI) / tfinal;

//     std::cout << "tfinal: " << tfinal << std::endl;
//     std::cout << "Omega: " << omega << std::endl;
//     vector<double> times;
//     for (double t = timestep; t <= tfinal; t += timestep) {
//         times.push_back(t);
//         //print times
//         // std::cout << "Time: " << t << std::endl;
//     }

//     // //calculate the time evolution
//     double t;

//     // vector of all transition amplitudes for all times: [time][target_state][amp]
//     vector<vector<vector<complex_Ex>>> all_transition_amplitudes(2);
//     vector<vector<complex<double>>> all_transition_amplitudes_g;
//     vector<vector<complex<double>>> all_transition_amplitudes_e;
    
//     //define start state
//     // int start_state = 0;
//     int total_steps = 100;

//     // Create a 3D vector to store the results
//     vector< //target_step
//     vector< //q
//     vector< //permutations
//     std::tuple<vector<DivdiffElement>, double, int, vector<int>>> // coefficients
//     >> all_coefficients;
    
//     // load_all_coefficients(all_coefficients, start_state, total_steps, MAX_TARGET, q_max);
//     // for (int qubit = 0; qubit < 2; qubit++){//TODO: qubit needs to be passed as an argument

//         // //loop over all times
//         // // for (int i = 0; i < times.size(); i++) { 
//         for (int i = 0; i < time; i++) {
//             t = times[i];
//         //     // t = 96;
            
//             //print time
//             std::cout << std::endl << "Main-Time: " << t << std::endl;

//             //vector of all transition amplitude for a given time
//             auto transition_amplitudes = vector<complex_Ex> {};
            
//             //loop over all target steps
//             for (int target = 0; target < max_target; target++) {
//                 complex_Ex total_amp;
//                 total_amp.real = 0;
//                 total_amp.imag = 0;

//                 // //print total_amp
//                 // std::cout <<  "init Total amp: "; 
//                 // print_complex_Ex(total_amp);

//                 //print target
//                 std::cout << "Target: " << target << std::endl;

//                 for(int q = 0; q <= q_max; q++){
//                     complex_Ex q_amp; 
//                     q_amp.real = 0;
//                     q_amp.imag = 0;

//                     double energy_coefficient_amp = std::pow(amplitude, q);

//                     if (target == start_state && q == 0){
//                         complex<double> exponent = std::exp(-i_num * t * (start_state * chi * qubit));
//                         total_amp.real = exponent.real();
//                         total_amp.imag = exponent.imag();
//                         continue;
//                     }

//                     std::string q_amp_filename = "./include/q_amp/dispersive_const_new/dispersive_q_amp_start_state_" + std::to_string(start_state) 
//                          + "_qubit_" + std::to_string(qubit) + "_t_" + std::to_string((int)t) + "_target_" + std::to_string(target) + "_q" + std::to_string(q) + ".bin";
//                     std::ifstream infile_q_amp(q_amp_filename);

//                     if (infile_q_amp) {
//                         q_amp = load_binary_complexEx(q_amp_filename);
//                         q_amp.real *= energy_coefficient_amp; 
//                         q_amp.imag *= energy_coefficient_amp; 

//                         divdiff_init();
//                         total_amp.real += q_amp.real.get_double();
//                         total_amp.imag += q_amp.imag.get_double();
//                         divdiff_clear_up();
//                         continue;
//                     }

//                     auto permutations = vector<vector<int>> {};
//                     std::string perm_filename = "./include/permutations/dispersive_const_perms/start_state" +
//                     std::to_string(start_state) + "_target_" + std::to_string(target) + "_q" + std::to_string(q) + ".txt";
                    
//                     std::ifstream infile_perm(perm_filename);
//                     if (infile_perm) {
//                         permutations = load_binary_perm(perm_filename);
//                     } else {
//                         permutations = ladder_permutations(start_state, MAX_STEP, target, q);
//                         if (permutations.empty()) {
//                             continue;
//                         } else {
//                             save_binary_perm(permutations, perm_filename);
//                         }
//                     }
                       
//                     vector<std::tuple<double, int, vector<int>>> all_coefficients_const_pulse;
//                     std::string all_coef_filename = "./include/parms/dispersive_const/start_state" +
//                     std::to_string(start_state) + "_target_" + std::to_string(target) + "_q" + std::to_string(q) + ".txt";
                   
//                     std::ifstream infile_coef(all_coef_filename);
//                     if (infile_coef) {
//                         all_coefficients_const_pulse = load_binary_coef(all_coef_filename);
//                     } else {
//                         for (const auto& permutation : permutations) {
//                             auto coefficients = cal_coefficient_const_pulse(permutation, start_state);
//                             all_coefficients_const_pulse.push_back(coefficients);
//                         }
//                         save_binary_coef(all_coefficients_const_pulse, all_coef_filename);
//                     }

                    
//                     divdiff_init();
//                     for (const auto& coefficient : all_coefficients_const_pulse) {
//                         complex_Ex res = cal_divdiff_const_amp(coefficient, t, q, qubit, chi);
//                         q_amp.real += res.real;
//                         q_amp.imag += res.imag;
//                     }
//                     divdiff_clear_up();
                    
//                     save_binary_complexEx(q_amp, q_amp_filename);
                    
//                     q_amp.real *= energy_coefficient_amp; 
//                     q_amp.imag *= energy_coefficient_amp; 

//                     divdiff_init();
//                     total_amp.real += q_amp.real.get_double();
//                     total_amp.imag += q_amp.imag.get_double();
//                     divdiff_clear_up();
//                 }

//                 // //print the Total amplitude
//                 // // std::cout << "Total amplitude: " << total_amp << std::endl << std::endl;
//                 // std::cout <<  "Total amp: "; 
//                 // print_complex_Ex(total_amp);
//                 transition_amplitudes.push_back(total_amp);
//             }

//             // // //print transition amplitudes
//             // std::cout << "Transition amplitudes: " << std::endl;
//             // int target_nums = 0;
//             // for (auto& amp : transition_amplitudes) {
//             //     std::cout << target_nums << ": "; 
//             //     print_complex_Ex(amp);
//             //     target_nums++;
//             // }

//             // normalize transition amplitude 
//             double norm = 0;
//             ExExFloat Ex_norm;
//             Ex_norm = 0;
//             divdiff_init();
//             for (auto& amp : transition_amplitudes) {
//                 Ex_norm += amp.real * amp.real + amp.imag * amp.imag;
//             }
//             norm = sqrt(Ex_norm.get_double());
//             for (auto& amp : transition_amplitudes) {
//                 amp.real = amp.real / norm;
//                 amp.imag = amp.imag / norm;
//             }
//             divdiff_clear_up();
            
//             all_transition_amplitudes[qubit].push_back(transition_amplitudes);
            
//         //end of time loop
//         }

//         // //print all transition amplitudes with location in the vector
//         // std::cout << "All transition amplitudes: " << std::endl;
//         // for (int i = 0; i < all_transition_amplitudes[qubit].size(); i++) {
//         //     std::cout << "Time: " << times[i] << std::endl;
//         //     for (int j = 0; j < all_transition_amplitudes[qubit][i].size(); j++) {
//         //         std::cout << "Target: " << j << " Amplitude: ";
//         //         print_complex_Ex(all_transition_amplitudes[qubit][i][j]);
//         //     }
//         // }
    
//     // } //end of qubit loop

//     std::cout << "Exiting solve_for_t" << std::endl;

//     return all_transition_amplitudes;
// }






int main(int argc, char* argv[]) {
    //get the parameters
    // initial values
    double orig_amplitude = 0.002;
    int q_max = 10; 
    int max_target = 8;

    int num_pulses = 2;

    int qubit = 0; // default qubit to 0
    // std::string final_state_filename = "amplitudes_g_chg_pulse3amp0.002000_q10.json"; 

    
    if(argc > 1){
        q_max = std::stoi(argv[1]);
        max_target = std::stoi(argv[2]);
        num_pulses = std::stoi(argv[3]);
        // qubit = std::stoi(argv[5]); // get the qubit number from command line arguments
        // if (argc > 6) {
        //     final_state_filename = argv[6]; // get the final state filename from command line arguments
        // }
    }


    double final_t = 192;
    double sample_t = final_t / 2 / num_pulses;
    int init_start_state = 0;
    

    // create start state
    vector<complex_Ex> start_state_orig = vector<complex_Ex>(max_target);
    for (int i = 0; i < max_target; i++) {
        complex_Ex c;
        if(i == 0){
            c.real = 1; // set the first state to 1
            c.imag = 0;
        } else {
            c.real = 0; // set the rest to 0
            c.imag = 0;
        }
        start_state_orig[i] = c;
    }
    

    for(int qubit = 0; qubit < 2; qubit++){
     
        // cal coef_matrices
        std::cout << "Calculating coefficient matrices..." << std::endl;
        vector<vector<vector<vector<complex_Ex>>>> all_coefs = coef_matrices(sample_t, q_max, max_target, qubit);
        std::cout << "Coefficient matrices calculated." << std::endl;
        
        // save all_coefs with values to json file
    std::string all_coefs_filename = (qubit == 0) ? "all_coefs_g_q" + std::to_string(q_max) + "_n" + std::to_string(max_target) + ".json" 
        : "all_coefs_e_q" + std::to_string(q_max) + "_n" + std::to_string(max_target) + ".json";
    nlohmann::json j_all_coefs;
    for (auto& time_vec : all_coefs) {
        nlohmann::json j_time;
        for (auto& q_vec : time_vec) {
            nlohmann::json j_q;
            for (auto& start_state_vec : q_vec) {
                nlohmann::json j_start_state;
                for (auto& target_state : start_state_vec) {
                    j_start_state.push_back({{"real", target_state.real.get_double()}, {"imag", target_state.imag.get_double()}});
                }
                j_q.push_back(j_start_state);
            }
            j_time.push_back(j_q);
        }
        j_all_coefs.push_back(j_time);
    }
    std::ofstream out_all_coefs(all_coefs_filename);
    if (out_all_coefs.is_open()) {
        out_all_coefs << j_all_coefs.dump(4); // Pretty print with an indent of 4 spaces
        out_all_coefs.close();
        std::cout << "Coefficient matrices saved to " << all_coefs_filename << std::endl;
    } else {
        std::cerr << "Error opening file for saving coefficient matrices: " << all_coefs_filename << std::endl;
        return 1; // Exit if the file cannot be opened
    }
    
    //print all_coefs 
    std::cout << "All coefficients: " << std::endl;
    for ( auto& time_vec : all_coefs) {
        std::cout << "Time step: " << &time_vec - &all_coefs[0] << std::endl; // Print the index of the time step
        for ( auto& q_vec : time_vec) {
            cout << endl;
            std::cout << "Q: " << &q_vec - &time_vec[0] << std::endl; // Print the index of q
            for ( auto& start_state_vec : q_vec) {
                std::cout << "Start state: " << &start_state_vec - &q_vec[0] << std::endl; // Print the index of start state
                for ( auto& target_state : start_state_vec) {
                    std::cout << "Target state: " << &target_state - &start_state_vec[0] << std::endl; // Print the index of target state
                    std::cout << "Real: " << target_state.real.get_double() << ", Imag: " << target_state.imag.get_double() << std::endl;
                    cout << endl;
                }
            }
        }
    }
    std::cout << "All coefficients saved to " << all_coefs_filename << std::endl;
}
    
    

    return 0;
}
