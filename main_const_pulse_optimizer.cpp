#include <iostream>
#include <vector>
#include <functional>
#include <cmath>
#include <complex>
#include "include/divdiffcomplex.h"
#include "include/permutation.h"
#include "include/Serializer.h"
#include <fstream>
#include "include/json-develop/single_include/nlohmann/json.hpp"
#include <string>
// #include "main_orig.h"

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


// int main(int argc, char* argv[]) {
vector<vector<vector<complex_Ex>>> solve_for_t(double time, int q_max, double amplitude, int start_state, int max_target, int qubit){
    complex<double> i_num(0, 1);
    
    //original parameters
    double multiplier = 4;
    float sigma = 48;
    float timestep = 1;
    double chi = -279e-6 * 2 * M_PI; 
    // double detuning = chi *1e-2;


    //TODO: notice that amp is double of the value to put in antiSymSim for sin pulse
    
    //print args
    std::cout << "Amplitude: " << amplitude << std::endl;
    std::cout << "q_max: " << q_max << std::endl;

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

    // vector of all transition amplitudes for all times: [time][target_state][amp]
    vector<vector<vector<complex_Ex>>> all_transition_amplitudes(2);
    vector<vector<complex<double>>> all_transition_amplitudes_g;
    vector<vector<complex<double>>> all_transition_amplitudes_e;
    
    //define start state
    // int start_state = 0;
    int total_steps = 100;

    // Create a 3D vector to store the results
    vector< //target_step
    vector< //q
    vector< //permutations
    std::tuple<vector<DivdiffElement>, double, int, vector<int>>> // coefficients
    >> all_coefficients;
    
    // load_all_coefficients(all_coefficients, start_state, total_steps, MAX_TARGET, q_max);
    // for (int qubit = 0; qubit < 2; qubit++){//TODO: qubit needs to be passed as an argument

        // //loop over all times
        // // for (int i = 0; i < times.size(); i++) { //TODO: i=0
        for (int i = 0; i < time; i++) { //TODO: i=0
            t = times[i];
        //     // t = 96;
            
            //print time
            std::cout << std::endl << "Main-Time: " << t << std::endl;

            //vector of all transition amplitude for a given time
            auto transition_amplitudes = vector<complex_Ex> {};
            
            //loop over all target steps
            for (int target = 0; target < max_target; target++) {
                complex_Ex total_amp;
                total_amp.real = 0;
                total_amp.imag = 0;

                // //print total_amp
                // std::cout <<  "init Total amp: "; 
                // print_complex_Ex(total_amp);

                //print target
                std::cout << "Target: " << target << std::endl;

                for(int q = 0; q <= q_max; q++){
                    complex_Ex q_amp; 
                    q_amp.real = 0;
                    q_amp.imag = 0;

                    double energy_coefficient_amp = std::pow(amplitude, q);

                    if (target == start_state && q == 0){
                        complex<double> exponent = std::exp(-i_num * t * (start_state * chi * qubit));
                        total_amp.real = exponent.real();
                        total_amp.imag = exponent.imag();
                        continue;
                    }

                    std::string q_amp_filename = "./include/q_amp/dispersive_const_new/dispersive_q_amp_start_state_" + std::to_string(start_state) 
                         + "_qubit_" + std::to_string(qubit) + "_t_" + std::to_string((int)t) + "_target_" + std::to_string(target) + "_q" + std::to_string(q) + ".bin";
                    std::ifstream infile_q_amp(q_amp_filename);

                    if (infile_q_amp) {
                        q_amp = load_binary_complexEx(q_amp_filename);
                        q_amp.real *= energy_coefficient_amp; 
                        q_amp.imag *= energy_coefficient_amp; 

                        divdiff_init();
                        total_amp.real += q_amp.real.get_double();
                        total_amp.imag += q_amp.imag.get_double();
                        divdiff_clear_up();
                        continue;
                    }

                    auto permutations = vector<vector<int>> {};
                    std::string perm_filename = "./include/permutations/dispersive_const_perms/start_state" +
                    std::to_string(start_state) + "_target_" + std::to_string(target) + "_q" + std::to_string(q) + ".txt";
                    
                    std::ifstream infile_perm(perm_filename);
                    if (infile_perm) {
                        permutations = load_binary_perm(perm_filename);
                    } else {
                        permutations = ladder_permutations(start_state, MAX_STEP, target, q);
                        if (permutations.empty()) {
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
                        q_amp.real += res.real;
                        q_amp.imag += res.imag;
                    }
                    divdiff_clear_up();
                    
                    save_binary_complexEx(q_amp, q_amp_filename);
                    
                    q_amp.real *= energy_coefficient_amp; 
                    q_amp.imag *= energy_coefficient_amp; 

                    divdiff_init();
                    total_amp.real += q_amp.real.get_double();
                    total_amp.imag += q_amp.imag.get_double();
                    divdiff_clear_up();
                }

                // //print the Total amplitude
                // // std::cout << "Total amplitude: " << total_amp << std::endl << std::endl;
                // std::cout <<  "Total amp: "; 
                // print_complex_Ex(total_amp);
                transition_amplitudes.push_back(total_amp);
            }

            // // //print transition amplitudes
            // std::cout << "Transition amplitudes: " << std::endl;
            // int target_nums = 0;
            // for (auto& amp : transition_amplitudes) {
            //     std::cout << target_nums << ": "; 
            //     print_complex_Ex(amp);
            //     target_nums++;
            // }

            // normalize transition amplitude 
            double norm = 0;
            ExExFloat Ex_norm;
            Ex_norm = 0;
            divdiff_init();
            for (auto& amp : transition_amplitudes) {
                Ex_norm += amp.real * amp.real + amp.imag * amp.imag;
            }
            norm = sqrt(Ex_norm.get_double());
            for (auto& amp : transition_amplitudes) {
                amp.real = amp.real / norm;
                amp.imag = amp.imag / norm;
            }
            divdiff_clear_up();
            
            all_transition_amplitudes[qubit].push_back(transition_amplitudes);
            
        //end of time loop
        }

        // //print all transition amplitudes with location in the vector
        // std::cout << "All transition amplitudes: " << std::endl;
        // for (int i = 0; i < all_transition_amplitudes[qubit].size(); i++) {
        //     std::cout << "Time: " << times[i] << std::endl;
        //     for (int j = 0; j < all_transition_amplitudes[qubit][i].size(); j++) {
        //         std::cout << "Target: " << j << " Amplitude: ";
        //         print_complex_Ex(all_transition_amplitudes[qubit][i][j]);
        //     }
        // }
    
    // } //end of qubit loop

    std::cout << "Exiting solve_for_t" << std::endl;

    return all_transition_amplitudes;
}


// calculate error between current state and target state
double cal_error(vector<vector<vector<complex_Ex>>>& curr_state, vector<vector<complex_Ex>>& target_state, int qubit) {
    double error = 0.0;
    //print current state
 
    int size_curr = curr_state[qubit].size();
    int size_target = target_state.size();

    // print the states for last time step
    std::cout << "Current state for qubit " << qubit << " at last time step: " << std::endl;
    for (size_t state = 0; state < curr_state[qubit][size_curr-1].size(); ++state) {
        std::cout << "State " << state << ": ";
        print_complex_Ex(curr_state[qubit][size_curr-1][state]);
    }
    std::cout << "Target state at last time step: " << std::endl;
    for (size_t state = 0; state < target_state[size_target-1].size(); ++state) {
        std::cout << "State " << state << ": ";
        print_complex_Ex(target_state[size_target-1][state]);
    }

    // for (size_t qubit = 0; qubit < curr_state.size(); ++qubit) { //TODO: qubit is not used, so we can remove it
        for (size_t state = 0; state < min(curr_state[qubit][size_curr-1].size(), target_state[size_target-1].size()); ++state) {
            double real_diff = curr_state[qubit][size_curr-1][state].real.get_double() - target_state[size_target-1][state].real.get_double();
            double imag_diff = curr_state[qubit][size_curr-1][state].imag.get_double() - target_state[size_target-1][state].imag.get_double();
            // print the differences
            std::cout << "State " << state << ": ";
            std::cout << "Real diff: " << real_diff << ", Imag diff: " << imag_diff << std::endl;
            // print accumulated error
            // print curr addition
            std::cout << "Current addition to error: " << std::pow(real_diff, 2) + std::pow(imag_diff, 2) << std::endl;
            
            error += real_diff * real_diff + imag_diff * imag_diff;
            // std::cout << "Accumulated error: " << error << std::endl;
        }
    // }//end of qubit loop
    return error;
}




int main(int argc, char* argv[]) {
    //get the parameters
    // initial values
    double orig_amplitude = 0.002;
    int q_max = 10; 
    int max_target = 8;

    int num_pulses = 2;

    int qubit = 1; // default qubit to 0
    std::string final_state_filename = "amplitudes_e_chg_pulse3amp0.002000_q10.json"; 

    // Check if command line arguments are provided
    if(argc > 1 && argc < 3){
        final_state_filename = argv[1]; // get the final state filename from command line arguments
    }

    if(argc > 2){
        orig_amplitude = std::stod(argv[2]);
        q_max = std::stoi(argv[3]);
        max_target = std::stoi(argv[4]);
        num_pulses = std::stoi(argv[5]);
        if (argc > 6) {
            qubit = std::stoi(argv[6]); // get the qubit number from command line arguments
        }
    }

    vector<double> amplitudes = {orig_amplitude*0.5, orig_amplitude, orig_amplitude*2, -1*orig_amplitude*2, -1 * orig_amplitude, -orig_amplitude*0.5};    // double  amplitudes[] = {};
    // for(int i = 0; i < num_pulses*2; i++){
    //     if(i < num_pulses){
    //         amplitudes[i] =  amplitude * (i + 1);
    //     }else{
    //         amplitudes[i] =  -1 * amplitude * (num_pulses*2 - i);
    //     }
    // }

    double final_t = 192;
    double sample_t = final_t / 2 / num_pulses;
    int init_start_state = 0;
    // double init_t = 192;
    // double mid_t = 96;
    //init final_state_ta to [2][1][8]
    
    
    //init all_state_ta 
    vector<vector<vector<vector<complex_Ex>>>> all_state_ta;
    vector<vector<vector<vector<complex_Ex>>>> new_all_state_ta;    
    for(int pulse=0; pulse < num_pulses*2; pulse++){
        vector<vector<vector<complex_Ex>>> final_state_ta;
        // for (int qubit = 0; qubit < 2; qubit++){//TODO: qubit is not used, so we can remove it
            vector<vector<complex_Ex>> qubit_state;
            for (int i = 0; i < sample_t; i++){// times size is 1
                vector<complex_Ex> state;
                for (int j = 0; j < max_target; j++){
                    complex_Ex c;
                    c.real = 0;
                    c.imag = 0;
                    state.push_back(c);
                }
                qubit_state.push_back(state);
            }
            final_state_ta.push_back(qubit_state);
        // } // end of qubit loop
        all_state_ta.push_back(final_state_ta);
        new_all_state_ta.push_back(final_state_ta);
    }

    vector<double> amplitudes_list(amplitudes.size());// to store the amplitudes that give the minimum distance
    double min_error = 1; // initialize to a large value

    
    //load final_state_ta from file
    // std::string final_state_filename = "amplitudes_e_chg_amp0.002000_q10.json"; 
    std::ifstream final_state_file(final_state_filename);
    vector<vector<complex_Ex>> final_state_ta;
    if (final_state_file.is_open()) {
        nlohmann::json j;
        final_state_file >> j;
        for (const auto& vec : j) {
            vector<complex_Ex> state;
            for (const auto& c : vec) {
                complex_Ex comp;
                comp.real = c["real"].get<double>();
                comp.imag = c["imag"].get<double>();
                state.push_back(comp);
            }
            final_state_ta.push_back(state);
        }
    } else {
        std::cerr << "Error opening file: " << final_state_filename << std::endl;
        return 1; // Exit if the file cannot be opened
    }

    // //print final_state_ta
    // std::cout << "Final state from file: " << std::endl;
    // for (int i = 0; i < final_state_ta.size(); i++) {
    //     std::cout << "Target: " << i << std::endl;
    //     for (const auto& c : final_state_ta[i]) {
    //         print_complex_Ex(c);
    //     }
    // }

    for(double& amp : amplitudes) {
        std::cout << "Amplitude: " << amp << std::endl;
        new_all_state_ta[0] = solve_for_t(sample_t, q_max, amp, 0, max_target, qubit);
        // cal error between new_all_state_ta and final_state_ta
        double error = cal_error(new_all_state_ta[0], final_state_ta, qubit);
        //print error
        std::cout << "Error: " << error << std::endl;
        if (error < min_error) {
            min_error = error;
            all_state_ta[0] = new_all_state_ta[0];
            amplitudes_list[0] = amp;
            std::cout << "New minimum error: " << min_error << std::endl;
        }
    }

    cout << "past first solve_for_t" << std::endl;
    cout << endl << endl << endl;


    // vector<vector<vector<complex_Ex>>> curr_state_ta;

    for(int curr_mul = 1; curr_mul < (num_pulses*2) ; curr_mul++){ 
        //run for time 192
       

        min_error = 1; // reset min_error for each time_frame
        // go over all amplitudes
        for(int amp = 0; amp < amplitudes.size(); amp++){
            double amplitude = amplitudes[amp];
            auto curr_state_ta = vector<vector<vector<vector<complex_Ex>>>> {};
            auto new_state_ta = vector<vector<vector<complex_Ex>>> {};
            // set new_state_ta to zeros
            new_state_ta.resize(2); //resize for qubits
            for (int qubit = 0; qubit < 2; qubit++){ //TODO: qubit is not used, so we can remove it
                new_state_ta[qubit].resize(int(sample_t)); //resize for time
                for (int i = 0; i < int(sample_t); i++){
                    new_state_ta[qubit][i].resize(max_target);//resize for target states
                    for (int j = 0; j < max_target; j++){
                        new_state_ta[qubit][i][j].real = 0;
                        new_state_ta[qubit][i][j].imag = 0;
                    }
                }
            }
            //print new_state_ta rested
            std::cout << "New state ta reset for current multiplication: " << curr_mul << std::endl;


            // print amplitude
            std::cout << "Amplitude: " << amplitude << std::endl;

            for (int start_state = 0; start_state < max_target; start_state++){
                std::cout << "Start state: " << start_state << std::endl;
                std::cout << std::endl;
                
                auto curr_state = vector<vector<vector<complex_Ex>>> {};

                curr_state = solve_for_t(sample_t, q_max, amplitude, start_state, max_target, qubit);//TODO: save the entire function result to save time?

                //insert curr state into curr_state_ta
                curr_state_ta.push_back(curr_state);
                               
                //save to final state
                // for(int qubit = 0; qubit < 2; qubit++){ //TODO: qubit is not used, so we can remove it
                    std::cout << "qubit: " << qubit << std::endl;
                    
                    for(int i=0; i < int(sample_t); i++){
                        cout << "Time: " << i << std::endl;
                        for (int j = 0; j < max_target; j++){//curr_state_ta[qubit][0] refers to time 96
                            
                            // //print target
                            // std::cout << "Target: " << j << std::endl;
                            // std::cout << "starting " << start_state << std::endl;

                            // //print curr_state_ta
                            // std::cout << "previos all_state_ta:" << std::endl;
                            // print_complex_Ex(all_state_ta[curr_mul-1][qubit][sample_t - 1][start_state]);
                            
                            // //print curr_state_ta
                            // std::cout << "curr_state_ta:" << std::endl;
                            // print_complex_Ex(curr_state_ta[start_state][qubit][i][j]);
                            
                            //multiply curr_state_ta with all_state_ta
                            curr_state_ta[start_state][qubit][i][j] = complex_mult(curr_state_ta[start_state][qubit][i][j], all_state_ta[curr_mul-1][qubit][sample_t - 1][start_state]);
                            
                            // //print curr_state_ta
                            // std::cout << "curr_state_ta:" << std::endl;
                            // print_complex_Ex(curr_state_ta[start_state][qubit][i][j]);
                            
                            // std::cout << "pre new_all_state_ta:" << std::endl;
                            // print_complex_Ex(new_state_ta[qubit][i][j]);
                            
                            divdiff_init();
                            new_state_ta[qubit][i][j].real += curr_state_ta[start_state][qubit][i][j].real;
                            new_state_ta[qubit][i][j].imag += curr_state_ta[start_state][qubit][i][j].imag;
                            divdiff_clear_up();
                            
                            // //print new_all_state_ta
                            // std::cout << "new_all_state_ta:" << std::endl;
                            // print_complex_Ex(new_state_ta[qubit][i][j]);
                            // // std::cout << "finished" << std::endl;
                                
                                
                        }
                    }
                // }//end of qubit loop
            }
            // calculate error between new_all_state_ta and final_state_ta
            double error = cal_error(new_state_ta, final_state_ta, qubit);
            // double error = cal_error(new_all_state_ta[curr_mul], final_state_ta, qubit);
            //print error
            std::cout << "Error for amplitude " << amplitude << ": " << error << std::endl;
            // if error is less than min_error, then save new_all_state_ta to all_state_ta
            if (error < min_error) {
                min_error = error;
                // save new_all_state_ta to all_state_ta
                all_state_ta[curr_mul] = new_state_ta;
                // all_state_ta[curr_mul] = new_all_state_ta[curr_mul];
                std::cout << "New minimum error: " << min_error << std::endl;
                amplitudes_list[curr_mul] = amplitude;
            }
        }
    }

    // //print final state
    // std::cout << "Final state!!!" << std::endl;
    // for(int qubit = 0; qubit < 1; qubit++){
    //     std::cout << "qubit: " << qubit << std::endl;
    //     for (int j = 0; j < MAX_TARGET; j++){
    //         std::cout << "Target: " << j << std::endl;
    //         print_complex_Ex(final_state_ta[qubit][0][j]);
    //     }
    // }
    
    //save final state to a file
    std::cout << "saving to file" << std::endl;

    //print amplitudes_list
    std::cout << "Amplitudes list: " << std::endl;
    for (const auto& amp : amplitudes_list) {
        std::cout << amp << " ";
    }
    
    // Assuming all_transition_amplitudes is a vector<vector<complex<double>>>
    for(int state=0; state < num_pulses*2; state++){
        nlohmann::json j;
        for(int qubit = 0; qubit < 2; qubit++){
            //clear j
            j.clear();
            for (auto& vec : all_state_ta[state][qubit]) {
                nlohmann::json sub_j;
                for (auto& c : vec) {
                    sub_j.push_back({{"real", c.real.get_double()}, {"imag", c.imag.get_double()}});
                }
                j.push_back(sub_j);
            }
    
            std::string filename;
            if (qubit == 0){
                filename = "amplitudes_g_chg_pulse" + std::to_string(state) + "amp" + std::to_string(orig_amplitude) + "_q" + std::to_string(q_max) + ".json";
            }else{
                filename = "amplitudes_e_chg_pulse" + std::to_string(state) + "amp" + std::to_string(orig_amplitude) + "_q" + std::to_string(q_max) + ".json";
            }
            std::ofstream o(filename);
            o << j << std::endl;
        }

    }

    // // save mid_ta to a file
    // nlohmann::json j_mid;
    // for(int qubit = 0; qubit < 2; qubit++){
    //     j_mid.clear();
    //     for (auto& vec : mid_ta[qubit]) {
    //         nlohmann::json sub_j;
    //         for (auto& c : vec) {
    //             sub_j.push_back({{"real", c.real.get_double()}, {"imag", c.imag.get_double()}});
    //         }
    //         j_mid.push_back(sub_j);
    //     }

    //     std::string filename;
    //     if (qubit == 0){
    //         filename = "mid_amplitudes_g_chg_amp" + std::to_string(amplitude) + "_q" + std::to_string(q_max) + ".json";
    //     }else{
    //         filename = "mid_amplitudes_e_chg_amp" + std::to_string(amplitude) + "_q" + std::to_string(q_max) + ".json";
    //     }
    //     std::ofstream o(filename);
    //     o << j_mid << std::endl;
    // }

    return 0;
}
