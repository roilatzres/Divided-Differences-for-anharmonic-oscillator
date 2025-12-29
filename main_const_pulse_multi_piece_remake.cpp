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


// int main(int argc, char* argv[]) {
vector<vector<vector<complex_Ex>>> solve_for_t(double tfinal, int t_k, double sample_t, int q_max, double amplitude, int start_state, int max_target, double chi, double omega){
    complex<double> i_num(0, 1);
    
    //TODO: notice that amp is double of the value to put in antiSymSim for sin pulse (epsilon/2) = amp_antiSymSim
    
    //print args
    std::cout << "Amplitude: " << amplitude << std::endl;
    std::cout << "q_max: " << q_max << std::endl;

    vector<double> times;
    for (double t = 1; t <= tfinal; t += 1) {
        times.push_back(t);
        //print times
        // std::cout << "Time: " << t << std::endl;
    }

    // //calculate the time evolution
    double delta_t;

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
    for (int qubit = 0; qubit < 2; qubit++){

        // //loop over all times
        for (int i = 0; i < sample_t; i++) {
            delta_t = times[i];
            //print time
            std::cout << std::endl << "sample_t: " << delta_t << std::endl;
            std::cout << std::endl << "t_k: " << t_k << std::endl;


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
                    //print q
                    std::cout << "q: " << q << std::endl;
                    complex_Ex q_amp; 
                    q_amp.real = 0;
                    q_amp.imag = 0;
                    
                    complex<double> energy_coefficient_amp = cal_neg_i_pow(q) * std::pow(amplitude, q);

                    if (target == start_state && q == 0){
                        complex<double> exponent = std::exp(-i_num * delta_t * (start_state * chi * qubit)); 
                        total_amp.real = exponent.real();
                        total_amp.imag = exponent.imag();
                        continue;
                    }

                    // Ensure the folder exists before using q_amp_filename
                    std::string q_amp_folder = "./include/q_amp/multipiece_td/qubit_"  + std::to_string(qubit) 
                        + "_time_" + std::to_string((int)t_k)
                        + "_sample_" + std::to_string((int)sample_t);
                    
                    if (!std::filesystem::exists(q_amp_folder)) {
                        std::filesystem::create_directories(q_amp_folder);
                    }


                    std::string q_amp_filename = q_amp_folder + "/dispersive_q_amp_start_state_" 
                        + std::to_string(start_state) 
                        + "_target_" + std::to_string(target) 
                        + "_q" + std::to_string(q) + ".bin";
                    
                    std::ifstream infile_q_amp(q_amp_filename);
                    if (infile_q_amp) {
                        // cout << "Loading q_amp from file: " << q_amp_filename << std::endl;
                        
                        divdiff_init();

                        // auto q_amp_vec_load = load_binary_q_amp_tuple_disp(q_amp_filename);
                        q_amp = load_binary_complexEx(q_amp_filename);
                        print_complex_Ex(q_amp);
                        


                        
                        complex_mult(q_amp, energy_coefficient_amp);
                        // q_amp.real *= energy_coefficient_amp; 
                        // q_amp.imag *= energy_coefficient_amp; 

                        // print q_amp
                        std::cout << "full q_amp: ";
                        print_complex_Ex(q_amp);
                        
                        double q_real = q_amp.real.get_double();
                        double q_imag = q_amp.imag.get_double();
                        total_amp.real += q_real;
                        total_amp.imag += q_imag;
                        divdiff_clear_up();
                        continue;
                    }

                    auto permutations = vector<vector<int>> {};
                    std::string perm_folder = "./include/permutations/multipiece_td";
                    if (!std::filesystem::exists(perm_folder)) {
                        std::filesystem::create_directories(perm_folder);
                    }

                    std::string perm_filename = perm_folder + "/start_state" +
                    std::to_string(start_state) + "_target_" + std::to_string(target) + "_q" + std::to_string(q) + ".bin";
                    
                    std::ifstream infile_perm(perm_filename);
                    if (infile_perm) {
                        // cout << "Loading permutations from file: " << perm_filename << std::endl;
                        permutations = load_binary_perm(perm_filename);
                    } else {
                        // cout << "Generating permutations and saving to file: " << perm_filename << std::endl;
                        permutations = ladder_permutations(start_state, MAX_STEP, target, q);
                        if (permutations.empty()) {
                            cout << "No permutations found for start_state " << start_state 
                                 << ", target " << target << ", q " << q << ". Skipping." << std::endl;
                            continue;
                        } else {
                            save_binary_perm(permutations, perm_filename);
                        }
                    }
                       
                    vector<std::tuple<vector<DivdiffElement>, double, int, vector<int>>> all_coefficients_td;
                    std::string all_coef_folder = "./include/parms/multipiece_td";
                    if (!std::filesystem::exists(all_coef_folder)) {
                        std::filesystem::create_directories(all_coef_folder);
                    }

                    std::string all_coef_filename = all_coef_folder 
                    + "/start_state" + std::to_string(start_state) 
                    + "_target_" + std::to_string(target) 
                    + "_q" + std::to_string(q) + ".bin";
                   
                    std::ifstream infile_coef(all_coef_filename);
                    if (infile_coef) {
                        cout << "Loading all coefficients from file: " << all_coef_filename << std::endl;
                        all_coefficients_td = load_binary_orig(all_coef_filename);
                    } else {
                        cout << "Calculating all coefficients and saving to file: " << all_coef_filename << std::endl;
                        for (const auto& permutation : permutations) {
                            auto coefficients = cal_coefficient(permutation, start_state);
                            all_coefficients_td.push_back(coefficients);
                        }
                        save_binary_orig(all_coefficients_td, all_coef_filename);
                    }

                    vector<std::tuple<complex_Ex, int>> q_amp_vec;
                    
                    divdiff_init();
                    for (const auto& coefficient : all_coefficients_td) {
                        cout << "Calculating divdiff for a permutation..." << std::endl;
                        complex_Ex res = cal_divdiff_td(coefficient, delta_t, q, qubit, chi, omega, target, t_k);//TODO: check this function
                        // auto final_state = std::get<2>(coefficient);
                        // q_amp_vec.emplace_back(res, final_state);
                        // cout << "Intermediate res: ";
                        // print_complex_Ex(res);
                        // res = complex_mult(res, start_time_phase);
                        // cout << "After adding start_time_phase, res: ";
                        // print_complex_Ex(res);                        

                        q_amp.real += res.real;
                        q_amp.imag += res.imag;
                        cout << "Current q_amp: ";
                        print_complex_Ex(q_amp);
                    }
                    divdiff_clear_up();
                    std::cout << "Finished calculating q_amp" << std::endl;
                    //print q_amp
                    std::cout << "q_amp before time phase: ";
                    print_complex_Ex(q_amp);
                    save_binary_complexEx(q_amp, q_amp_filename);
                    
                    
                    print_complex_Ex(q_amp);
                    // save_binary_q_amp_tuple_disp(q_amp_vec, q_amp_filename);
                    
                    complex_mult(q_amp, energy_coefficient_amp);
                    // q_amp.real *= energy_coefficient_amp; 
                    // q_amp.imag *= energy_coefficient_amp; 

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

        //print all transition amplitudes with location in the vector
        std::cout << "All transition amplitudes: " << std::endl;
        for (int i = 0; i < all_transition_amplitudes[qubit].size(); i++) {
            std::cout << "Time: " << times[i] << std::endl;
            for (int j = 0; j < all_transition_amplitudes[qubit][i].size(); j++) {
                std::cout << "Target: " << j << " Amplitude: ";
                print_complex_Ex(all_transition_amplitudes[qubit][i][j]);
            }
        }
    
    }

    std::cout << "Exiting solve_for_t" << std::endl;

    return all_transition_amplitudes;
}

// // condisp_pulse1 = amp *  np.sin((times) * detuning) * np.exp(-1j*omega_shift*times) 
// vector<double> cal_amps(double amplitude, double chi, double t, double detuning, double sample_t){
//     complex<double> i_num(0, 1);
    
//     vector<double> amps;

//     for(int i = 1; i <= t+1; i+=sample_t){
//         double sin_omega_t = amplitude * std::sin(detuning * i );
//         // //print sin_omega_t
//         // std::cout << "sin_omega_t: " << sin_omega_t << std::endl;
//         // complex<double> expo1 = std::exp(-i_num * chi * (double)i);
//         // complex<double> expo2 = std::exp(i_num * chi * (double)i);
//         // //print expo
//         // std::cout << "expo: " << expo << std::endl;
//         // complex<double> amp = sin_omega_t * (expo1+expo2);
//         double amp = sin_omega_t;
//         amps.push_back(amp);
//     }
//     return amps;
// }



int main(int argc, char* argv[]) {
    //get the parameters
    // initial values
    double amplitude = 0.002;
    int q_max = 10; 
    int max_target = 8;

    int num_pulses = 1;
    
    std::cout << "Program started\n";
    // Load parameters from file
    nlohmann::json params;
    std::ifstream params_file("params_td.json");
    if (!params_file.is_open()) {
        std::cerr << "Could not open params.json!" << std::endl;
        return 1;
    }
    params_file >> params;


    // Assign parameters
    double final_t = params["final_t"];
    double multiplier = params["multiplier"];
    float sigma = params["sigma"];
    float timestep = params["timestep"];
    double chi = params["chi"]; // multiply by 2*M_PI below
    amplitude = params["amplitude"];
    q_max = params["q_max"];
    max_target = params["max_target"];
    num_pulses = params["num_pulses"];
    std::string out_dir = params["out_dir"];

    // Apply 2*M_PI where needed
    chi *= 2 * M_PI;
    
    // final_t = sigma * multiplier;
    double omega = sigma * multiplier;
    double sample_t = final_t / num_pulses;
    
    // vector<double> all_amps = cal_amps(amplitude, chi, final_t,  2 * M_PI / final_t, sample_t);


   
    //init final_state_ta to [2][1][8]
    
    //init all_state_ta 
    vector<vector<vector<vector<complex_Ex>>>> all_state_ta;
    for(int state=0; state < num_pulses; state++){
        vector<vector<vector<complex_Ex>>> final_state_ta;
        for (int qubit = 0; qubit < 2; qubit++){
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
        }
        all_state_ta.push_back(final_state_ta);
    }


    
    all_state_ta[0] = solve_for_t(final_t, 0, sample_t, q_max, amplitude, 0, max_target, chi, omega);

    // run all states for time 192
    // vector<vector<vector<complex_Ex>>> curr_state_ta;

    for(int curr_mul = 1; curr_mul < (num_pulses) ; curr_mul++){
        //run for time 192
        cout << "Running for time " << curr_mul << std::endl;

        auto curr_state_ta = vector<vector<vector<vector<complex_Ex>>>> {};
        for (int start_state = 0; start_state < max_target; start_state++){
            std::cout << "Start state: " << start_state << std::endl;
            std::cout << std::endl;
    
            auto curr_state = vector<vector<vector<complex_Ex>>> {};
            curr_state = solve_for_t(final_t, curr_mul, sample_t, q_max,  amplitude, start_state, max_target, chi, omega);
            
            //insert curr state into curr_state_ta
            curr_state_ta.push_back(curr_state);
            
            //save to final state
            for(int qubit = 0; qubit < 2; qubit++){
                // std::cout << "qubit: " << qubit << std::endl;
                
                for(int i=0; i < int(sample_t); i++){
                    // cout << "Time: " << i << std::endl;
                    for (int j = 0; j < max_target; j++){//curr_state_ta[qubit][0] refers to time 96
                        // //print target
                        // std::cout << "Target: " << j << std::endl;
                        // std::cout << "starting " << start_state << std::endl;
    
                        // //print curr_state_ta
                        // std::cout << "curr_state_ta:" << std::endl;
                        // print_complex_Ex(curr_state_ta[start_state][qubit][i][j]);
    
                        divdiff_init();
                        curr_state_ta[start_state][qubit][i][j] = complex_mult(curr_state_ta[start_state][qubit][i][j], all_state_ta[curr_mul-1][qubit][sample_t - 1][start_state]);
                        
                        // //print curr_state_ta
                        // std::cout << "new curr_state_ta:" << std::endl;
                        // print_complex_Ex(curr_state_ta[start_state][qubit][i][j]);
                        
                        all_state_ta[curr_mul][qubit][i][j].real += curr_state_ta[start_state][qubit][i][j].real;
                        all_state_ta[curr_mul][qubit][i][j].imag += curr_state_ta[start_state][qubit][i][j].imag;
                        divdiff_clear_up();
    
                        // //print final_state_ta
                        // std::cout << "final_state_ta:" << std::endl;
                        // print_complex_Ex(all_state_ta[curr_mul][qubit][i][j]);
                        
                        // std::cout << "finished" << std::endl;
                    }
                }
            }
        }
    }



    
    // // Check if  folder exists, create if not
    // std::string folder_name = "dipsersive_multi_piece_runs";
    if (!std::filesystem::exists(out_dir)) {
        std::filesystem::create_directory(out_dir);
    }

    //save final state to a file
    std::cout << "saving to file" << std::endl;

    // Create one JSON array for each qubit
    for(int qubit = 0; qubit < 2; qubit++){
        nlohmann::json j_all_states;
        for(int state=0; state < num_pulses; state++){
            nlohmann::json j_state;
            for (auto& vec : all_state_ta[state][qubit]) {
                nlohmann::json sub_j;
                for (auto& c : vec) {
                    sub_j.push_back({{"real", c.real.get_double()}, {"imag", c.imag.get_double()}});
                }
                j_state.push_back(sub_j);
            }
            j_all_states.push_back(j_state);
        }

        std::string filename;
        if (qubit == 0){
            filename = out_dir + "/amplitudes_g_chg_amp" + std::to_string(amplitude) + "_q" + std::to_string(q_max) + "_n" + std::to_string(max_target) + "_p" + std::to_string(num_pulses) + ".json";
        }else{
            filename = out_dir + "/amplitudes_e_chg_amp" + std::to_string(amplitude) + "_q" + std::to_string(q_max) + "_n" + std::to_string(max_target) + "_p" + std::to_string(num_pulses) + ".json";
        }
        std::ofstream o(filename);
        o << j_all_states.dump(4) << std::endl;
    }



    return 0;
}
