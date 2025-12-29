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


// //function to print complex_ex
// void print_complex_Ex(complex_Ex a){
//     a.real.print();
//     std::cout << "+ ";
//     a.imag.print();
//     std::cout << "i" << std::endl; 
//     std::cout << std::endl;
// }


// int main(int argc, char* argv[]) {
vector<vector<vector<complex_Ex>>> solve_for_t(double time, int q_max, double amplitude, int start_state, int max_target){
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
    
    // //define start state
    // // int start_state = 0;
    // int total_steps = 100;

    // Create a 3D vector to store the results
    vector< //target_step
    vector< //q
    vector< //permutations
    std::tuple<vector<DivdiffElement>, double, int, vector<int>>> // coefficients
    >> all_coefficients;
    
    // load_all_coefficients(all_coefficients, start_state, total_steps, MAX_TARGET, q_max);
    for (int qubit = 0; qubit < 2; qubit++){

        // //loop over all times
        // // for (int i = 0; i < times.size(); i++) { 
        for (int i = 0; i < time; i++) { 
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

                // //print target
                // std::cout << "Target: " << target << std::endl;

                for(int q = 0; q <= q_max; q++){
                    complex_Ex q_amp; 
                    q_amp.real = 0;
                    q_amp.imag = 0;

                    double energy_coefficient_amp = std::pow(amplitude, q); //TODO: check if complex amplitude is valid for pow

                    if (target == start_state && q == 0){
                        complex<double> exponent = std::exp(-i_num * t * (start_state * chi * qubit));
                        total_amp.real = exponent.real();
                        total_amp.imag = exponent.imag();
                        continue;
                    }

                    std::string q_amp_filename = "./include/q_amp/dispersive_sin/dispersive_q_amp_start_state_" + std::to_string(start_state) 
                         + "_qubit_" + std::to_string(qubit) + "_t_" + std::to_string((int)t) + "_target_" + std::to_string(target) + "_q" + std::to_string(q) + ".bin";
                    std::ifstream infile_q_amp(q_amp_filename);

                    if (infile_q_amp) {
                        q_amp = load_binary_complexEx(q_amp_filename);
                        q_amp.real *= energy_coefficient_amp; 
                        q_amp.imag *= energy_coefficient_amp; 

                        // // take into account the imag part of the drive - epsilon * a + epsilon^star * a^dagger
                        // complex<double> chi_exp = std::exp(i_num * (double)(target - start_state) * chi * t);
                        // q_amp.real *= chi_exp.real();
                        // q_amp.imag *= chi_exp.imag();

                        divdiff_init();
                        total_amp.real += q_amp.real.get_double();
                        total_amp.imag += q_amp.imag.get_double();
                        divdiff_clear_up();
                        continue;
                    }
                    
                    //TODO: check if same for const pulse and dont need extra folder
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
                       
                    //TODO: check if same for const pulse and dont need extra folder
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
                    
                    // // take into account the imag part of the drive - epsilon * a + epsilon^star * a^dagger
                    // complex<double> chi_exp = std::exp(i_num * (double)(target - start_state) * chi * t);
                    // q_amp.real *= chi_exp.real();
                    // q_amp.imag *= chi_exp.imag();
                    

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
    
    }

    std::cout << "Exiting solve_for_t" << std::endl;

    return all_transition_amplitudes;
}

// # condisp_pulse1 = amp *  np.sin((times) * detuning) * np.exp(-1j*omega_shift*times) 
vector<double> cal_amps(double amplitude, double chi, double t, double detuning, double sample_t){
    complex<double> i_num(0, 1);
    
    vector<double> amps;

    for(int i = 1; i <= t+1; i+=sample_t){
        double sin_omega_t = amplitude * std::sin(detuning * i );
        // //print sin_omega_t
        // std::cout << "sin_omega_t: " << sin_omega_t << std::endl;
        // complex<double> expo1 = std::exp(-i_num * chi * (double)i);
        // complex<double> expo2 = std::exp(i_num * chi * (double)i);
        // //print expo
        // std::cout << "expo: " << expo << std::endl;
        // complex<double> amp = sin_omega_t * (expo1+expo2);
        double amp = sin_omega_t;
        amps.push_back(amp);
    }
    return amps;
}




int main(int argc, char* argv[]) {
    //get the parameters
    // initial values
    double amplitude = 0.002;
    int q_max = 10; 
    int max_target = 8;

    int num_pulses = 1;
    
    if(argc > 1){
        amplitude = std::stod(argv[1]);
        q_max = std::stoi(argv[2]);
        max_target = std::stoi(argv[3]);
        num_pulses = std::stoi(argv[4]);
    }

    //original parameters
    double multiplier = 4;
    float sigma = 48;
    float timestep = 1;
    double chi = -279e-6 * 2 * M_PI; 
    
    
    double init_t = sigma * multiplier;
    double sample_t = init_t / num_pulses;
    double final_t = init_t;
    // vector<double> all_amps = {};
    vector<double> all_amps = cal_amps(amplitude, chi, final_t,  2 * M_PI / final_t, sample_t);


    printf("Amplitudes: \n");
    for (int i = 0; i < all_amps.size(); i++) {
        std::cout << "Amplitude " << i << ": " << all_amps[i] << std::endl;
    }

    //init final_state_ta to [2][1][8]
    
    //init all_state_ta 
    vector<vector<vector<vector<complex_Ex>>>> all_state_ta;
    for(int state=0; state <= num_pulses; state++){
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


    
    all_state_ta[0] = solve_for_t(sample_t, q_max, all_amps[0], 0, max_target);

    // run all states for time 192
    // vector<vector<vector<complex_Ex>>> curr_state_ta;

    for(int curr_mul = 1; curr_mul <= num_pulses ; curr_mul++){
        //run for time 192
        cout << std::endl;
        cout << "Running for time " << curr_mul << std::endl;
        auto curr_state_ta = vector<vector<vector<vector<complex_Ex>>>> {};
        for (int start_state = 0; start_state < max_target; start_state++){
            // std::cout << "Start state: " << start_state << std::endl;
            std::cout << std::endl;
    
            auto curr_state = vector<vector<vector<complex_Ex>>> {};
            curr_state = solve_for_t(sample_t, q_max, all_amps[curr_mul], start_state, max_target);//TODO: add total_t for the cal of the divdiff
            
            //insert curr state into curr_state_ta
            curr_state_ta.push_back(curr_state);
            
            //save to final state
            for(int qubit = 0; qubit < 2; qubit++){
                // std::cout << "qubit: " << qubit << std::endl;
                
                for(int i=0; i < int(sample_t); i++){
                    // cout << "Time: " << i << std::endl;
                    for (int j = 0; j < max_target; j++){//curr_state_ta[qubit][0] refers to time 96
                        // //print target
                        std::cout << "Target: " << j << std::endl;
                        std::cout << "starting " << start_state << std::endl;
    
                        // // //print curr_state_ta
                        // std::cout << "curr_state_ta:" << std::endl;
                        // print_complex_Ex(curr_state_ta[start_state][qubit][i][j]);
    
                        curr_state_ta[start_state][qubit][i][j] = complex_mult(curr_state_ta[start_state][qubit][i][j], all_state_ta[curr_mul-1][qubit][sample_t - 1][start_state]);
                        
                        // //print curr_state_ta
                        // std::cout << "new curr_state_ta:" << std::endl;
                        // print_complex_Ex(curr_state_ta[start_state][qubit][i][j]);
                        
                        divdiff_init();
                        all_state_ta[curr_mul][qubit][i][j].real += curr_state_ta[start_state][qubit][i][j].real;
                        all_state_ta[curr_mul][qubit][i][j].imag += curr_state_ta[start_state][qubit][i][j].imag;
                        divdiff_clear_up();
    
                        //print final_state_ta
                        std::cout << "final_state_ta:" << std::endl;
                        print_complex_Ex(all_state_ta[curr_mul][qubit][i][j]);
                        
                        std::cout << "finished" << std::endl;
                    }
                }
            }
        }
    }


    
    //save final state to a file
    std::cout << "saving to file" << std::endl;
    
    // Assuming all_transition_amplitudes is a vector<vector<complex<double>>>
    for(int state=0; state < num_pulses; state++){
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
                filename = "sin_test4/amplitudes_g_chg_pulse" + std::to_string(state) + "amp" + std::to_string(amplitude) + "_q" + std::to_string(q_max) + "_n" + std::to_string(max_target) + "_p" + std::to_string(num_pulses)+ ".json";
            }else{
                filename = "sin_test4/amplitudes_e_chg_pulse" + std::to_string(state) + "amp" + std::to_string(amplitude) + "_q" + std::to_string(q_max) + "_n" + std::to_string(max_target) + "_p" + std::to_string(num_pulses)+ ".json";
            }
            std::ofstream o(filename);
            o << j << std::endl;
        }

    }



    return 0;
}
