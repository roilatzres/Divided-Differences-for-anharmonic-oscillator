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

#define  MAX_STEP 100
#define  MAX_TARGET 6
#define MAX_Q 10





// int main(int argc, char* argv[]) {
vector<vector<vector<complex_Ex>>> solve_for_t(double t, int qubit, int q_max, double amplitude, int start_state){
    complex<double> i_num(0, 1);
    
    //original parameters
    double multiplier = 4;
    float sigma = 48;
    float timestep = 1;
    double chi = -279e-6 * 2 * M_PI; 
    // double detuning = chi *1e-2;


    //TODO: notice that amp is double of the value to put in antiSymSim for sin pulse
    
    // // initial values
    // double amplitude = 0.02;
    // int q_max = 12; 
    
    // if(argc > 1){
    //     amplitude = std::stod(argv[1]);
    //     q_max = std::stoi(argv[2]);
    // }

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
    // double t;

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
    

    // int qubit = 1;
    for (int qubit = 0; qubit <=1; qubit++){

        // //loop over all times
        // // for (int i = 0; i < times.size(); i++) { //TODO: i=0
        // for (int i = 0; i < 97; i++) { //TODO: i=0
        //     t = times[i];
        //     // t = 96;
            
            //print time
            std::cout << std::endl << "Main-Time: " << t << std::endl;

            //vector of all transition amplitude for a given time
            vector<complex_Ex> transition_amplitudes;
            
            //loop over all target steps
            for (int target = 0; target < MAX_TARGET; target++) {
                complex_Ex total_amp;
                total_amp.real = 0;
                total_amp.imag = 0;
                // total_amp.imag() = 0;

                //print total_amp
                // std::cout << "Total amplitude: " << total_amp << std::endl << std::endl;
                std::cout <<  "init Total amp: "; 
                total_amp.real.print();
                std::cout << "+ ";
                total_amp.imag.print();
                std::cout << "i" << std::endl; 
                std::cout << std::endl;

                //print target
                std::cout << "Target: " << target << std::endl;

                for(int q = 0; q <= q_max; q++){
                    //print q
                    std::cout << "Q: " << q << std::endl;

                    complex_Ex q_amp; 
                    q_amp.real = 0;
                    q_amp.imag = 0;
                    std::cout <<  "init q_amp: "; 
                    q_amp.real.print();
                    std::cout << "+ ";
                    q_amp.imag.print();
                    std::cout << "i" << std::endl; 
                    std::cout << std::endl;

                    double energy_coefficient_amp = std::pow(amplitude, q);

                    if (target == 0 && q == 0){
                        total_amp.real = 1;
                        std::cout <<  "Total amp target,q = 0: "; 

                        total_amp.real.print();
                        std::cout << "+ ";
                        total_amp.imag.print();
                        std::cout << "i" << std::endl; 
                        std::cout << std::endl;

                        continue;
                    }

                    // auto coefficients = all_coefficients[target][q];
                    vector<vector<int>> permutations = ladder_permutations(start_state, MAX_STEP, target, q);//TODO:validate start_state
                    vector<std::tuple<double, int, vector<int>>> all_coefficients_const_pulse;

                     for (const auto& permutation : permutations) {
                        auto coefficients = cal_coefficient_const_pulse(permutation, start_state);
                        all_coefficients_const_pulse.push_back(coefficients);
                     }

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
                    
                    // // check if file for q_amp already exists
                    // std::string q_amp_filename = "./include/parms/dispersive/q_amp_recursive/dispersive_q_amp_qubit_" 
                    //     + std::to_string(qubit) + "_t_" + std::to_string((int)t) + "_target_" + std::to_string(target) + "_q" + std::to_string(q) + ".bin";
                    // std::ifstream infile(q_amp_filename);
                    // if (infile) {//TODO: remove q != 12
                    //     // std::cout << "File " << q_amp_filename << " already exists." << std::endl;
                        
                    //     // load q_amp from file
                    //     q_amp = load_complex_from_file(q_amp_filename);

                    //     q_amp = q_amp * energy_coefficient_amp;
                        
                    //     // if(q > 12){ *************************************
                    //         // print the q amplitude
                    //         std::cout << q << " Original q_amp: " << q_amp << std::endl;
                    //     // }***********************************************
                    //     //add to total amp
                    //     total_amp += q_amp; 
                    //     continue;
                    // }

                    std::cout << "pre DIVDIFF" << std::endl;
                    // Need to calculate q_amp
                    divdiff_init();
                    // std::cout << coefficients.size() << std::endl; 
                    for (const auto& coefficient : all_coefficients_const_pulse) {
                        complex_Ex res = cal_divdiff_const_amp(coefficient, t, q, qubit, chi, omega, target);
                        q_amp.real += res.real;
                        q_amp.imag += res.imag;
                    }
                    //clear divdiff
                    divdiff_clear_up();
                    std::cout << "post DIVDIFF" << std::endl;

                    // std::cout << "will i get stuck?" << std::endl;

                    // // save q_amp to a file
                    // save_complex_to_file(q_amp, q_amp_filename);
                    
                    // print energy_coefficient_amp
                    // std::cout << "Energy coefficient: " << energy_coefficient_amp.real() << " + " << energy_coefficient_amp.imag() << "i" << std::endl;
                    std::cout << "Energy coefficient: " << energy_coefficient_amp << std::endl;
                    //print q_amp
                    std::cout << "Q: " << q << std::endl;
                    std::cout << "Q_amp: " << q_amp.real.get_double() << " + " << q_amp.imag.get_double() << "i" << std::endl;


                    q_amp.real *= energy_coefficient_amp; 
                    q_amp.imag *= energy_coefficient_amp; 

                    // print the q amplitude
                    cout << q << " Original q_amp: " << endl;
                    q_amp.real.print();
                    cout << "+ ";
                    q_amp.imag.print();
                    cout << "i" << endl;

                    cout << std::endl;
                    
                    //add to total amp
                    //print hi 
                    std::cout << "HI" << std::endl;
                    total_amp.real += q_amp.real;
                    total_amp.imag += q_amp.imag;
                    
                    std::cout << "HI23" << std::endl;
                    
                    // total_amp.real() = q_amp_real;
                    // total_amp.imag() += q_amp_imag;
                    
                    // total_amp.real() += q_amp.real();
                    // total_amp.imag() += q_amp.imag();
                    // std::cout << "HI" << std::endl;
                    
                    //print the Total amplitude
                    std::cout <<  "Total amp: "; 
                    total_amp.real.print();
                    std::cout << "+ ";
                    total_amp.imag.print();
                    std::cout << "i" << std::endl; 
                    std::cout << std::endl;
                }

                //print the Total amplitude
                // std::cout << "Total amplitude: " << total_amp << std::endl << std::endl;
                std::cout <<  "Total amp: "; 
                total_amp.real.print();
                std::cout << "+ ";
                total_amp.imag.print();
                std::cout << "i" << std::endl; 
                std::cout << std::endl;
                transition_amplitudes.push_back(total_amp);
            }

            // //print transition amplitudes
            std::cout << "Transition amplitudes: " << std::endl;
            int target_nums = 0;
            for (auto& amp : transition_amplitudes) {
                std::cout << target_nums << ": "; 
                amp.real.print();
                std::cout << "+ ";
                amp.imag.print();
                std::cout << "i" << std::endl; 
                std::cout << std::endl;
                target_nums++;
            }

            // // save the unnormailzed transition amplitudes
            // reg_trans_amp = transition_amplitudes;
            
            // normalize transition amplitude 
            double norm = 0;
            ExExFloat Ex_norm;
            Ex_norm = 0;
            for (auto& amp : transition_amplitudes) {
                Ex_norm += amp.real * amp.real + amp.imag * amp.imag;
            }
            norm = sqrt(Ex_norm.get_double());
            for (auto& amp : transition_amplitudes) {
                amp.real = amp.real / norm;
                amp.imag = amp.imag / norm;
            }
            
            all_transition_amplitudes[qubit].push_back(transition_amplitudes);
            // if (qubit == 0){
            //     all_transition_amplitudes_g.push_back(transition_amplitudes);
            // }else{
            //     all_transition_amplitudes_e.push_back(transition_amplitudes);
            // }
            // all_transition_amplitudes.push_back(transition_amplitudes);

            std::cout << std::endl << std::endl;
        // //end of time loop
        // }


        // if (qubit == 0){
        //     all_transition_amplitudes = all_transition_amplitudes_g;
        // }else{
        //     all_transition_amplitudes = all_transition_amplitudes_e;
        // }

        //print all transition amplitudes with location in the vector
        std::cout << "All transition amplitudes: " << std::endl;
        for (int i = 0; i < all_transition_amplitudes[qubit].size(); i++) {
            std::cout << "Time: " << times[i] << std::endl;
            for (int j = 0; j < all_transition_amplitudes[qubit][i].size(); j++) {
                std::cout << "Target: " << j << " Amplitude: ";
                all_transition_amplitudes[qubit][i][j].real.print();
                std::cout << "+ ";
                all_transition_amplitudes[qubit][i][j].imag.print();
                std::cout << "i" << std::endl;
                cout << std::endl;
            }
        }
    }
    return all_transition_amplitudes;
}




int main(int argc, char* argv[]) {
    //get the parameters
    // initial values
    double amplitude = 0.02;
    int q_max = 12; 
    
    if(argc > 1){
        amplitude = std::stod(argv[1]);
        q_max = std::stoi(argv[2]);
    }
    double mid_t = 96;
    double final_t = 192;

    //init final_state_ta to [2][1][8]
    vector<vector<vector<complex_Ex>>> final_state_ta;
    for (int qubit = 0; qubit < 2; qubit++){
        vector<vector<complex_Ex>> qubit_state;
        for (int i = 0; i < 1; i++){// times size is 1
            vector<complex_Ex> state;
            for (int j = 0; j < MAX_TARGET; j++){
                complex_Ex c;
                c.real = 0;
                c.imag = 0;
                state.push_back(c);
            }
            qubit_state.push_back(state);
        }
        final_state_ta.push_back(qubit_state);
    }

    std::cout << "Mid state!!!" << std::endl;
    //run for time 96
    vector<vector<vector<complex_Ex>>> mid_ta = solve_for_t(mid_t, 0, q_max, amplitude, 0);
    //get state

    // final_state_ta = mid_ta;//TODO: check how to init to 0

    // run all states for time 192
    for (int start_state = 0; start_state < MAX_TARGET; start_state++){
        //print start state
        std::cout << "Start state: " << start_state << std::endl;
        std::cout << std::endl;
            
        vector<vector<vector<complex_Ex>>> curr_state_ta;

        curr_state_ta = solve_for_t(final_t - mid_t, 0, q_max, (-1) * amplitude, start_state);//TODO: check t - 192/96/95?
        //save to final state
        for(int qubit = 0; qubit < 2; qubit++){
            std::cout << "qubit: " << qubit << std::endl;

            for (int j = 0; j < MAX_TARGET; j++){//curr_state_ta[qubit][0] refers to time 96
                //print target
                std::cout << "Target: " << j << std::endl;
                std::cout << "starting " << start_state << std::endl;

                //print mid_ta 
                std::cout << "mid_ta:" << std::endl;
                    mid_ta[qubit][0][start_state].real.print();
                    std::cout << "+ ";
                    mid_ta[qubit][0][start_state].imag.print();
                    std::cout << "i" << std::endl;
                    cout << std::endl;

                //print curr_state_ta
                std::cout << "curr_state_ta:" << std::endl;
                    curr_state_ta[qubit][0][j].real.print();
                    std::cout << "+ ";
                    curr_state_ta[qubit][0][j].imag.print();
                    std::cout << "i" << std::endl;
                    cout << std::endl;

                curr_state_ta[qubit][0][j] = complex_mult(curr_state_ta[qubit][0][j], mid_ta[qubit][0][start_state]);
                
                //print curr_state_ta
                std::cout << "new curr_state_ta:" << std::endl;
                    curr_state_ta[qubit][0][j].real.print();
                    std::cout << "+ ";
                    curr_state_ta[qubit][0][j].imag.print();
                    std::cout << "i" << std::endl;
                    cout << std::endl;
                
                final_state_ta[qubit][0][j].real += curr_state_ta[qubit][0][j].real;
                final_state_ta[qubit][0][j].imag += curr_state_ta[qubit][0][j].imag;
                std::cout << "finished" << std::endl;
            }
        }
    }

    //print final state
    std::cout << "Final state!!!" << std::endl;
    for(int qubit = 0; qubit < 2; qubit++){
        std::cout << "qubit: " << qubit << std::endl;
        for (int j = 0; j < MAX_TARGET; j++){
            std::cout << "Target: " << j << std::endl;
            final_state_ta[qubit][0][j].real.print();
            std::cout << "+ ";
            final_state_ta[qubit][0][j].imag.print();
            std::cout << "i" << std::endl;
            cout << std::endl;
        }
    }
    
    //save final state to a file
    std::cout << "saving to file" << std::endl;
    
    // Assuming all_transition_amplitudes is a vector<vector<complex<double>>>
    
    nlohmann::json j;
    for(int qubit = 0; qubit < 2; qubit++){
        // for (auto& vec : mid_ta[qubit]) {
        for (auto& vec : final_state_ta[qubit]) {
            nlohmann::json sub_j;
            for (auto& c : vec) {
                sub_j.push_back({{"real", c.real.get_double()}, {"imag", c.imag.get_double()}});
            }
            j.push_back(sub_j);
        }

        std::string filename;
        if (qubit == 0){
            filename = "amplitudes_g_chg_amp" + std::to_string(amplitude) + "_q" + std::to_string(q_max) + ".json";
        }else{
            filename = "amplitudes_e_chg_amp" + std::to_string(amplitude) + "_q" + std::to_string(q_max) + ".json";
        }
        std::ofstream o(filename);
        o << j << std::endl;
    }

    return 0;
}
