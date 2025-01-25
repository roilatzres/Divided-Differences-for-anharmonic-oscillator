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

// #define  MAX_STEP 14
#define  MAX_TARGET 8
#define MAX_Q 10



int main(int argc, char* argv[]) {
    complex<double> i_num(0, 1);
    
    //original parameters
    double multiplier = 4;
    float sigma = 48;
    float timestep = 1;
    double chi = -279e-6 * 2 * M_PI; 
    // double detuning = chi *1e-2;


    //TODO: notice that amp is double of the value to put in antiSymSim
    // initial values
    double amplitude = 0.02;
    int q_max = 12; 
    
    if(argc > 1){
        amplitude = std::stod(argv[1]);
        q_max = std::stoi(argv[2]);
    }

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
    for (double t = timestep; t <= tfinal; t += timestep) {//TODO: check if t=1 is correct (t=0 gives 0 amplitude)
        times.push_back(t);
        //print times
        // std::cout << "Time: " << t << std::endl;
    }


    //calculate the time evolution
    double t;

    // vector of all transition ampli tudes for all times: [time][target_state][amp]
    vector<vector<vector<complex<ExExFloat>>>> all_transition_amplitudes(2);
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
    load_all_coefficients(all_coefficients, start_state, total_steps, MAX_TARGET, q_max);


    // int qubit = 1;
    for (int qubit = 0; qubit <=1; qubit++){

        //loop over all times
        for (int i = 191; i < times.size(); i++) { //TODO: i=0
            // t = times[i];
            t = 192;
            
            //print time
            std::cout << std::endl << "Main-Time: " << t << std::endl;

            //vector of all transition amplitude for a given time
            vector<complex<ExExFloat>> transition_amplitudes;
            
            //loop over all target steps
            for (int target = 0; target <= MAX_TARGET; target++) {
                complex<ExExFloat> total_amp;
                total_amp *= 0;
                // total_amp.imag() = 0;

                //print total_amp
                // std::cout << "Total amplitude: " << total_amp << std::endl << std::endl;
                std::cout <<  "First Total amp: "; 
                total_amp.real().print();
                std::cout << "+ ";
                total_amp.imag().print();
                std::cout << "i" << std::endl; 
                std::cout << std::endl;

                //print target
                std::cout << "Target: " << target << std::endl;

                for(int q = 0; q <= q_max; q++){
                    //print q
                    std::cout << "Q: " << q << std::endl;

                    complex<ExExFloat> q_amp; //TODO: change to EXEXFLOAT
                    q_amp *= 0;
                    std::cout <<  "First q_amp: "; 
                    total_amp.real().print();
                    std::cout << "+ ";
                    total_amp.imag().print();
                    std::cout << "i" << std::endl; 
                    std::cout << std::endl;

                    complex<double> energy_coefficient_amp = cal_neg_i_pow(q) * std::pow(amplitude, q);

                    if (target == 0 && q == 0){
                        total_amp = 1;
                        std::cout <<  "second Total amp: "; 

                        total_amp.real().print();
                        std::cout << "+ ";
                        total_amp.imag().print();
                        std::cout << "i" << std::endl; 
                        std::cout << std::endl;

                        continue;
                    }

                    // std::cout << "hi" << std::endl; 
                    auto coefficients = all_coefficients[target][q];
                    // std::cout << "hi" << std::endl; 

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
                    for (const auto& coefficient : coefficients) {
                        
                        q_amp += cal_divdiff(coefficient, t, q, qubit, chi, omega, target);
                    }
                    //clear divdiff
                    divdiff_clear_up();
                    std::cout << "post DIVDIFF" << std::endl;

                    // std::cout << "will i get stuck?" << std::endl;

                    // // save q_amp to a file
                    // save_complex_to_file(q_amp, q_amp_filename);
                    
                    // print energy_coefficient_amp
                    std::cout << "Energy coefficient: " << energy_coefficient_amp.real() << " + " << energy_coefficient_amp.imag() << "i" << std::endl;
                    //print q_amp
                    std::cout << "Q: " << q << std::endl;
                    std::cout << "Q_amp: " << q_amp.real().get_double() << " + " << q_amp.imag().get_double() << "i" << std::endl;

                    q_amp *= energy_coefficient_amp; 

                    // print the q amplitude
                    cout << q << " Original q_amp: " << endl;
                    q_amp.real().print();
                    cout << "+ ";
                    q_amp.imag().print();
                    cout << "i" << endl;

                    cout << std::endl;
                    
                    //add to total amp
                    total_amp += q_amp;
                    // total_amp.imag() += q_amp.imag();
                }

                //print the Total amplitude
                // std::cout << "Total amplitude: " << total_amp << std::endl << std::endl;
                std::cout <<  "Total amp: "; 
                total_amp.real().print();
                std::cout << "+ ";
                total_amp.imag().print();
                std::cout << "i" << std::endl; 
                std::cout << std::endl;
                transition_amplitudes.push_back(total_amp);
            }

            // //print transition amplitudes
            std::cout << "Transition amplitudes: " << std::endl;
            int target_nums = 0;
            for (const auto& amp : transition_amplitudes) {
                std::cout << target_nums << ": "; 
                amp.real().print();
                std::cout << "+ ";
                amp.imag().print();
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
            for (const auto& amp : transition_amplitudes) {
                Ex_norm += amp.real() * amp.real() + amp.imag() * amp.imag();
            }
            norm = sqrt(Ex_norm.get_double());
            for (auto& amp : transition_amplitudes) {
                amp = amp / norm;
            }
            
            all_transition_amplitudes[qubit].push_back(transition_amplitudes);
            // if (qubit == 0){
            //     all_transition_amplitudes_g.push_back(transition_amplitudes);
            // }else{
            //     all_transition_amplitudes_e.push_back(transition_amplitudes);
            // }
            // all_transition_amplitudes.push_back(transition_amplitudes);

            std::cout << std::endl << std::endl;
        
        //end of time loop
        }


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
                all_transition_amplitudes[qubit][i][j].real().print();
                std::cout << "+ ";
                all_transition_amplitudes[qubit][i][j].imag().print();
                std::cout << "i" << std::endl;
                cout << std::endl;
            }
        }

        // Assuming all_transition_amplitudes is a vector<vector<complex<double>>>
        nlohmann::json j;
        for (const auto& vec : all_transition_amplitudes[qubit]) {
            nlohmann::json sub_j;
            for (const auto& c : vec) {
                sub_j.push_back({{"real", c.real().get_double()}, {"imag", c.imag().get_double()}});
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
