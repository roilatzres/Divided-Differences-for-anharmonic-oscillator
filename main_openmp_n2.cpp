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

// #include <mpi.h>
#include <omp.h> 
#include "include/openmp_utils.h"
#include <chrono>


using std::vector;
using std::complex;

// #define  MAX_STEP 14
#define  MAX_TARGET 4
#define MAX_Q 10



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
    int target_max = 4;

    if(argc > 1){
        amplitude = std::stod(argv[1]);
        q_max = std::stoi(argv[2]);
        target_max = std::stoi(argv[3]);
        
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
    for (int target_step = 0; target_step <= MAX_TARGET; ++target_step) {//TODO: change to target_max?
        all_coefficients.push_back(vector<vector<std::tuple<vector<DivdiffElement>, double, int, vector<int>>>>());

        for (int q = 0; q <= MAX_Q; ++q) {//TODO: change to q_max?
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

    

    // // vector of all transition amplitudes for all times: [qubit][time][target_state][q]
    // // in cell q=0 (i.e [qubit][time][target][0]) we get the sum of all q's in the specified target
    // vector<vector<vector<vector<complex<double>>>>> all_transition_amplitudes(2,//qubit = 0/1
    //        vector<vector<vector<complex<double>>>>(times.size(),                         //times 0-191
    //               vector<vector<complex<double>>>(target_max, 
    //                      vector<complex<double>>(q_max) )));

    //test for double nested
    vector<vector<vector<complex<double>>>> all_transition_amplitudes(2,//qubit = 0/1
                  vector<vector<complex<double>>>(target_max+1, 
                         vector<complex<double>>(q_max) ));
    
    
    omp_set_nested(1);

    // Start measuring time
    auto start = std::chrono::high_resolution_clock::now();


    // int qubit = 1;
    #pragma omp parallel for schedule(dynamic)
    for (int qubit = 0; qubit <=1; qubit++){

        // //loop over all times
        // #pragma omp parallel for schedule(dynamic)
        // for (int i = 0; i < times.size(); i++) { 
            // t = times[i];
            t = 192; //TODO: create parallelizem here

            // //print time
            // std::cout << std::endl << "Main-Time: " << t << std::endl;

            // //vector of all transition amplitude for a given time
            // vector<complex<double>> transition_amplitudes(target_max);

            double norm = 0;

            #pragma omp parallel for reduction(+:norm) schedule(dynamic)
            for(int target=0; target <= target_max; target++){
                
                double total_amp_real = 0.0;
                double total_amp_imag = 0.0;
            


                #pragma omp parallel for reduction(+: total_amp_real, total_amp_imag) schedule(static)
                for(int q = 0; q <= q_max; q++){
                    // all_transition_amplitudes[qubit][target][q] = mpi_compute_q_amplitude(target, q, qubit, amplitude, t, chi, omega, total_amp, all_coefficients);
                     complex<double> local_amp = mpi_compute_q_amplitude(target, q, qubit, amplitude, t, chi, omega, all_coefficients);

                    // Store the amplitude in the shared array
                    all_transition_amplitudes[qubit][target][q] = local_amp;

                    // Accumulate the total amplitude
                    #pragma omp atomic
                    total_amp_real += local_amp.real();
                    #pragma omp atomic
                    total_amp_imag += local_amp.imag();
                }
                //TODO: add sync mech here and sum all q
                // }//pragma q
                complex<double> total_amp(total_amp_real, total_amp_imag);

                
                all_transition_amplitudes[qubit][target][0] = total_amp;
                norm += std::norm(total_amp);//TODO: check if norm works with complex<double>
                                            //FIXME: watch out for the line above if debug is required
            }

            //normalize transition amp for each time
            norm = sqrt(norm);

            for(int target=0; target <= target_max; target++){
                all_transition_amplitudes[qubit][target][0] /= norm;
            }
            // }//pragma target
    
        // }//end times loop
        // }//pragma times
    }//end qubit loop
    // }//pragma qubit
    


    // Stop measuring time
    auto end = std::chrono::high_resolution_clock::now();

    // Calculate the duration
    std::chrono::duration<double> elapsed = end - start;

    // Print the result
    std::cout << "Execution time: " << elapsed.count() << " seconds" << std::endl;
    

    //print all transition amplitudes with location in the vector and save to json
    std::cout << "All transition amplitudes: " << std::endl;
    for (int qubit_p = 0; qubit_p < all_transition_amplitudes.size(); qubit_p++) {
        std::cout << "qubit: " << qubit_p  << std::endl;
        
        nlohmann::json j;
        // for (int time_p = 0; time_p < all_transition_amplitudes[qubit_p].size(); time_p++) {
            int time_p = 192;
            std::cout << "Time: " << time_p << std::endl;

            nlohmann::json sub_j;

            
            for(int target_p = 0; target_p <= target_max; target_p++){
                complex<double> final_amp = all_transition_amplitudes[qubit_p][target_p][0];
                std::cout << "Target: " << target_p << " Amplitude: " << final_amp << std::endl;
                sub_j.push_back({{"real", final_amp.real()}, {"imag", final_amp.imag()}});
                // if(qubit_p == 0){
                //     sub0.push_back({{"real", final_amp.real()}, {"imag", final_amp.imag()}});
                // }else{
                //     sub1.push_back({{"real", final_amp.real()}, {"imag", final_amp.imag()}});
                    
                // }

                // //for debug
                // for(int q_amp_p = 0; q_amp_p <= q_max; q++){
                //     if(target_p == 0){
                //         if(q_amp_p == 0){
                //             cout << "q_amp: " << 1 << endl;
                //         }else{
                //             cout << "q_amp: " << 0 << endl;
                //         }
                //         cout << "Q_amp: " << all_transition_amplitudes[qubit_p][time_p][q_amp_p] << endl;

                //     }
                // }
                cout <<"hi" << endl;
            }
            j.push_back(sub_j);
            // if(qubit_p ==0){
            //     j0.push_back(sub0);
            // }else{
            //     j1.push_back(sub1);
            // }
        
        // }//time loop
        // Determine the filename based on the qubit value
        std::string filename = (qubit_p == 0 ? "amplitudes_g_chg_amp" : "amplitudes_e_chg_amp") + 
                            std::to_string(amplitude) + "_q" + std::to_string(q_max) + ".json";
        
        cout <<"hi" << endl;
        // Write the JSON object to the file
        std::ofstream o(filename);
        if (o.is_open()) {
            o << j << std::endl;
            // if(qubit_p == 0){
            // }else{
            //     o << sub1 << std::endl;
            // }
        } else {
            std::cerr << "Failed to open file: " << filename << std::endl;
        }
        // cout <<"bye" << endl;

    }

    cout << "finished" << endl;
    return 0;
}





