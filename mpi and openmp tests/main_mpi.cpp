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
#include <mpi.h>
#include "include/mpi_utils.h"

using std::vector;
using std::complex;

// #define  MAX_STEP 14
#define  MAX_TARGET 6
#define MAX_Q 8



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
    int target_max = 10;

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

    // vector of all transition ampli tudes for all times: [time][target_state][amp]
    vector<vector<vector<complex<double>>>> all_transition_amplitudes;
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

    MPI_Init(&argc, &argv);

    int world_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);


    // int qubit = 1;
    // for (int qubit = 0; qubit <=1; qubit++){

    int qubit_group = world_rank / 2;
    MPI_Comm qubit_comm;
    MPI_Comm_split(MPI_COMM_WORLD, qubit_group, world_rank, &qubit_comm);

    int qubit_rank, qubit_size;
    MPI_Comm_rank(qubit_comm, &qubit_rank);
    MPI_Comm_size(qubit_comm, &qubit_size);

    // //loop over all times
    // for (int i = 0; i < times.size(); i++) { 
    //     t = times[i];
    t = 192; //TODO: create parallelizem here

            
        // //print time
        // std::cout << std::endl << "Main-Time: " << t << std::endl;

    //vector of all transition amplitude for a given time
    vector<complex<double>> transition_amplitudes(target_max);

    int target_group = qubit_rank / target_max;
    MPI_Comm target_comm;
    MPI_Comm_split(qubit_comm, target_group, world_rank, &target_comm);

    int target_rank, target_size;
    MPI_Comm_rank(target_comm, &target_rank);
    MPI_Comm_size(target_comm, &target_size);


    int q_group = target_group / q_max;
    MPI_Comm q_comm;
    MPI_Comm_split(target_comm, q_group, world_rank, &q_comm);
    
    int q_rank, q_size;
    MPI_Comm_rank(q_comm, &q_rank);
    MPI_Comm_size(q_comm, &q_size);


    complex<double> total_amp;
    complex<double> q_amp = mpi_compute_q_amplitude(target_rank, q_rank, qubit_rank, amplitude, t, chi, omega, total_amp, all_coefficients);
    
    MPI_Barrier(q_comm);

    MPI_Reduce(&q_amp, &total_amp, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, q_comm);
    // if(q_rank == 0){
    //     transition_amplitudes[target_rank] = total_amp;
    // }
    
    MPI_Barrier(target_comm);

    // for(int i=0; i < target_group; i++){
    //     if(target_rank == 0){
    //         //print transition amplitudes
    //         std::cout << "Transition amplitudes: " << std::endl;
    //     }
    //     if(target_rank == i){
    //         std::cout << target_rank << ": " << total_amp << std::endl;

    //     }
    //     // Ensure the processes proceed in order
    //     MPI_Barrier(target_comm);
    // }

    if(target_rank != 0){
        MPI_Send(&total_amp, 1, MPI_DOUBLE_COMPLEX, 0, 0, target_comm);
        // MPI_Gather(&total_amp, 1, MPI_DOUBLE_COMPLEX, nullptr, 0, MPI_DOUBLE_COMPLEX, 0, target_comm)
    }else{
        transition_amplitudes[0] = total_amp;
        for(int i=1; i < target_group; i++){
            MPI_Recv(&transition_amplitudes[i], 1, MPI_DOUBLE_COMPLEX, i, 0, target_comm, MPI_STATUS_IGNURE);
        }
        // MPI_Gather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, results.data(), 1, MPI_INT, split_rank, split_comm);
    }

    MPI_Barrier(target_comm);

    if(target_rank == 0){
        // normalize transition amplitude 
        double norm = 0;
        for (const auto& amp : transition_amplitudes) {
            norm += std::norm(amp);//TODO: check if norm works with complex<double>
        }
        norm = sqrt(norm);
        for (auto& amp : transition_amplitudes) {
            amp = amp / norm;
        }

        // all_transition_amplitudes[qubit_rank] = transition_amplitudes
        // all_transition_amplitudes[qubit_rank][t]

        // Assuming all_transition_amplitudes is a vector<vector<complex<double>>>
        nlohmann::json j;
        for (const auto& c : transition_amplitudes) {
            j.push_back({{"real", c.real()}, {"imag", c.imag()}});
        }
        if (qubit_rank == 0){
            std::string filename = "amplitudes_g_chg_amp" + std::to_string(amplitude) + "_q" + std::to_string(q_max) + ".json";
        }else{
            std::string filename = "amplitudes_e_chg_amp" + std::to_string(amplitude) + "_q" + std::to_string(q_max) + ".json";
        }
        
        std::ofstream o(filename);
        o << j << std::endl;
    }

    //TODO: add support for more than a single t
    

    







    

    
            //loop over all target steps
            for (int target = 0; target <= MAX_TARGET; target++) {
                complex<double> total_amp;

                //print target
                std::cout << "Target: " << target << std::endl;
                
                
                mpi_compute_q_amplitude(target, q_max, qubit, amplitude, t, chi, omega, total_amp, all_coefficients);
                
                // Synchronize all processes before moving to the next target
                MPI_Barrier(MPI_COMM_WORLD);

                // // Print the Total amplitude (only from the root process)
                // if (q_rank == 0) {
                //     cout << "Total amplitude for target " << target << ": " << total_amp << endl << endl;
                // }
                

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
                norm += std::norm(amp);//TODO: check if norm works with complex<double>
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

            // Synchronize all processes before moving to the next target
            MPI_Barrier(MPI_COMM_WORLD);

        //end of time loop
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

        // Synchronize all processes before moving to the next target
        MPI_Barrier(MPI_COMM_WORLD);
    
    // } //end of qubit loop

    MPI_Finalize();


    return 0;
}





