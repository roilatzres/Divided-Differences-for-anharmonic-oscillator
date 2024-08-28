#ifndef OPENMP_UTILS_H
#define OPENMP_UTILS_H

// #include <mpi.h>
#include <omp.h>
#include <iostream>
#include <vector>
#include <complex>
#include <fstream>
#include <string>

using std::vector;

//TODO: check up on all_coeefiicient refrencing

complex<double> mpi_compute_q_amplitude(int target, int q, int qubit, double amplitude, double t, double chi, double omega, 
                            vector<vector<vector<std::tuple<vector<DivdiffElement>, double, int, vector<int>>>>> &all_coefficients) {

    complex<double> q_amp = 0;
    complex<double> energy_coefficient_amp = cal_neg_i_pow(q) * std::pow(amplitude, q);

    if (target == 0 && q == 0) {
        return 1;
    }
        

    // std::cout << "Permutations:" << std::endl;
    // for (const auto& combination : permutations) {
    //     for (int step : combination) {
    //         std::cout << step << " ";
    //     }
    //     std::cout << std::endl;
    // }
    auto coefficients = all_coefficients[target][q];
    
    //there are no permutaions - skip
    if (coefficients.size() == 0) {
        return 0;
    }
    
    // Check if file for q_amp already exists
    std::string q_amp_filename = "./include/parms/dispersive/q_amp_recursive/dispersive_q_amp_qubit_" 
        + std::to_string(qubit) + "_t_" + std::to_string((int)t) + "_target_" + std::to_string(target) + "_q" + std::to_string(q) + ".bin";
    
    std::ifstream infile(q_amp_filename);
    if (infile) {
        // Load q_amp from file
        q_amp = load_complex_from_file(q_amp_filename);
        q_amp = q_amp * energy_coefficient_amp;
        
        #pragma omp critical
        std::cout  << "qubit: " << qubit << ", target: " << target << ", " << q << " Original q_amp: " << q_amp << std::endl;

        
    }else{// Calculate q_amp
        
        divdiff_init();
        for (const auto& coefficient : coefficients) {
            q_amp += cal_divdiff(coefficient, t, q, qubit, chi, omega, target);
        }
        divdiff_clear_up();

        // Save q_amp to file
        save_complex_to_file(q_amp, q_amp_filename);
        q_amp = q_amp * energy_coefficient_amp;

        // // Print the q amplitude
        // std::cout << q << " Original q_amp: " << q_amp << std::endl;

    }
 
    return q_amp;
}


// overload of reduction operator
struct complex_double_reduction {
    typedef std::complex<double> value_type;
    
    static void init(value_type& val) {
        val = value_type(0, 0);
    }
    
    static void combine(value_type& dest, const value_type& src) {
        dest += src;
    }
    
    static void final(value_type& val) {
        // No action needed
    }
};


#endif // OPENMP_UTILS_H