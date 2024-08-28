#include <iostream>
#include <vector>
#include <functional>
#include <cmath>
#include <complex>
#include "include\divdiffcomplex.h"
#include "include\permutation.h"
#include <fstream>
#include "include\json-develop\single_include\nlohmann\json.hpp"



int main() {
    //original parameters
    double multiplier = 4;
    float sigma = 48;
    float timestep = 1;
    double chi = -279e-6 * 2 * M_PI; 
    // double detuning = chi *1e-2;

    // //print vars
    // std::cout << "Multiplier: " << multiplier << std::endl;
    // std::cout << "Sigma: " << sigma << std::endl;
    // std::cout << "Timestep: " << timestep << std::endl;
    // std::cout << "Chi: " << chi << std::endl;
    // // std::cout << "Detuning: " << detuning << std::endl;


    double tfinal = multiplier * sigma;
    double omega = (2 * M_PI) / tfinal;


    int t = 192;

    for(int q = 0; q < 8; q++){

        std::complex<double> q_amp = 0;
        int max_step = 10;
        //ladder_permutations(max_step, target, Q = total moves)
        auto permutations = ladder_permutations(max_step, 0, q);
        
        //there are no permutaions - skip
        if (permutations.size() == 0) {
            if (target == 0 && q == 0){
                q_amp += 1;
            }
            continue;
        }

        //print all permutations
        std::cout << "Possible permutations:" << std::endl;
        for (const auto& combination : permutations) {
            for (int step : combination) {
                std::cout << step << " ";
            }
            std::cout << std::endl;
        }

        // loop over all permutations
        for (const auto& permutation : permutations) {
            std::complex<double> perm_amp = 0;

        }
    }
}