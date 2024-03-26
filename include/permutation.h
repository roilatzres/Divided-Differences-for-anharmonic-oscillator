#include <iostream>
#include <vector>
#include <functional>
#include <cmath>
#include <complex>



double factorial(unsigned int n) {
    if (n == 0)
        return 1;
    else
        return n * factorial(n - 1);
}


// ****************************************************************************************
// ***************** permutation generation ***********************************************
// ****************************************************************************************


// Function to compute all possible permutations of steps to reach the target step on a ladder
std::vector<std::vector<int>> ladder_permutations(int total_steps, int target_step, int moves) {
    std::vector<std::vector<int>> all_combinations;

    // Backtrack function to generate permutations recursively
    std::function<void(int, int, std::vector<int>&)> backtrack = [&](int curr_step, int remaining_moves, std::vector<int>& combination) {
        // Base case: If remaining moves is 0 and current step is target step, add current combination
        if (remaining_moves == 0 && curr_step == target_step) {
            all_combinations.push_back(combination);
            return;
        }

        // Base case: If remaining moves is 0 but current step is not target step, return
        if (remaining_moves == 0 || curr_step >= total_steps) {
            return;
        }

        // Try moving up one step
        combination.push_back(1);
        backtrack(curr_step + 1, remaining_moves - 1, combination);
        combination.pop_back();

        // Try moving down one step
        if (curr_step > 0) {
            combination.push_back(-1);
            backtrack(curr_step - 1, remaining_moves - 1, combination);
            combination.pop_back();
        }
    };

    // Initialize combination list
    std::vector<int> combination;
    // Start backtracking from step 0 with given number of moves
    backtrack(0, moves, combination);

    return all_combinations;
}



// ***************************************************************************************
// ***************** coefficient calculation *********************************************
// ***************************************************************************************

// Struct to represent divdiff element
struct DivdiffElement {
    int energy;
    int omega;
    int chi;
    int coefficient; // either 1 or -1
};

// Calculate i^power
std::complex<double> cal_i_pow(int power) {
    if(power%4 == 0){
        return 1;
    }
    else if(power%4 == 1){
        return 1j;
    }
    else if(power%4 == 2){
        return -1;
    }
    else{//power%4 == 3
        return -1j;
    }
}

// Calculate -i^power
std::complex<double> cal_neg_i_pow(int power) {
    if(power%4 == 0){
        return 1;
    }
    else if(power%4 == 1){
        return -1j;
    }
    else if(power%4 == 2){
        return -1;
    }
    else{//power%4 == 3
        return 1j;
    }
}



// Function to compute coefficients and divdiff elements for a given permutation
std::tuple<std::vector<std::vector<DivdiffElement>>, std::complex<double>, int> cal_coefficient(std::vector<int> permutation, int start_state) {
    // Initialize divdiff elements, current state, and coefficient
    double energy_coefficient_real = 1.0;
    int current_state = start_state;

    // initialize vector for all states in the permutation
    std::vector<int> states;
    states.push_back(start_state);

    // Iterate through each step in the permutation
    for (int i = 0; i < permutation.size(); ++i) {
        // Update current state
        int new_state = current_state + permutation[i];
        states.push_back(new_state);
        
        // Update energy_coefficient_real
        if (permutation[i] == 1) {
            energy_coefficient_real = sqrt(current_state + 1) * energy_coefficient_real;
        } else {
            energy_coefficient_real = sqrt(current_state) * energy_coefficient_real;
        }
        current_state = new_state;
    }

    std::complex<double> energy_coefficient_imag = cal_neg_i_pow(permutation.size()) * std::pow(0.0005, permutation.size());
    std::complex<double> energy_coefficient = energy_coefficient_imag * energy_coefficient_real;
    int final_state = current_state;

    // Initialize beta list with the final state
    std::vector<std::vector<DivdiffElement>> beta = {{{states.back(), 0, 0, 1}}};

    // iterate from end to start over the states and update the divdiff elements starting from the penultimate state
    for (int i = states.size() - 2; i >= 0; --i) {
        std::vector<std::vector<DivdiffElement>> new_divdiff;

        // Update divdiff elements
        for (int j = 0; j < beta.size(); ++j) {
            // Extract the first element from beta[j]
            DivdiffElement first_element = beta[j].front();

            // prep new divdiff elements 
            int p_mat = states[i] - states[i + 1]; // going backwards on states, checking previous state
            
            // create new divdiff elements according to the permutation matrix
            int chi_factor;
            if(p_mat == 1){
                chi_factor = 1;
            }
            else{//p_mat == -1
                chi_factor = -1;
                // DivdiffElement plus = {states[i], first_element.omega + 1, first_element.chi + 1, first_element.coefficient * (1i / 2)};
                // DivdiffElement minus = {states[i], first_element.omega - 1, first_element.chi + 1, first_element.coefficient * (-1i / 2)};
            }
            
            DivdiffElement plus = {states[i], first_element.omega + 1, first_element.chi - chi_factor, (-first_element.coefficient)};
            DivdiffElement minus = {states[i], first_element.omega - 1, first_element.chi - chi_factor, first_element.coefficient};

            // Copy the list and insert new elements in the front
            std::vector<DivdiffElement> new_divdiff_plus = beta[j];
            new_divdiff_plus.insert(new_divdiff_plus.begin(), plus);

            std::vector<DivdiffElement> new_divdiff_minus = beta[j];
            new_divdiff_minus.insert(new_divdiff_minus.begin(), minus);
            

            // Append new lists to new_divdiff
            new_divdiff.push_back(new_divdiff_plus);
            new_divdiff.push_back(new_divdiff_minus);
        }

        beta = new_divdiff;
    }

    return std::make_tuple(beta, energy_coefficient, final_state);
}

