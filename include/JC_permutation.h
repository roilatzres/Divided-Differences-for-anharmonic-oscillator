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




// Function to reverse the permutation
std::vector<std::vector<std::pair<int, int>>> revesre_permutaitons(std::vector<std::vector<std::pair<int, int>>> all_combinations) {
    
    std::vector<std::vector<std::pair<int, int>>> all_permutation;
    
    for (int i = 0; i < all_combinations.size(); i++) {
        std::vector<std::pair<int, int>> permutation;
        for (int j = 0; j < all_combinations[i].size(); j++) {
            permutation.push_back({-all_combinations[i][j].first, -all_combinations[i][j].second});
        }
        all_permutation.push_back(permutation);
    }
    return all_permutation;
}


void backtrack(std::pair<int,int> curr_step, std::pair<int,int> target_step, int remaining_moves, std::pair<int, int> max_step, std::vector<std::pair<int, int>>& combination, std::vector<std::vector<std::pair<int, int>>>& all_comb) {
//print current step
        std::cout << "Current step: " << curr_step.first << " " << curr_step.second << std::endl;
        //print remaining moves
        std::cout << "Remaining moves: " << remaining_moves << std::endl;
        //print target step
        std::cout << "Target step: " << target_step.first << " " << target_step.second << std::endl;

        std::cout << "Combination: ";
        for (std::pair<int, int> step : combination) {
            std::cout << "{" << step.first << ", " << step.second << "} ";
        }

        std::cout << std::endl << std::endl;

        // Base case: If current step is target step, add the combination to the list
        if (curr_step.first == target_step.first && curr_step.second == target_step.second && remaining_moves == 0) {
            all_comb.push_back(combination);
            //print combination
            std::cout << "Combination successfully added to all_comb" << std::endl << std::endl ;
            return;
        }

        // Base case: If remaining moves is 0 but current step is not target step, return
        if (remaining_moves == 0) {
            return;
        }

        // Try moving (1, 0) one step up
        if(curr_step.first < max_step.first){
            combination.push_back({1, 0});
            curr_step.first += 1;
            backtrack(curr_step, target_step, remaining_moves - 1, max_step, combination, all_comb);
            curr_step.first -= 1;
            combination.pop_back();
        }

        // Try moving (-1, 0) one step
        if (curr_step.first > 0) {
            combination.push_back({-1, 0});
            curr_step.first -= 1;
            backtrack(curr_step, target_step, remaining_moves - 1, max_step, combination, all_comb);
            curr_step.first += 1;
            combination.pop_back();
        }

        
        // Try moving (1, -1) one step up
        //print attempt
        std::cout << "try (1, -1)" << std::endl;
        if(curr_step.first < max_step.first && curr_step.second > 0){//TODO: add print
            combination.push_back({1, -1});
            curr_step.first += 1;
            curr_step.second -= 1;
            backtrack(curr_step, target_step, remaining_moves - 1, max_step, combination, all_comb);
            curr_step.first -= 1;
            curr_step.second += 1;
            combination.pop_back();
        }

        // Try moving (-1, 1) one step down
        //print attempt
        std::cout << "try (-1, 1)" << std::endl;
        if (curr_step.first > 0 && curr_step.second < max_step.second) {
            combination.push_back({-1, 1});
            curr_step.first -= 1;
            curr_step.second += 1;
            backtrack(curr_step, target_step, remaining_moves - 1, max_step, combination, all_comb);
            curr_step.first += 1;
            curr_step.second -= 1;
            combination.pop_back();
        }
}


// Function to compute all possible permutations of steps to reach the target step on a ladder
std::tuple<std::vector<std::vector<std::pair<int, int>>>, std::vector<std::vector<std::pair<int, int>>>> ladder_permutations(std::pair<int, int> start_step, std::pair<int, int> mid_step, std::pair<int, int> target_step, std::pair<int, int> max_step, int steps) {
    std::vector<std::vector<std::pair<int, int>>> all_combinations_up;
    std::vector<std::vector<std::pair<int, int>>> all_combinations_down;


    // // Backtrack function to generate permutations recursively
    // std::function<void(std::pair<int,int> , std::pair<int,int>, int, std::vector<std::pair<int, int>>&, std::vector<std::vector<std::pair<int, int>>>&)> backtrack = 
    // [&](std::pair<int,int> curr_step, std::pair<int,int> target_step, int remaining_moves, std::vector<std::pair<int, int>>& combination, std::vector<std::vector<std::pair<int, int>>>& all_comb) {

    //     //print current step
    //     std::cout << "Current step: " << curr_step.first << " " << curr_step.second << std::endl;
    //     //print remaining moves
    //     std::cout << "Remaining moves: " << remaining_moves << std::endl;

    //     // Base case: If current step is target step, add the combination to the list
    //     if (curr_step.first == target_step && remaining_moves == 0) {
    //         all_comb.push_back(combination);
    //         //print combination
    //         std::cout << "Combination: ";
    //         for (std::pair<int, int> step : combination) {
    //             std::cout << "{" << step.first << ", " << step.second << "} ";
    //         }
    //         return;
    //     }

    //     // Base case: If remaining moves is 0 but current step is not target step, return
    //     if (remaining_moves == 0) {
    //         return;
    //     }

    //     // Try moving (1, 0) one step up
    //     if(curr_step.first < max_step.first){
    //         combination.push_back({1, 0});
    //         curr_step.first += 1;
    //         backtrack(curr_step, target_step, remaining_moves - 1, combination, all_comb);
    //         combination.pop_back();
    //     }

    //     // Try moving (-1, 0) one step
    //     if (curr_step.first > 0) {
    //         combination.push_back({-1, 0});
    //         curr_step.first -= 1;
    //         backtrack(curr_step, target_step, remaining_moves - 1, combination, all_comb);
    //         combination.pop_back();
    //     }

        
    //     // Try moving (1, -1) one step up
    //     //print attempt
    //     if(curr_step.first < max_step.first && curr_step.second > 0){//TODO: add print
    //         combination.push_back({1, -1});
    //         curr_step.first += 1;
    //         curr_step.second -= 1;
    //         backtrack(curr_step, target_step, remaining_moves - 1, combination, all_comb);
    //         combination.pop_back();
    //     }

    //     // Try moving (-1, 1) one step down
    //     if (curr_step.first > 0 && curr_step.second < max_step.second) {
    //         combination.push_back({-1, 1});
    //         curr_step.first -= 1;
    //         curr_step.second += 1;
    //         backtrack(curr_step, target_step, remaining_moves - 1, combination, all_comb);
    //         combination.pop_back();
    //     }
    // };

    // Initialize combination list
    std::vector<std::pair<int, int>> combinations;

    // Start backtracking from start_step with given number of moves
    backtrack(start_step, mid_step, steps, max_step, combinations, all_combinations_up);
    
    // if target_step equals start_step, reverse the permutation
    if (target_step == start_step &&  mid_step == max_step) {
        all_combinations_down = revesre_permutaitons(all_combinations_up);
    }else{
        combinations.clear();
        // Start backtracking from target to start step with given number of moves
        backtrack(mid_step, target_step, steps, max_step, combinations, all_combinations_down);
    }
    
    return std::make_tuple(all_combinations_up, all_combinations_down);
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

