#include <iostream>
#include <vector>
#include <functional>
#include <cmath>


// **********************************************************
// ***************** permutation generation *****************
// **********************************************************


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



// ***********************************************************
// ***************** coefficient calculation *****************
// ***********************************************************

// Struct to represent divdiff element
struct DivdiffElement {
    int energy;
    int omega;
};

// Function to compute coefficients and divdiff elements for a given permutation
std::tuple<std::vector<std::vector<DivdiffElement>>, double, int> cal_coefficient(std::vector<int> permutation, int start_state) {
    // Initialize divdiff elements, current state, and coefficient
    std::vector<std::vector<DivdiffElement>> beta = {{{start_state, 0}}};
    double coefficient = 1.0;
    int current_state = start_state;

    // Iterate through each step in the permutation
    for (int i = 0; i < permutation.size(); ++i) {
        std::vector<std::vector<DivdiffElement>> new_divdiff;

        // Update current state
        int new_state = current_state + permutation[i];

        // Update coefficient
        if (permutation[i] == 1) {
            coefficient = sqrt(current_state + 1) * coefficient;
        } else {
            coefficient = sqrt(current_state) * coefficient;
        }

        // Update divdiff elements
        for (int j = 0; j < beta.size(); ++j) {
            // Extract the last element from beta[j]
            DivdiffElement last_element = beta[j].back();

            // Create new divdiff elements
            DivdiffElement plus = {last_element.energy + permutation[i], last_element.omega + 1};
            DivdiffElement minus = {last_element.energy + permutation[i], last_element.omega - 1};

            // Copy the list and append new elements
            std::vector<DivdiffElement> new_divdiff_plus = beta[j];
            new_divdiff_plus.push_back(plus);

            std::vector<DivdiffElement> new_divdiff_minus = beta[j];
            new_divdiff_minus.push_back(minus);

            // Append new lists to new_divdiff
            new_divdiff.push_back(new_divdiff_plus);
            new_divdiff.push_back(new_divdiff_minus);
        }

        beta = new_divdiff;
        current_state = new_state;
    }

    return std::make_tuple(beta, coefficient, current_state);
}

