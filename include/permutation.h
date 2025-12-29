#ifndef PERMUTATION_H
#define PERMUTATION_H

#include <iostream>
#include <vector>
#include <functional>
#include <cmath>
#include <complex>
// #include "dlib/dlib/optimization.h"
// #include "dlib/dlib/matrix.h"


using std::vector;
using std::complex;
using std::cout;
using std::endl;

//struct consisting of 2 ExExFloats
struct complex_Ex{
    ExExFloat real;
    ExExFloat imag;
    
    complex_Ex(ExExFloat r = 0, ExExFloat i = 0) {
        real = r;
        imag = i;
    }

};

//function to print complex_ex
void print_complex_Ex(complex_Ex a){
    a.real.print();
    std::cout << "+ ";
    a.imag.print();
    std::cout << "i" << std::endl; 
    std::cout << std::endl;
}

complex_Ex complex_mult(complex_Ex a, complex_Ex b){
    complex_Ex res;
    ExExFloat temp1 = a.real * b.real;
    ExExFloat temp2 = a.imag * b.imag;
    res.real = temp1 - temp2;

    ExExFloat temp3 = a.real * b.imag;
    ExExFloat temp4 = a.imag * b.real;
    res.imag = temp3 + temp4;


    return res;
}

complex_Ex complex_mult(complex_Ex a, complex<double> b){
    divdiff_init();
    complex_Ex res;
    ExExFloat temp1 = a.real * b.real();
    ExExFloat temp2 = a.imag * b.imag();
    res.real = temp1 - temp2;

    ExExFloat temp3 = a.real * b.imag();
    ExExFloat temp4 = a.imag * b.real();
    res.imag = temp3 + temp4;


    return res;
}

double factorial(unsigned int n) {
    if (n == 0)
        return 1;
    else
        return n * factorial(n - 1);
}

//factorial for ExExFloat
ExExFloat ExFactorial(ExExFloat base, int n) {
    if (n == 1)
        return (ExExFloat) 1;
    else{
        // //print base and n 
        // std::cout << "base: " << std::endl;
        // base.print();
        // std::cout << std::endl;
        // std::cout << "n: " << n << std::endl;
        // ExExFloat one = 1;
        // //print one
        // std::cout << "one: " << std::endl;
        // one.print();
        // std::cout << std::endl;

        // ExExFloat new_base;
        // cout << "base: " << endl;
        // new_base = base - one;
        // //print new_base
        // std::cout << "new_base: " << std::endl;
        // new_base.print();
        // std::cout << std::endl;
        
        //calculate result
        ExExFloat result = base * ExFactorial(base - 1, n - 1);
        // // //print result
        // std::cout << "result: "  << std::endl;
        // result.print();
        // std::cout << std::endl;
        return result;
    }
}

//power for ExExFloat
ExExFloat ExPow(ExExFloat base, int exp) {
    if (exp == 0)
        return 1;
    else{
        ExExFloat result = base * ExPow(base, exp - 1);
        // //print result
        // std::cout << "result: "  << std::endl;
        // result.print(); 
        return result;
    }
}

// ****************************************************************************************
// ***************** permutation generation ***********************************************
// ****************************************************************************************
struct step {
    int cavity;
    int qubit;
};

// Function to compute all possible permutations of steps to reach the target step on a ladder
vector<vector<step>> jc_ladder_permutations(step start_step, step top_step, step target, int moves) {
    vector<vector<step>> all_combinations;

    // Backtrack function to generate permutations recursively
    std::function<void(step, int, vector<step>&)> jc_backtrack = [&](step curr_step, int remaining_moves, vector<step>& combination) {
        
        // // Print the current path
        // std::cout << "Backtracking: ";
        // std::cout << "curr step: (" << curr_step.cavity << "," << curr_step.qubit << ") |" << endl;
        // for (const auto& s : combination) {
        //     std::cout << " (" << s.cavity << "," << s.qubit << ") -> ";
        // }
        // std::cout << "| Moves left: " << remaining_moves << std::endl;

        // Base case: If remaining moves is 0 and current step is target step, add current combination
        if (remaining_moves == 0) {
            if (curr_step.cavity == target.cavity && curr_step.qubit == target.qubit) {
                // std::cout << "Found valid permutation: ";
                all_combinations.push_back(combination);
            }
            return;
        }


        // Try moving cavity + 1, qubit
        if( curr_step.cavity < top_step.cavity) {
            combination.push_back({curr_step.cavity + 1, curr_step.qubit});
            jc_backtrack({curr_step.cavity + 1, curr_step.qubit}, remaining_moves - 1, combination);
            combination.pop_back();
        }
        // Try moving cavity - 1, qubit
        if( curr_step.cavity > 0) {
            combination.push_back({curr_step.cavity - 1, curr_step.qubit});
            jc_backtrack({curr_step.cavity - 1, curr_step.qubit}, remaining_moves - 1, combination);
            combination.pop_back();
        }

        // Try moving cavity + 1, qubit - 1
        if( curr_step.qubit > 0 && curr_step.cavity < top_step.cavity) {
            combination.push_back({curr_step.cavity + 1, curr_step.qubit - 1});
            jc_backtrack({curr_step.cavity + 1, curr_step.qubit - 1}, remaining_moves - 1, combination);
            combination.pop_back();
        }   

        // Try moving cavity - 1, qubit + 1
        if( curr_step.cavity > 0 && curr_step.qubit < top_step.qubit) {
            combination.push_back({curr_step.cavity - 1, curr_step.qubit + 1});
            jc_backtrack({curr_step.cavity - 1, curr_step.qubit + 1}, remaining_moves - 1, combination);
            combination.pop_back();
        }   
    };

    // Initialize combination list
    vector<step> combination;
    // Start backtracking start_step with given number of moves
    jc_backtrack(start_step, moves, combination);

    return all_combinations;
}

// Function to compute all possible permutations of steps to reach the target step on a ladder
vector<vector<int>> ladder_permutations(int start_step, int top_step, int target_step, int moves) {
    vector<vector<int>> all_combinations;

    // Backtrack function to generate permutations recursively
    std::function<void(int, int, vector<int>&)> backtrack = [&](int curr_step, int remaining_moves, vector<int>& combination) {
        // Base case: If remaining moves is 0 and current step is target step, add current combination
        if (remaining_moves == 0 && curr_step == target_step) {
            all_combinations.push_back(combination);
            return;
        }

        // Base case: If remaining moves is 0 but current step is not target step, return
        if (remaining_moves == 0 || curr_step > top_step) {
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
    vector<int> combination;
    // Start backtracking from the given start step with given number of moves
    backtrack(start_step, moves, combination);

    return all_combinations;
}



// ***************************************************************************************
// ***************** coefficient calculation *********************************************
// ***************************************************************************************

// Struct to represent divdiff element
struct DivdiffElement {
    vector<int> omega;
    int coefficient; // either 1 or -1
};

// Calculate i^power
complex<double> cal_i_pow(int power) {
    complex<double> j(0, 1);

    if(power%4 == 0){
        return 1;
    }
    else if(power%4 == 1){
        return j;
    }
    else if(power%4 == 2){
        return -1;
    }
    else{//power%4 == 3
        return -j;
    }
}

// Calculate -i^power
complex<double> cal_neg_i_pow(int power) {
    complex<double> j(0, 1);
    if(power%4 == 0){
        return 1;
    }
    else if(power%4 == 1){
        return -j;
    }
    else if(power%4 == 2){
        return -1;
    }
    else{//power%4 == 3
        return j;
    }
}



// Function to compute coefficients and divdiff elements for a given permutation
std::tuple<vector<DivdiffElement>, double, int, vector<int>> cal_coefficient(vector<int> permutation, int start_state) {
    // Initialize divdiff elements, current state, and coefficient
    double energy_coefficient_real = 1;
    int current_state = start_state;

    // initialize vector for all states in the permutation
    vector<int> states;
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


    int final_state = current_state;

    // Initialize beta list with the final state
    DivdiffElement first_elem;
    first_elem.omega = {0};
    first_elem.coefficient = 1;
    vector<DivdiffElement> beta;
    beta.push_back(first_elem); 

    // iterate from end to start over the states and update the divdiff elements starting from the penultimate state
    for (int i = states.size() - 2; i >= 0; --i) {//TODO: check direction of omega sum
        vector<DivdiffElement> temp_divdiff;

        // Update divdiff elements
        for (int j = 0; j < beta.size(); ++j) {
            // Extract the first element from beta[j]
            DivdiffElement beta_element = beta[j];
            
            DivdiffElement plus, minus;
            plus.coefficient = -beta_element.coefficient;
            plus.omega = beta_element.omega; 
            
            //access last element of omega and add 1
            int last_element_plus = beta_element.omega.back();
            plus.omega.push_back(last_element_plus + 1);

            minus.coefficient = beta_element.coefficient;
            minus.omega = beta_element.omega;

            //access last element of omega and add 1
            int last_element_minus = beta_element.omega.back();
            minus.omega.push_back(last_element_minus - 1);

            //insert new elements to temp_divdiff
            temp_divdiff.push_back(plus);
            temp_divdiff.push_back(minus);

        }

        beta = temp_divdiff;
    }

    return std::make_tuple(beta, energy_coefficient_real, final_state, states);
}



// Function to compute coefficients and divdiff elements for a given permutation
std::tuple<double, int, vector<int>> cal_coefficient_const_pulse(vector<int> permutation, int start_state) {
    // Initialize divdiff elements, current state, and coefficient
    double energy_coefficient_real = 1;
    int current_state = start_state;

    // initialize vector for all states in the permutation
    vector<int> states;
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

    
    int final_state = current_state;
    return std::make_tuple(energy_coefficient_real, final_state, states);
    
}

// Function to compute coefficients and divdiff elements for a given permutation with the jc hamiltonaian
std::tuple<double, step, vector<step>, step> cal_coefficient_const_pulse_jc(vector<step> permutation, step start_state, double amp, double coupling) {
    // Initialize divdiff elements, current state, and coefficient
    double energy_coefficient_real = 1;
    step current_state = start_state;
    step count = {0,0};
    
    // initialize vector for all states in the permutation
    vector<step> states;
    states.push_back(start_state);

    for(int i = 0; i < permutation.size(); ++i) {
                
        // Update energy_coefficient_real
        if (permutation[i].cavity - current_state.cavity == 1 && permutation[i].qubit == current_state.qubit) { // a^dagger
            energy_coefficient_real *= sqrt(current_state.cavity + 1);
            count.cavity += 1;
        }
        else if (permutation[i].cavity - current_state.cavity == -1 && permutation[i].qubit == current_state.qubit) { // a
            energy_coefficient_real *= sqrt(current_state.cavity);
            count.cavity += 1;
        }
        else if (permutation[i].cavity - current_state.cavity == 1 && permutation[i].qubit - current_state.qubit == -1) { // a^dagger q
            energy_coefficient_real *= sqrt(current_state.qubit * (current_state.cavity + 1)) ;
            count.qubit += 1;
        } 
        else if (permutation[i].cavity - current_state.cavity == -1 && permutation[i].qubit - current_state.qubit == 1) { // a q^dagger
            energy_coefficient_real *= sqrt((current_state.qubit + 1) * current_state.cavity) ;
            count.qubit += 1;
        }
        
        states.push_back(permutation[i]);
        current_state = permutation[i];
    }

    return std::make_tuple(energy_coefficient_real, current_state, states, count);

}



/************** calculate divided differences for specific q and t****************/

void cal_divdiff_aux_td(divdiffcomplex &d, const vector<DivdiffElement> &divdiff, const vector<int> &states, 
    int *CurrentLength, int *num_elem, double t, int q, int qubit, double chi, double omega, int t_k, complex_Ex &sum_final_elements){
    complex<double> j(0, 1);
    
    //first state
    if (*CurrentLength == 0){
        cout << "First state. CurrentLength: " << *CurrentLength << endl;
        cout << "num_elem: " << *num_elem << endl;  
        cout << "states size: " << states.size() << endl;

        const auto curr_elem = divdiff[*num_elem];
        
        const double& beta_coefficient = curr_elem.coefficient;
        const auto& element_omega = curr_elem.omega[*CurrentLength];
        const auto& state_n = states[*CurrentLength];

        // std::cout << "(" << element.energy << ", " << element.omega << ", " << element.chi << ", " << element.coefficient << ") ";

        complex<double> n;
        n = (-j * t * (state_n * chi * (qubit) + element_omega * omega ));
        
        // print n
        std::cout << "n: " << n << std::endl;
        cout << "Adding element to divdiffcomplex" << std::endl;
        cout << "CurrentLength before adding: " << d.CurrentLength << std::endl;
        
        d.AddElement(n);
        *CurrentLength += 1;

        cout << "CurrentLength after adding: " << d.CurrentLength << endl;

        cal_divdiff_aux_td(d, divdiff, states, CurrentLength, num_elem, t, q, qubit, chi, omega, t_k, sum_final_elements);
        return;
    }
    
   
    for(int i=0; i < 2; i++){
        const auto curr_elem = divdiff[*num_elem];
            
        const double& beta_coefficient = curr_elem.coefficient;
        const auto& element_omega = curr_elem.omega[*CurrentLength];
        const auto& state_n = states[*CurrentLength];
        cout << "num_elem: " << *num_elem << endl;
        cout << "CurrentLength: " << *CurrentLength << endl;
        cout << "states size: " << states.size() << endl;
        cout << "element_omega: " << element_omega << endl;
        cout << "state_n: " << state_n << endl;
        cout << "beta_coefficient: " << beta_coefficient << endl;
        cout << "qubit: " << qubit << endl;
        cout << "chi: " << chi << endl;

        complex<double> n;
        n = (-j * t * (state_n * chi * (qubit) + element_omega * omega ));
        
        // print n
        std::cout << "n: " << n << std::endl;
        cout << n.imag() << endl;  
        
        d.AddElement(n);
        *CurrentLength += 1;

        //final state
        if(*CurrentLength == (states.size())){
            cout << "Final state reached. CurrentLength: " << *CurrentLength << endl;
            //cal z_n 
            complex<ExExFloat> z_n = d.divdiffs[*CurrentLength - 1];
            // z_n.real().print();
            // cout << endl;
            // z_n.imag().print();
            // cout << endl;
            // double real = z_n.real().get_double();
            // double imag = z_n.imag().get_double();
            //print real and imag
            // std::cout << "Divdiff = " << real << " + " << imag << "i" << std::endl;
            // //print divdiffs elements
            // d.PrintList(d.divdiffs, d.CurrentLength, "Divdiffs");
            
            complex<ExExFloat> final_element = (z_n * cal_neg_i_pow(q) * pow(t, q)) /factorial(q) * beta_coefficient;//TODO: switch to complex_Ex?
            
            complex_Ex start_time_phase = complex_Ex(cos(omega * element_omega * t_k), sin(omega * element_omega * t_k));
            //print start_time_phase
            cout << "Start time phase: ";
            print_complex_Ex(start_time_phase);
            complex_Ex res = complex_mult(complex_Ex(final_element.real(), final_element.imag()), start_time_phase);
            cout << "Resulting element: ";
            print_complex_Ex(res);
            // complex<ExExFloat> norm_real = z_n.real() * cal_neg_i_pow(q) * pow(t, q) / factorial(q) ;
            // complex<ExExFloat> norm_imag = z_n.imag() * (-1) * cal_neg_i_pow(q+1) * pow(t, q) / factorial(q);// (-it)^q / q! * Imag = (-it)^q / q! * (i) * Imag.real = (-i)^q+1 * t^q * -Imag.real / q!
            // complex<ExExFloat> final_element = (norm_real + norm_imag) * beta_coefficient;

            //add the divdiff element to the total sum
            // *sum_final_elements = *sum_final_elements + res; 
            sum_final_elements.real += res.real;
            sum_final_elements.imag += res.imag;
            // print sum_final_elements
            cout << "Sum final elements so far: ";
            print_complex_Ex(sum_final_elements);
            cout << endl;
            *num_elem += 1;

        }else{// middle states
            cout << "Middle state. CurrentLength: " << *CurrentLength << endl;
            cal_divdiff_aux_td(d, divdiff, states, CurrentLength, num_elem, t, q, qubit, chi, omega, t_k, sum_final_elements);
        }

        d.RemoveElement();
        *CurrentLength -= 1;
    }
    return;
}

// set divdiff calculation
complex_Ex cal_divdiff_td(const std::tuple<std::vector<DivdiffElement>, double, int, std::vector<int>> &coefficient,
                             double t, int q, int qubit, double chi, double omega, int target_step, int t_k) {
    auto divdiff = std::get<0>(coefficient);
    auto energy_coefficient = std::get<1>(coefficient);
    auto final_state = std::get<2>(coefficient);
    auto states = std::get<3>(coefficient);
    
    //calculate divided differences
    // divdiff_init();

    divdiffcomplex d(14,500);// max step is 14)

    // complex<ExExFloat> sum_final_elements;
    // sum_final_elements = 0;

    complex_Ex sum_final_elements = complex_Ex(0,0);

    //print divdiff size
    int num_elem = 0;
    int CurrentLength = 0;

    cal_divdiff_aux_td(d, divdiff, states, &CurrentLength, &num_elem, t, q, qubit, chi, omega, t_k, sum_final_elements);
    
    // print sum_final_elements
    // if(q == 13 && target_step == 13)
    // std::cout << "Sum final elements: " << sum_final_elements << std::endl;
    
    // return sum_final_elements * energy_coefficient;

    sum_final_elements.real = sum_final_elements.real * energy_coefficient;
    sum_final_elements.imag = sum_final_elements.imag * energy_coefficient;

    return sum_final_elements;
        
}




void cal_divdiff_aux(divdiffcomplex &d, const vector<DivdiffElement> &divdiff, const vector<int> &states, int *CurrentLength, int *num_elem, double t, int q, int qubit, double chi, double omega, complex<ExExFloat> *sum_final_elements){
    complex<double> j(0, 1);
    //first state
    if (*CurrentLength == 0){
        const auto curr_elem = divdiff[*num_elem];
        
        const double& beta_coefficient = curr_elem.coefficient;
        const auto& element_omega = curr_elem.omega[*CurrentLength];
        const auto& state_n = states[*CurrentLength];
        // std::cout << "(" << element.energy << ", " << element.omega << ", " << element.chi << ", " << element.coefficient << ") ";

        complex<double> n;
        n = (-j * t * (state_n * chi * (qubit - 1) + element_omega * omega ));
        
        // // print n
        // std::cout << "n: " << n << std::endl;
        
        d.AddElement(n);
        *CurrentLength += 1;

        cal_divdiff_aux(d, divdiff, states, CurrentLength, num_elem, t, q, qubit, chi, omega, sum_final_elements);
        return;
    }
    
    //final state
    if(*CurrentLength == (states.size() - 1)){
        for(int i=0; i < 2; i++){
            //load plus/minus
            const auto curr_elem = divdiff[*num_elem];
                
            const double& beta_coefficient = curr_elem.coefficient;
            const auto& element_omega = curr_elem.omega[*CurrentLength];
            const auto& state_n = states[*CurrentLength];
            // std::cout << "(" << element.energy << ", " << element.omega << ", " << element.chi << ", " << element.coefficient << ") ";

            complex<double> n;
            n = (-j * t * (state_n * chi * (qubit - 1) + element_omega * omega ));
            
            // // print n
            // std::cout << "n: " << n << std::endl;
            
            d.AddElement(n);
            *CurrentLength += 1;
            
            //cal z_n 
            complex<ExExFloat> z_n = d.divdiffs[*CurrentLength - 1];
            // z_n.real().print();
            // cout << endl;
            // z_n.imag().print();
            // cout << endl;

            // double real = z_n.real().get_double();
            // double imag = z_n.imag().get_double();
            //print real and imag
            // std::cout << "Divdiff = " << real << " + " << imag << "i" << std::endl;

            // //print divdiffs elements
            // d.PrintList(d.divdiffs, d.CurrentLength, "Divdiffs");
            
            complex<ExExFloat> final_element = (z_n * cal_neg_i_pow(q) * pow(t, q)) /factorial(q) * beta_coefficient;
            // complex<double> final_element_double = final_element.real().get_double() + 1j * final_element.imag().get_double() ;  
            
            // complex<ExExFloat> norm_real = z_n.real() * cal_neg_i_pow(q) * pow(t, q) / factorial(q) ;
            // complex<ExExFloat> norm_imag = z_n.imag() * (-1) * cal_neg_i_pow(q+1) * pow(t, q) / factorial(q);// (-it)^q / q! * Imag = (-it)^q / q! * (i) * Imag.real = (-i)^q+1 * t^q * -Imag.real / q!
            // complex<ExExFloat> final_element = (norm_real + norm_imag) * beta_coefficient;

            // //normalize the divided differences by (-it)^q / q! * beta_coefficient
            // complex<double> norm_real = cal_neg_i_pow(q) * std::pow(t, q) * real / factorial(q) ;
            // complex<double> norm_imag = cal_neg_i_pow(q+1) * std::pow(t, q) * (-imag) / factorial(q);// (-it)^q / q! * Imag = (-it)^q / q! * (i) * Imag.real = (-i)^q+1 * t^q * -Imag.real / q!
            // complex<double> final_element = (norm_real + norm_imag) * beta_coefficient;

            //add the divdiff element to the total sum
            *sum_final_elements = *sum_final_elements + final_element; 
            *num_elem += 1;

            d.RemoveElement();
            *CurrentLength -= 1;
        }
        return;
    }
    
    //middle states
    for(int i=0; i < 2; i++){
        const auto curr_elem = divdiff[*num_elem];
            
        const double& beta_coefficient = curr_elem.coefficient;
        const auto& element_omega = curr_elem.omega[*CurrentLength];
        const auto& state_n = states[*CurrentLength];
        // std::cout << "(" << element.energy << ", " << element.omega << ", " << element.chi << ", " << element.coefficient << ") ";

        complex<double> n;
        n = (-j * t * (state_n * chi * (qubit - 1) + element_omega * omega ));
        
        // // print n
        // std::cout << "n: " << n << std::endl;
        
        d.AddElement(n);
        *CurrentLength += 1;

        cal_divdiff_aux(d, divdiff, states, CurrentLength, num_elem, t, q, qubit, chi, omega, sum_final_elements);

        d.RemoveElement();
        *CurrentLength -= 1;
    }
    return;
}

// set divdiff calculation
complex_Ex cal_divdiff(const std::tuple<std::vector<DivdiffElement>, double, int, std::vector<int>> &coefficient,
                             double t, int q, int qubit, double chi, double omega, int target_step) {
    auto divdiff = std::get<0>(coefficient);
    auto energy_coefficient = std::get<1>(coefficient);
    auto final_state = std::get<2>(coefficient);
    auto states = std::get<3>(coefficient);
    
    //calculate divided differences
    // divdiff_init();

    divdiffcomplex d(14,500);// max step is 14)

    complex<ExExFloat> sum_final_elements;
    sum_final_elements = 0;

    //print divdiff size
    int num_elem = 0;
    int CurrentLength = 0;

    cal_divdiff_aux(d, divdiff, states, &CurrentLength, &num_elem, t, q, qubit, chi, omega, &sum_final_elements);
    
    // print sum_final_elements
    // if(q == 13 && target_step == 13)
    // std::cout << "Sum final elements: " << sum_final_elements << std::endl;
    

    complex_Ex res;
    res.real = sum_final_elements.real().get_double() * energy_coefficient;
    res.imag = sum_final_elements.imag().get_double() * energy_coefficient;

    return res;
    // return sum_final_elements * energy_coefficient;

}

/************** calculate divided differences for specific q and t****************/

void cal_divdiff_aux_orig(divdiffcomplex &d, const vector<DivdiffElement> &divdiff, const vector<int> &states, int *CurrentLength, int *num_elem, double t, int q, int qubit, double chi, double omega, complex<ExExFloat> *sum_final_elements){
    complex<double> j(0, 1);
    //first state
    if (*CurrentLength == 0){
        const auto curr_elem = divdiff[*num_elem];
        
        const double& beta_coefficient = curr_elem.coefficient;
        const auto& element_omega = curr_elem.omega[*CurrentLength];
        const auto& state_n = states[*CurrentLength];
        // std::cout << "(" << element.energy << ", " << element.omega << ", " << element.chi << ", " << element.coefficient << ") ";

        complex<double> n;
        n = (-j * t * (state_n * chi * (qubit -1) + element_omega * omega));
        
        // // print n
        // std::cout << "n: " << n << std::endl;
        
        d.AddElement(n);
        *CurrentLength += 1;

        cal_divdiff_aux_orig(d, divdiff, states, CurrentLength, num_elem, t, q, qubit, chi, omega, sum_final_elements);
        return;
    }
    
    //final state
    if(*CurrentLength == (states.size() - 1)){
        for(int i=0; i < 2; i++){
            //load plus/minus
            const auto curr_elem = divdiff[*num_elem];
                
            const double& beta_coefficient = curr_elem.coefficient;
            const auto& element_omega = curr_elem.omega[*CurrentLength];
            const auto& state_n = states[*CurrentLength];
            // std::cout << "(" << element.energy << ", " << element.omega << ", " << element.chi << ", " << element.coefficient << ") ";

            complex<double> n;
            n = (-j * t * (state_n * chi * (qubit -1) + element_omega * omega ));
            
            // // print n
            // std::cout << "n: " << n << std::endl;
            
            d.AddElement(n);
            *CurrentLength += 1;
            
            //cal z_n 
            complex<ExExFloat> z_n = d.divdiffs[*CurrentLength - 1];
            // z_n.real().print();
            // cout << endl;
            // z_n.imag().print();
            // cout << endl;

            // double real = z_n.real().get_double();
            // double imag = z_n.imag().get_double();
            //print real and imag
            // std::cout << "Divdiff = " << real << " + " << imag << "i" << std::endl;

            // //print divdiffs elements
            // d.PrintList(d.divdiffs, d.CurrentLength, "Divdiffs");
            
            complex<ExExFloat> final_element = (z_n * cal_neg_i_pow(q) * pow(t, q)) /factorial(q) * beta_coefficient;
            // complex<double> final_element_double = final_element.real().get_double() + 1j * final_element.imag().get_double() ;  
            
            // complex<ExExFloat> norm_real = z_n.real() * cal_neg_i_pow(q) * pow(t, q) / factorial(q) ;
            // complex<ExExFloat> norm_imag = z_n.imag() * (-1) * cal_neg_i_pow(q+1) * pow(t, q) / factorial(q);// (-it)^q / q! * Imag = (-it)^q / q! * (i) * Imag.real = (-i)^q+1 * t^q * -Imag.real / q!
            // complex<ExExFloat> final_element = (norm_real + norm_imag) * beta_coefficient;

            // //normalize the divided differences by (-it)^q / q! * beta_coefficient
            // complex<double> norm_real = cal_neg_i_pow(q) * std::pow(t, q) * real / factorial(q) ;
            // complex<double> norm_imag = cal_neg_i_pow(q+1) * std::pow(t, q) * (-imag) / factorial(q);// (-it)^q / q! * Imag = (-it)^q / q! * (i) * Imag.real = (-i)^q+1 * t^q * -Imag.real / q!
            // complex<double> final_element = (norm_real + norm_imag) * beta_coefficient;

            //add the divdiff element to the total sum
            *sum_final_elements = *sum_final_elements + final_element; 
            *num_elem += 1;

            d.RemoveElement();
            *CurrentLength -= 1;
        }
        return;
    }
    
    //middle states
    for(int i=0; i < 2; i++){
        const auto curr_elem = divdiff[*num_elem];
            
        const double& beta_coefficient = curr_elem.coefficient;
        const auto& element_omega = curr_elem.omega[*CurrentLength];
        const auto& state_n = states[*CurrentLength];
        // std::cout << "(" << element.energy << ", " << element.omega << ", " << element.chi << ", " << element.coefficient << ") ";

        complex<double> n;
        n = (-j * t * (state_n * chi * (qubit-1) + element_omega * omega ));
        
        // // print n
        // std::cout << "n: " << n << std::endl;
        
        d.AddElement(n);
        *CurrentLength += 1;

        cal_divdiff_aux_orig(d, divdiff, states, CurrentLength, num_elem, t, q, qubit, chi, omega, sum_final_elements);

        d.RemoveElement();
        *CurrentLength -= 1;
    }
    return;
}

// set divdiff calculation
complex<ExExFloat> cal_divdiff_orig(const std::tuple<std::vector<DivdiffElement>, double, int, std::vector<int>> &coefficient,
                             double t, int q, int qubit, double chi, double omega, int target_step) {
    auto divdiff = std::get<0>(coefficient);
    auto energy_coefficient = std::get<1>(coefficient);
    auto final_state = std::get<2>(coefficient);
    auto states = std::get<3>(coefficient);
    
    //calculate divided differences
    // divdiff_init();

    divdiffcomplex d(14,500);// max step is 14)

    complex<ExExFloat> sum_final_elements;
    sum_final_elements = 0;

    //print divdiff size
    int num_elem = 0;
    int CurrentLength = 0;

    cal_divdiff_aux_orig(d, divdiff, states, &CurrentLength, &num_elem, t, q, qubit, chi, omega, &sum_final_elements);
    
    // print sum_final_elements
    // if(q == 13 && target_step == 13)
    // std::cout << "Sum final elements: " << sum_final_elements << std::endl;
    
    return sum_final_elements * energy_coefficient;

        
}




complex_Ex cal_divdiff_const_amp(const std::tuple<double, int, std::vector<int>> &coefficient,
                             double t, int q, int qubit, double chi) {
    auto energy_coefficient = std::get<0>(coefficient);
    auto final_state = std::get<1>(coefficient);
    auto states = std::get<2>(coefficient);
    complex<double> j(0, 1);

    
    //calculate divided differences
    // divdiff_init();

    divdiffcomplex d(14,500);// max step is 14)


    //print divdiff size
    int num_elem = 0;
    int CurrentLength = 0;

    for(int i=0; i < states.size(); i++){
        const auto& state_n = states[i];
        complex<double> n;
        // n = (-j * t * (state_n * chi * (qubit - 1)));
        n = (-j * t * (state_n * chi * (qubit)));
        
        d.AddElement(n);
    }
    complex<ExExFloat> z_n = d.divdiffs[states.size() -1];
    complex<ExExFloat> final_element = (z_n * cal_neg_i_pow(q) * pow(t, q)) /factorial(q);

    complex_Ex res;
    res.real = final_element.real().get_double() * energy_coefficient;
    res.imag = final_element.imag().get_double() * energy_coefficient;

    return res;
}

complex_Ex cal_divdiff_sin_amp(const std::tuple<double, int, std::vector<int>> &coefficient,
                             double t, int q, int qubit, double chi) {
    auto energy_coefficient = std::get<0>(coefficient);
    auto final_state = std::get<1>(coefficient);
    auto states = std::get<2>(coefficient);
    complex<double> j(0, 1);

    
    //calculate divided differences
    // divdiff_init();

    divdiffcomplex d(14,500);// max step is 14)


    //print divdiff size
    int num_elem = 0;
    int CurrentLength = 0;
    int state_chi = 0;

    for(int i=0; i < states.size(); i++){
        const auto& state_n = states[i];
        complex<double> n;
        if(i == 0){
            // n = (-j * t * (state_n * chi * (qubit - 1)));
            n = (-j * t * (state_n * chi * (qubit)));
        }else{
            const auto& state_n1 = states[i-1];
            if (state_n1 < state_n){//a^dagger
                state_chi += 1;
            }else{ //a
                state_chi -= 1;
            }
            n = (-j * t * (state_n * chi * (qubit) + state_chi * chi));
        }
        d.AddElement(n);
    }
    complex<ExExFloat> z_n = d.divdiffs[states.size() -1];
    complex<ExExFloat> final_element = (z_n * cal_neg_i_pow(q) * pow(t, q)) /factorial(q);

    complex_Ex res;
    res.real = final_element.real().get_double() * energy_coefficient;
    res.imag = final_element.imag().get_double() * energy_coefficient;

    return res;
}

complex_Ex cal_divdiff_jc(const std::tuple<double, step, std::vector<step>, step> &coefficient,
                             double t, int q, int qubit, double chi, double alpha, double delta) {
    // auto energy_coefficient = std::get<0>(coefficient);
    // auto final_state = std::get<1>(coefficient);
    auto states = std::get<2>(coefficient);
    complex<double> j(0, 1);

    //calculate divided differences
    // divdiff_init();
    divdiffcomplex d(14,500);// max step is 14)

    //print divdiff size
    int num_elem = 0;
    int CurrentLength = 0;
    int curr_chi = 0;
    int curr_delta = 0;

    for(int i=0; i < states.size(); i++){
        complex<double> n;
        
        double diagonal_element = (alpha) * states[i].qubit * (states[i].qubit - 1);
        
        if(i == 0){
            // n = (-j * t * (state_n * chi * (qubit - 1)));
            n = (-j * t * diagonal_element);
        }else{
             // Update energy_coefficient_real
             if (states[i].cavity - states[i-1].cavity == 1 && states[i].qubit == states[i-1].qubit) { // a^dagger
                curr_chi += 1;
            }
            else if (states[i].cavity - states[i-1].cavity == -1 && states[i].qubit == states[i-1].qubit) { // a
                curr_chi -= 1; 
            }
            else if (states[i].cavity - states[i-1].cavity == 1 && states[i].qubit - states[i-1].qubit == -1) { // a^dagger q
                curr_delta += 1;
            } 
            else if (states[i].cavity - states[i-1].cavity == -1 && states[i].qubit - states[i-1].qubit == 1) { // a q^dagger
                curr_delta -= 1;
            }
            
            n = (-j * t * (curr_chi * chi + curr_delta * delta + diagonal_element));
        }
        d.AddElement(n);
    }
    complex<ExExFloat> z_n = d.divdiffs[states.size() - 1];
    complex<ExExFloat> final_element = (z_n * cal_neg_i_pow(q) * pow(t, q)) /factorial(q);

    complex_Ex res;//TODO: check if need complex_Ex or complex<ExFloat>
    res.real = final_element.real().get_double();
    res.imag = final_element.imag().get_double();

    return res;
}






// clear up the divdiff memory


complex_Ex operator+(complex_Ex a, complex_Ex b) {
    divdiff_init();
    complex_Ex res;

    ExExFloat temp1;
    temp1 = a.real + b.real;
    ExExFloat temp2 = b.imag + b.imag;

    res.real = temp1;
    res.imag = temp2;
    divdiff_clear_up();

    return res;
}

complex_Ex operator*(complex_Ex a, complex_Ex b) {
   divdiff_init();
    complex_Ex res;
    ExExFloat temp1 = a.real * b.real;
    ExExFloat temp2 = a.imag * b.imag;

    res.real = temp1 - temp2;
    

    ExExFloat temp3 = a.real * b.imag;
    ExExFloat temp4 = a.imag * b.real;

    res.imag = temp3 + temp4;
    divdiff_clear_up();

    return res;
}

complex_Ex operator*( complex_Ex a, double x) {
    divdiff_init();
    complex_Ex res;

    res.real = a.real * x;
    res.imag = a.imag * x;
    divdiff_clear_up();

    return res;
}

double abs2(complex_Ex z) {
    ExExFloat res1 = z.real * z.real;
    ExExFloat res2 = z.imag * z.imag;
    ExExFloat res = res1 + res2;
    return res.get_double();
}

complex_Ex conj(complex_Ex z) {
    complex_Ex res;

    res.real = z.real;
    res.imag = z.imag * (-1);

    return res;
}




#endif // PERMUTATION_H