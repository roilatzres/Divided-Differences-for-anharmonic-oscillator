#ifndef PERMUTATION_H
#define PERMUTATION_H

#include <iostream>
#include <vector>
#include <functional>
#include <cmath>
#include <complex>
#include "dlib/dlib/optimization.h"
#include "dlib/dlib/matrix.h"


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

complex_Ex complex_mult(complex_Ex a, complex_Ex b){
    divdiff_init();
    complex_Ex res;
    //print res.real
    ExExFloat temp1 = a.real * b.real;
    ExExFloat temp2 = a.imag * b.imag;
    // //print temp
    // std::cout << "temp1: " << temp1.get_double() << std::endl;
    // std::cout << "temp2: " << temp2.get_double() << std::endl;


    res.real = temp1 - temp2;
    // //print res.real
    // std::cout << "res.real: " << res.real.get_double() << std::endl;
    

    ExExFloat temp3 = a.real * b.imag;
    ExExFloat temp4 = a.imag * b.real;
    // //print temp
    // std::cout << "temp3: " << temp3.get_double() << std::endl;
    // std::cout << "temp4: " << temp4.get_double() << std::endl;

    res.imag = temp3 + temp4;
    // //print res.imag
    // std::cout << "res.imag: " << res.imag.get_double() << std::endl;

    divdiff_clear_up();

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
struct Step {
    int first;
    int second;
};

// Function to compute all possible permutations of steps to reach the target step on a ladder
vector<vector<Step>> jc_ladder_permutations(int target_first, int target_second, int moves) {
    vector<vector<Step>> all_combinations;

    // Backtrack function to generate permutations recursively
    std::function<void(Step, int, vector<Step>&)> jc_backtrack = [&](Step curr_step, int remaining_moves, vector<Step>& combination) {
        // Base case: If remaining moves is 0 and current step is target step, add current combination
        if (remaining_moves == 0) {
            if (curr_step.first == target_first && curr_step.second == target_second) {
                all_combinations.push_back(combination);
            }
            return;
        }

        // Try moving first + 1, second
        combination.push_back({curr_step.first + 1, curr_step.second});
        jc_backtrack({curr_step.first + 1, curr_step.second}, remaining_moves - 1, combination);
        combination.pop_back();

        // Try moving first - 1, second
        combination.push_back({curr_step.first - 1, curr_step.second});
        jc_backtrack({curr_step.first - 1, curr_step.second}, remaining_moves - 1, combination);
        combination.pop_back();

        // Try moving first + 1, second - 1
        combination.push_back({curr_step.first + 1, curr_step.second - 1});
        jc_backtrack({curr_step.first + 1, curr_step.second - 1}, remaining_moves - 1, combination);
        combination.pop_back();

        // Try moving first - 1, second + 1
        combination.push_back({curr_step.first - 1, curr_step.second + 1});
        jc_backtrack({curr_step.first - 1, curr_step.second + 1}, remaining_moves - 1, combination);
        combination.pop_back();
    };

    // Initialize combination list
    vector<Step> combination;
    Step start = {0, 0}; // Start from {0, 0}
    // Start backtracking from step {0, 0} with given number of moves
    jc_backtrack(start, moves, combination);

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

    // TODO: check why the E-3 factor is needed
    // FIXME: the factor is amplitude!!

    int final_state = current_state;

    // Initialize beta list with the final state
    DivdiffElement first_elem;
    first_elem.omega = {0};
    first_elem.coefficient = 1;
    vector<DivdiffElement> beta;
    beta.push_back(first_elem); 

    // iterate from end to start over the states and update the divdiff elements starting from the penultimate state
    for (int i = states.size() - 2; i >= 0; --i) {
        vector<DivdiffElement> temp_divdiff;

        // Update divdiff elements
        for (int j = 0; j < beta.size(); ++j) {
            // Extract the first element from beta[j]
            DivdiffElement beta_element = beta[j];
            
            DivdiffElement plus, minus;
            plus.coefficient = -beta_element.coefficient;
            plus.omega = beta_element.omega; // TODO: check if vector placement is correct
            
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

    // TODO: check why the E-3 factor is needed
    // FIXME: the factor is amplitude!!

    int final_state = current_state;
    return std::make_tuple(energy_coefficient_real, final_state, states);

}


/************** calculate divided differences for specific q and t****************/

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
complex<ExExFloat> cal_divdiff(const std::tuple<std::vector<DivdiffElement>, double, int, std::vector<int>> &coefficient,
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

    //TODO: remove d from memory
    complex_Ex res;
    res.real = final_element.real().get_double() * energy_coefficient;
    res.imag = final_element.imag().get_double() * energy_coefficient;

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
    divdiff_init();
    ExExFloat res1 = z.real * z.real;
    ExExFloat res2 = z.imag * z.imag;
    ExExFloat res = res1 + res2;
    divdiff_clear_up();
    return res.get_double();
}

complex_Ex conj(complex_Ex z) {
    divdiff_init();
    complex_Ex res;

    res.real = z.real;
    res.imag = z.imag * (-1);
    divdiff_clear_up();

    return res;
}

// ==== Typedefs ====
using mat_c = vector<vector<complex_Ex>>;
using vec_c = vector<complex_Ex>;

// ==== Polynomial evaluation ====
mat_c poly_eval(double x, const vector<mat_c>& C) {
    int N = C[0].size();
    mat_c result(N, vec_c(N, complex_Ex(0, 0)));

    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            for (size_t k = 0; k < C.size(); ++k)
                result[i][j] = result[i][j] + C[k][i][j] * pow(x, k);

    return result;
}

// ==== Matrix multiply ====
mat_c matmul(const mat_c& A, const mat_c& B) {
    int N = A.size();
    mat_c result(N, vec_c(N, complex_Ex(0, 0)));

    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            for (int k = 0; k < N; ++k)
                result[i][j] = result[i][j] + A[i][k] * B[k][j];

    return result;
}

// ==== Matrix-vector multiply ====
vec_c matvec(const mat_c& M, const vec_c& v) {
    int N = M.size();
    vec_c result(N, complex_Ex(0, 0));

    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            result[i] = result[i] + M[i][j] * v[j];

    return result;
}

// ==== Inner product ====
complex_Ex inner(const vec_c& a, const vec_c& b) {
    complex_Ex sum(0, 0);
    for (size_t i = 0; i < a.size(); ++i)
        sum = sum + conj(a[i]) * b[i];
    return sum;
}

// ==== Objective function ====
class Objective {
public:
    vector<mat_c> C;
    vec_c start, target;

    Objective(const vector<mat_c>& C_in, const vec_c& s, const vec_c& t)
        : C(C_in), start(s), target(t) {}

    double operator()(const dlib::matrix<double, 0, 1>& x) const {
        mat_c M = identity_matrix(C[0].size());

        for (size_t i = 0; i < x.size(); ++i)
            M = matmul(M, poly_eval(x(i), C));

        vec_c psi = matvec(M, start);

        // Normalize psi
        double norm_psi = 0.0;
        for (const auto& z : psi) norm_psi += abs2(z);
        norm_psi = sqrt(norm_psi);

        if (norm_psi > 1e-12) {
            for (auto& z : psi)
                z = z * (1.0 / norm_psi);
        }

        complex_Ex amp = inner(target, psi);
        return 1.0 - abs2(amp);
    }

    static mat_c identity_matrix(int N) {
        mat_c I(N, vec_c(N, complex_Ex(0, 0)));
        for (int i = 0; i < N; ++i)
            I[i][i] = complex_Ex(1, 0);
        return I;
    }
};





#endif // PERMUTATION_H