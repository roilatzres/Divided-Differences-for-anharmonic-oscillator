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
#include <chrono>


// #include "main_orig.h"

using std::vector;
using std::complex;

#define  MAX_STEP_CAV 300
#define  MAX_STEP_QUBIT 2
// #define  MAX_TARGET 15
#define MAX_Q 30
#define M_PI 3.14159265358979323846


//function to print complex_ex
void print_complex_Ex(complex_Ex a){
    a.real.print();
    std::cout << "+ ";
    a.imag.print();
    std::cout << "i" << std::endl; 
    std::cout << std::endl;
}


// int main(int argc, char* argv[]) {
vector<vector<vector<complex_Ex>>> solve_for_t(double sample_t, 
    int q_max,
    double amplitude, 
    step start_state, 
    step max_target, 
    double alpha, 
    double delta, 
    double coupling, 
    double chi, 
    double timestep,
    step max_total = {MAX_STEP_CAV, MAX_STEP_QUBIT}
    // int num_pulses
    ) {

    complex<double> i_num(0, 1);
    // //original parameters
    // double multiplier = 4;
    // float sigma = 48;
    // float timestep = 1;
    // double chi = -279e-6 * 2 * M_PI; 
    // double detuning = chi *1e-2;


    
    // //print args
    // std::cout << "Amplitude: " << amplitude << std::endl;
    // std::cout << "q_max: " << q_max << std::endl;

    // //print vars
    // std::cout << "Multiplier: " << multiplier << std::endl;
    // std::cout << "Sigma: " << sigma << std::endl;
    // std::cout << "Timestep: " << timestep << std::endl;
    // std::cout << "Chi: " << chi << std::endl;
    // // std::cout << "Detuning: " << detuning << std::endl;


    // double tfinal = multiplier * sigma;
    // double omega = (2 * M_PI) / tfinal;

    // std::cout << "tfinal: " << tfinal << std::endl;
    // std::cout << "Omega: " << omega << std::endl;
    vector<double> times;
    for (double t = 1; t <= sample_t; t += timestep) {
        times.push_back(t);
        //print times
        // std::cout << "Time: " << t << std::endl;
    }

    // //calculate the time evolution
    double t;

    // vector of all transition amplitudes for all times: [qubit_space][time][cavity_space][amp]
    // vector<vector<vector<complex_Ex>>> all_transition_amplitudes(max_target.qubit, vector<vector<complex_Ex>>(times.size(), vector<complex_Ex>(max_target.cavity)));
    vector<vector<vector<complex_Ex>>> all_transition_amplitudes(max_target.qubit);
    vector<vector<complex<double>>> all_transition_amplitudes_g;
    vector<vector<complex<double>>> all_transition_amplitudes_e;
    
    //define start state
    // int start_state = 0;
    int total_steps = 100;

    // Create a 3D vector to store the results
    vector< //target_step
    vector< //q
    vector< //permutations
    std::tuple<vector<DivdiffElement>, double, int, vector<int>>> // coefficients
    >> all_coefficients;
    
    // load_all_coefficients(all_coefficients, start_state, total_steps, MAX_TARGET, q_max);
    
    // //loop over all times
    for (int i = 0; i < sample_t; i++) { 
        t = times[i];
        //     // t = 96;
        //print time
        // std::cout << std::endl << "Main-Time: " << t << std::endl;
        
        for (int qubit = 0; qubit < max_target.qubit; qubit++){
            // print qubit
            // std::cout << "Qubit: " << qubit << std::endl;

            //vector of all transition amplitude for a given time
            auto transition_amplitudes = vector<complex_Ex> {};
            
            //loop over all target steps
            for (int target = 0; target < max_target.cavity; target++) {
                // print target
                // std::cout << "Target: " << target << std::endl;

                complex_Ex total_amp;
                total_amp.real = 0;
                total_amp.imag = 0;

                // //print total_amp
                // std::cout <<  "init Total amp: "; 
                // print_complex_Ex(total_amp);

                // //print target
                // std::cout << "Target: " << target << std::endl;

                for(int q = 0; q <= q_max; q++){
                    // print q
                    // std::cout << "Q: " << q << std::endl;

                    complex_Ex q_amp; 
                    q_amp.real = 0;
                    q_amp.imag = 0;

                    // double energy_coefficient_amp = std::pow(amplitude, q);

                    if (target == start_state.cavity && qubit == start_state.qubit && q == 0){
                        // cout << "Skipping q = 0 for target = " << target << " and qubit = " << qubit << std::endl;
                        complex<double> exponent = std::exp(-i_num * (double)(t * qubit * (qubit - 1) * alpha));
                        total_amp.real = exponent.real();
                        total_amp.imag = exponent.imag();
                        continue;
                    }

                    // std::string q_amp_filename = "./include/q_amp/jc/"+ std::to_string(num_pulses) 
                    //     +"/start_state_" + std::to_string(start_state.cavity) 
                    //     + "_qubit_" + std::to_string(start_state.qubit) 
                    //     + "_t_" + std::to_string((int)t) 
                    //     + "_q" + std::to_string(q)
                    //     + "_qubit_" + std::to_string(qubit) 
                    //     + "_target_" + std::to_string(target)  + ".bin";
                    // std::ifstream infile_q_amp(q_amp_filename);

                    // if (infile_q_amp) {
                    //     q_amp = load_binary_complexEx(q_amp_filename);
                    //     q_amp.real *= energy_coefficient_amp; 
                    //     q_amp.imag *= energy_coefficient_amp; 

                    //     divdiff_init();
                    //     total_amp.real += q_amp.real.get_double();
                    //     total_amp.imag += q_amp.imag.get_double();
                    //     divdiff_clear_up();
                    //     continue;
                    // }

                    
                    std::string perm_dir = "./include/permutations/jc/" + std::to_string(max_total.cavity) + "c_" + std::to_string(max_total.qubit) + "q/";
                    // cout<<"perm_dir: " << perm_dir << endl;
                    if (!std::filesystem::exists(perm_dir)) {
                        std::filesystem::create_directories(perm_dir);
                        cout << "Created directory: " << perm_dir << endl;
                    }

                    auto permutations = vector<vector<step>> {};
                    std::string perm_filename = perm_dir + "start_state" +
                    std::to_string(start_state.cavity) + "_qubit_" + std::to_string(start_state.qubit) + "_target_" + std::to_string(target) + "_q" + std::to_string(q) + ".txt";
                    
                    
                    std::ifstream infile_perm(perm_filename);
                    if (infile_perm) {
                        permutations = load_binary_perm_jc(perm_filename);
                        // print permutations 
                        
                        // std::cout << "Permutations loaded from file: " << perm_filename << std::endl;
                        // std::cout << "Number of permutations: " << permutations.size() << std::endl;
                        // for (const auto& perm : permutations) {
                        //     std::cout << "Permutation: ";
                        //     for (const auto& step : perm) {
                        //         std::cout << "(" << step.cavity << ", " << step.qubit << ") ";
                        //     }
                        // }
                    } else {
                        // step max_step = {300, 8};
                        permutations = jc_ladder_permutations(start_state, max_total, {target, qubit}, q);
                        // for (const auto& perm : permutations) {
                        //     std::cout << "Permutation: ";
                        //     for (const auto& step : perm) {
                        //         std::cout << "(" << step.cavity << ", " << step.qubit << ") ";
                        //     }
                        // }
                        if (permutations.empty()) {
                            // std::cerr << "No permutations found for target " << target << " and qubit " << qubit << " with q = " << q << std::endl;
                            continue;
                        } else {
                            save_binary_perm_jc(permutations, perm_filename);
                        }
                    }
                       
                    vector<std::tuple<double, step, vector<step>>> all_coefficients_const_pulse; //FIXME: decide whether to save coefficients - depends on amp and coupling
                    // std::string all_coef_filename = "./include/parms/dispersive_const/start_state_" +
                    // std::to_string(start_state.cavity) + "_qubit_" + std::to_string(start_state.qubit) + "_target_" + std::to_string(target) + "_q" + std::to_string(q) + ".txt";
                    
                    // std::ifstream infile_coef(all_coef_filename);
                    // if (infile_coef) {
                    //     all_coefficients_const_pulse = load_binary_coef_jc(all_coef_filename);
                    // } else {
                    //     for (const auto& permutation : permutations) {
                    //         auto coefficients = cal_coefficient_const_pulse_jc(permutation, start_state, amplitude, coupling); //TODO: amplitude and coupling can be calculated in advance but it depends on permutation
                    //         all_coefficients_const_pulse.push_back(coefficients);
                    //     }
                    //     save_binary_coef_jc(all_coefficients_const_pulse, all_coef_filename);
                    // }
                    // cout << "pre cal coefficients" << std::endl;
                    //cal time for loop
                    

                    //print permutations size
                    // std::cout << "Number of permutations: " << permutations.size() << std::endl;
                    for (const auto& permutation : permutations) {
                        auto coefficients = cal_coefficient_const_pulse_jc(permutation, start_state, amplitude, coupling); //TODO: amplitude and coupling can be calculated in advance but it depends on permutation
                        // // print coefficients
                        // std::cout << "Coefficient: " << std::get<0>(coefficients) << std::endl;
                        // std::cout << "Final state: " << std::get<1>(coefficients).cavity << ", " << std::get<1>(coefficients).qubit << std::endl;
                        // std::cout << "States: " << std::get<2>(coefficients).size() << std::endl;
                        // for (const auto& state : std::get<2>(coefficients)) {
                        //     std::cout << "State: " << state.cavity << ", " << state.qubit << std::endl;
                        // }
                        all_coefficients_const_pulse.push_back(coefficients);
                    }
                    
                    
                    
                    
                    // cout << "pre divdiff_init" << std::endl;

                    
                    divdiff_init();
                    //print all_coeeficients_const_pulse size
                    // std::cout << "Number of coefficients: " << all_coefficients_const_pulse.size() << std::endl;
                    for (const auto& coefficient : all_coefficients_const_pulse) {
                        complex_Ex res = cal_divdiff_sin_amp_jc(coefficient, t, q, qubit, chi, alpha, delta);
                        
                        // print_complex_Ex(res);

                        q_amp.real += res.real;
                        q_amp.imag += res.imag;
                    }
                    divdiff_clear_up();


                    // save_binary_complexEx(q_amp, q_amp_filename);
                    // q_amp.real *= energy_coefficient_amp; 
                    // q_amp.imag *= energy_coefficient_amp; 

                    divdiff_init();
                    total_amp.real += q_amp.real.get_double();
                    total_amp.imag += q_amp.imag.get_double();
                    divdiff_clear_up();
                }

                // //print the Total amplitude
                // // std::cout << "Total amplitude: " << total_amp << std::endl << std::endl;
                // std::cout <<  "Total amp: "; 
                // print_complex_Ex(total_amp);
                transition_amplitudes.push_back(total_amp);
            }

            // // //print transition amplitudes
            // std::cout << "Transition amplitudes: " << std::endl;
            // int target_nums = 0;
            // for (auto& amp : transition_amplitudes) {
            //     std::cout << target_nums << ": "; 
            //     print_complex_Ex(amp);
            //     target_nums++;
            // }

            //end of qubit loop
            all_transition_amplitudes[qubit].push_back(transition_amplitudes);
        }
        
        double norm = 0;
        ExExFloat Ex_norm;
        Ex_norm = 0;
        divdiff_init();
        for (int qubit = 0; qubit < max_target.qubit; qubit++) {
            auto& all_transition_amplitudes_qubit = all_transition_amplitudes[qubit][i];
            for (auto& amp : all_transition_amplitudes_qubit) {
                Ex_norm += amp.real * amp.real + amp.imag * amp.imag;
            }
        }
        norm = sqrt(Ex_norm.get_double());
        for(int qubit = 0; qubit < max_target.qubit; qubit++) {
            auto& transition_amplitudes = all_transition_amplitudes[qubit][i];
            for (auto& amp : transition_amplitudes) {
                amp.real = amp.real / norm;
                amp.imag = amp.imag / norm;
                // print_complex_Ex(amp);
            }
        }
        divdiff_clear_up();
        

        
        
        //end of time loop
    }
    // //print all transition amplitudes with location in the vector
    // std::cout << "All transition amplitudes: " << std::endl;
    // for (int qubit = 0; qubit < max_target.qubit; qubit++) {
    //     std::cout << "Qubit: " << qubit << std::endl;
    //     for (int i = 0; i < all_transition_amplitudes[qubit].size(); i++) {
    //         std::cout << "Time: " << times[i] << std::endl;
    //         for (int j = 0; j < all_transition_amplitudes[qubit][i].size(); j++) {
    //             std::cout << "Target: " << j << " Amplitude: ";
    //             print_complex_Ex(all_transition_amplitudes[qubit][i][j]);
    //         }
    //     }
    // }
    
    std::cout << "Exiting solve_for_t" << std::endl;

    return all_transition_amplitudes;
}

vector<double> cal_amps(double amplitude, double chi, double t, double detuning, double sample_t){
    complex<double> i_num(0, 1);
    
    vector<double> amps;

    for(int i = 1; i <= t+1; i+=sample_t){
        double sin_omega_t = amplitude * std::sin(detuning * i );
        // //print sin_omega_t
        // std::cout << "sin_omega_t: " << sin_omega_t << std::endl;
        // complex<double> expo1 = std::exp(-i_num * chi * (double)i);
        // complex<double> expo2 = std::exp(i_num * chi * (double)i);
        // //print expo
        // std::cout << "expo: " << expo << std::endl;
        // complex<double> amp = sin_omega_t * (expo1+expo2);
        double amp = sin_omega_t;
        amps.push_back(amp);
    }
    return amps;
}

void accumulate_states(
    vector<vector<vector<vector<complex_Ex>>>>& all_state_ta,
    const vector<double>& all_amps,
    int num_pulses,
    int q_max,
    step max_target,
    double sample_t,
    double alpha,
    double delta,
    double coupling,
    double chi,
    double timestep,
    step max_total = {MAX_STEP_CAV, MAX_STEP_QUBIT}
) {
    for(int curr_pulse = 1; curr_pulse < num_pulses ; curr_pulse++){
        auto start_pulse = std::chrono::high_resolution_clock::now();


        std::cout << "Current pulse: " << curr_pulse << std::endl;
        
        auto curr_state_ta = vector<vector<vector<vector<vector<complex_Ex>>>>>{};
        
        //create curr_state_ta with dimensions [start_cavity][start_qubit][time][qubit_space][cavity_space]
        curr_state_ta.resize(max_target.cavity);
        for(int i = 0; i < max_target.cavity; i++){
            curr_state_ta[i].resize(max_target.qubit);
            for(int j = 0; j < max_target.qubit; j++){
                curr_state_ta[i][j].resize(sample_t);
                for(int k = 0; k < sample_t; k++){
                    curr_state_ta[i][j][k].resize(max_target.qubit);
                    for(int l = 0; l < max_target.qubit; l++){
                        curr_state_ta[i][j][k][l].resize(max_target.cavity);
                        for(int m = 0; m < max_target.cavity; m++){
                            curr_state_ta[i][j][k][l][m] = complex_Ex(0,0);
                        }
                    }
                }
            }
        }

        for (int start_cavity = 0; start_cavity < max_target.cavity; start_cavity++){
            // cout << "Start state: " << start_cavity << endl;
            // cout << endl;

            //save to final state
            for(int start_qubit = 0; start_qubit < max_target.qubit; start_qubit++){
                // cout << "Start qubit: " << start_qubit << endl;
                // cout << endl;

                step curr_start_state = {start_cavity, start_qubit};
                auto curr_state = vector<vector<vector<complex_Ex>>>{};
                
                auto start_solve_t = std::chrono::high_resolution_clock::now();

                curr_state = solve_for_t(sample_t, q_max, all_amps[curr_pulse], curr_start_state, max_target, alpha, delta, coupling, chi, timestep, max_total);
                auto end_solve_t = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> elapsed = end_solve_t - start_solve_t;
                std::cout << "Elapsed time for solve_t: " << elapsed.count() << " seconds" << std::endl;
                // // Ensure curr_state_ta is large enough//TODO: create curr_state_ta with correct size from the start
                // if (curr_state_ta.size() <= start_cavity)
                //     curr_state_ta.resize(start_cavity + 1);
                // if (curr_state_ta[start_cavity].size() <= start_qubit)
                //     curr_state_ta[start_cavity].resize(start_qubit + 1);

                curr_state_ta[start_cavity][start_qubit] = curr_state;
                // cout << "starting calculation" << endl;

                for(int i=0; i < int(sample_t); i++){
                    for (int j = 0; j < max_target.cavity; j++){ 
                        for (int k = 0; k < max_target.qubit; k++){
                            // cout << "i: " << i << ", j: " << j << ", k: " << k << endl;
                            
                            curr_state_ta[start_cavity][start_qubit][k][i][j] = complex_mult(
                                curr_state_ta[start_cavity][start_qubit][k][i][j],
                                all_state_ta[curr_pulse-1][start_qubit][sample_t - 1][start_cavity]
                            );

                            // cout << "starting divdiff" << endl;
                            divdiff_init();
                            all_state_ta[curr_pulse][k][i][j].real += curr_state_ta[start_cavity][start_qubit][k][i][j].real;
                            all_state_ta[curr_pulse][k][i][j].imag += curr_state_ta[start_cavity][start_qubit][k][i][j].imag;
                            divdiff_clear_up();
                        }
                    }
                }
                // cout << "finished calculation" << endl;
            }
        }

        auto end_pulse = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end_pulse - start_pulse;
        std::cout << "Elapsed time for pulse "<< curr_pulse <<": " << elapsed.count() << " seconds" << std::endl;

    }
}
void save_all_transition_amplitudes_single_file(
    std::vector<std::vector<std::vector<std::vector<complex_Ex>>>>& all_state_ta,
    int num_pulses,
    int max_target_qubit,
    std::string& filename
) {
    nlohmann::json j;
    for (int pulse = 0; pulse < num_pulses; ++pulse) {
        nlohmann::json j_pulse;
        for (int qubit = 0; qubit < max_target_qubit; ++qubit) {
            nlohmann::json j_qubit;
            for (auto& vec : all_state_ta[pulse][qubit]) {
                nlohmann::json sub_j;
                for (auto& c : vec) {
                    sub_j.push_back({{"real", c.real.get_double()}, {"imag", c.imag.get_double()}});
                }
                j_qubit.push_back(sub_j);
            }
            j_pulse.push_back(j_qubit);
        }
        j.push_back(j_pulse);
    }
    std::ofstream o(filename);
    o << j.dump(4) << std::endl; // Pretty print with indent
}


void save_all_state_ta_to_json(
    vector<vector<vector<vector<complex_Ex>>>>& all_state_ta,
    int num_pulses,
    double amplitude,
    int q_max,
    int max_target_cavity,
    int max_target_qubit
) {
    for(int pulse = 0; pulse < num_pulses; pulse++){
        nlohmann::json j;
        for(int qubit = 0; qubit < max_target_qubit; qubit++){
            j.clear();
            for (auto& vec : all_state_ta[pulse][qubit]) {
                nlohmann::json sub_j;
                for (auto& c : vec) {
                    sub_j.push_back({{"real", c.real.get_double()}, {"imag", c.imag.get_double()}});
                }
                j.push_back(sub_j);
            }

            std::string filename;
            if (qubit == 0){
                filename = "jc_runs/amplitudes_g_chg_pulse" + std::to_string(pulse) + "amp" + std::to_string(amplitude) + "_q" + std::to_string(q_max) + "_nc" + std::to_string(max_target_cavity) + "_nq" + std::to_string(max_target_qubit) + "_p" + std::to_string(num_pulses) + ".json";
            }else if (qubit == 1){
                filename = "jc_runs/amplitudes_e_chg_pulse" + std::to_string(pulse) + "amp" + std::to_string(amplitude) + "_q" + std::to_string(q_max) + "_nc" + std::to_string(max_target_cavity) + "_nq" + std::to_string(max_target_qubit) + "_p" + std::to_string(num_pulses) + ".json";
            }
            std::ofstream o(filename);
            o << j << std::endl;
        }
    }
}


int main(int argc, char* argv[]) {
    std::cout << "Program started\n";
    // Load parameters from file
    nlohmann::json params;
    std::ifstream params_file("params.json");
    if (!params_file.is_open()) {
        std::cerr << "Could not open params.json!" << std::endl;
        return 1;
    }
    params_file >> params;

    // Assign parameters
    double final_t = params["final_t"];
    double multiplier = params["multiplier"];
    float sigma = params["sigma"];
    float timestep = params["timestep"];
    double chi = params["chi"]; // multiply by 2*M_PI below
    double alpha = params["alpha"]; // multiply by 2*M_PI below
    double delta = params["delta"]; // multiply by 2*M_PI below
    double amplitude = params["amplitude"];
    int q_max = params["q_max"];
    int max_target_cavity = params["max_target_cavity"];
    int max_target_qubit = params["max_target_qubit"];
    int start_state_cavity = params["start_state_cavity"];
    int start_state_qubit = params["start_state_qubit"];
    int num_pulses = params["num_pulses"];
    std::string out_dir = params["out_dir"];
    int max_total_cavity = params["max_total_cavity"];
    int max_total_qubit = params["max_total_qubit"];

    // Apply 2*M_PI where needed
    chi *= 2 * M_PI;
    alpha *= 2 * M_PI;
    delta *= 2 * M_PI;

    std::cout << "Parameters:" << std::endl;
    std::cout << "final_t: " << final_t << std::endl;
    std::cout << "multiplier: " << multiplier << std::endl;
    std::cout << "sigma: " << sigma << std::endl;
    std::cout << "timestep: " << timestep << std::endl;
    std::cout << "chi: " << chi << std::endl;
    std::cout << "alpha: " << alpha << std::endl;
    std::cout << "delta: " << delta << std::endl;
    std::cout << "amplitude: " << amplitude << std::endl;
    std::cout << "q_max: " << q_max << std::endl;
    std::cout << "max_target_cavity: " << max_target_cavity << std::endl;
    std::cout << "max_target_qubit: " << max_target_qubit << std::endl;
    std::cout << "start_state_cavity: " << start_state_cavity << std::endl;
    std::cout << "start_state_qubit: " << start_state_qubit << std::endl;
    std::cout << "num_pulses: " << num_pulses << std::endl;
    std::cout << "out_dir: " << out_dir << std::endl;
    std::cout << "max_total_cavity: " << max_total_cavity << std::endl;
    std::cout << "max_total_qubit: " << max_total_qubit << std::endl;
    
    
    
    // // initial values
    // double amplitude = 0.002;
    // int q_max = 10; 
    // int max_target_cavity = 8;
    // int max_target_qubit = 2;
    // int num_pulses = 1;
    
    // //get the parameters
    // if(argc > 1){
    //     amplitude = std::stod(argv[1]);
    //     q_max = std::stoi(argv[2]);
    //     max_target_cavity = std::stoi(argv[3]);
    //     max_target_qubit = std::stoi(argv[4]);
    //     num_pulses = std::stoi(argv[5]);
    // }
    // //original parameters TODO: recheck parameters with eliya
    // double multiplier = 4;
    // float sigma = 48;
    // float timestep = 1;
    // double chi = -279e-6 * 2 * M_PI; 
    // double alpha = -200e-3 * 2 * M_PI; 
    // double delta = 2 * 2 * M_PI; 

    double sample_t = final_t / num_pulses;
    step max_target = {max_target_cavity, max_target_qubit};
    step start_state = {start_state_cavity, start_state_qubit}; //start state is always 0,0
    step max_total = {max_total_cavity, max_total_qubit};

    double coupling = sqrt(abs(delta * (delta + alpha) * chi / alpha));

    vector<double> all_amps = cal_amps(amplitude, chi, final_t,  2 * M_PI / final_t, sample_t);
    
    //print amplitudes
    std::cout << "Amplitudes: " << std::endl;
    for (int i = 0; i < all_amps.size(); i++) {
        std::cout << "Amplitude " << i << ": " << all_amps[i] << std::endl;
    }
    
    //init all_state_ta 
    vector<vector<vector<vector<complex_Ex>>>> all_state_ta;
    for(int state=0; state < num_pulses; state++){
        vector<vector<vector<complex_Ex>>> final_state_ta;
        for (int qubit = 0; qubit < max_target_qubit; qubit++){
            vector<vector<complex_Ex>> qubit_state;
            for (int i = 0; i < sample_t; i++){// times size is 1
                vector<complex_Ex> state;
                for (int j = 0; j < max_target_cavity; j++){
                    complex_Ex c;
                    c.real = 0;
                    c.imag = 0;
                    state.push_back(c);
                }
                qubit_state.push_back(state);
            }
            final_state_ta.push_back(qubit_state);
        }
        all_state_ta.push_back(final_state_ta);
    }

    // first state is always 0,0
    all_state_ta[0] = solve_for_t(sample_t, q_max, all_amps[0], start_state, max_target, alpha, delta, coupling, chi, timestep, max_total);

    accumulate_states(all_state_ta, all_amps, num_pulses, q_max, max_target, sample_t, alpha, delta, coupling, chi, timestep, max_total);

    // //print final state
    // std::cout << "Final state!!!" << std::endl;
    // for(int qubit = 0; qubit < 1; qubit++){
    //     std::cout << "qubit: " << qubit << std::endl;
    //     for (int j = 0; j < MAX_TARGET; j++){
    //         std::cout << "Target: " << j << std::endl;
    //         print_complex_Ex(final_state_ta[qubit][0][j]);
    //     }
    // }
    
    //save final state to a file
    std::cout << "saving to file" << std::endl;
    std::string filename = out_dir + 
            "amplitudes_" + std::to_string(amplitude) +
            "_q" + std::to_string(q_max) +
            "_nc" + std::to_string(max_target_cavity) +
            "_nq" + std::to_string(max_target_qubit) + 
            "_p" + std::to_string(num_pulses) + ".json";
    
    save_all_transition_amplitudes_single_file(all_state_ta, num_pulses, max_target_qubit, filename);


    // save_all_state_ta_to_json(all_state_ta, num_pulses, amplitude, q_max, max_target_cavity, max_target_qubit);
    
    
    std::cout << "Done saving to file" << std::endl;

    return 0;
}
