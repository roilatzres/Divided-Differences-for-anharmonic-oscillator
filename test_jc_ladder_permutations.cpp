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

void print_permutations(const std::vector<std::vector<step>>& perms) {
    for (size_t i = 0; i < perms.size(); ++i) {
        std::cout << "Permutation " << i << ": ";
        for (const auto& s : perms[i]) {
            std::cout << "(" << s.cavity << "," << s.qubit << ") ";
        }
        std::cout << std::endl;
    }
}

int main(int argc, char* argv[]) {
    // Usage: test_jc_ladder_permutations <target_cavity> <target_qubit> <moves>
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " <target_cavity> <target_qubit> <moves>" << std::endl;
        return 1;
    }

    int target_cavity = std::atoi(argv[1]);
    int target_qubit = std::atoi(argv[2]);
    int moves = std::atoi(argv[3]);

    step start_step = {0, 0};
    step top_step = {20, 10};
    step target = {target_cavity, target_qubit};

    std::vector<std::vector<step>> perms = jc_ladder_permutations(start_step, top_step, target, moves);

    std::cout << "Number of permutations: " << perms.size() << std::endl;
    print_permutations(perms);

    // Simple check: all permutations should end at target
    bool all_end_at_target = true;
    for (const auto& perm : perms) {
        if (!perm.empty()) {
            const step& last = perm.back();
            if (last.cavity != target.cavity || last.qubit != target.qubit) {
                all_end_at_target = false;
                break;
            }
        }
    }
    if (all_end_at_target) {
        std::cout << "All permutations end at the target step." << std::endl;
    } else {
        std::cout << "Error: Some permutations do not end at the target step!" << std::endl;
    }

    return 0;
}