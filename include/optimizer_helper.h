#include <iostream>
#include <vector>
#include <cmath>
#include <dlib/dlib/optimization.h>
// #include "divdiffcomplex.h"
// #include "permutation.h"


using namespace std;
using namespace dlib;



complex_Ex operator+(const complex_Ex a, const complex_Ex b) {
    divdiff_init();
    complex_Ex res;

    ExExFloat temp1 = a.real + b.real;
    ExExFloat temp2 = b.imag + b.imag;

    res.real = temp1;
    res.imag = temp2;
    divdiff_clear_up();

    return res;
}

complex_Ex operator*(const complex_Ex a, const complex_Ex b) {
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

complex_Ex operator*(const complex_Ex a, double x) {
    divdiff_init();
    complex_Ex res;

    res.real = a.real * x;
    res.imag = a.imag * x;
    divdiff_clear_up();

    return res;
}

double abs2(const complex_Ex z) {
    divdiff_init();
    ExExFloat res1 = z.real * z.real
    ExExFloat res2 = z.imag * z.imag;
    ExExFloat res = res1 + res2;
    divdiff_clear_up();
    return res.get_double();
}

complex_Ex conj(const complex_Ex z) {
    divdiff_init();
    complex_Ex res;

    res.real = z.real;
    res.imag = -z.imag;
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

    double operator()(const matrix<double, 0, 1>& x) const {
        mat_c M = identity_matrix(C[0].size());

        for (long i = 0; i < x.size(); ++i)
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