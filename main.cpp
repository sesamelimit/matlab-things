#include <iostream>
#include "mh.h"
#include <fstream>

double ** fix_iter_m(double ** const A, size_t const n, double ** const b,double eps){
    if (!matrix_pos_def)
        A_t_A(A,n);

double mu = 1/(matrix_norm(A,n,1));
auto B = E_min_muA(mu,A,n);
print_matrix(B,n,n);
auto c = multMatrix(mu,b,n,1);
double coef = matrix_norm(B,n,1)/(1- matrix_norm(B,n,1));
auto x1 = dyn_array(n,1);
auto x2 = dyn_array(n,1);
Bx_plus_с(B,x1,x2,c,n);
size_t k=1;
while(coef*vector_norm(x1,x2,n)>=eps){
    Bx_plus_с(B,x1,x2,c,n);
    ++k;
}
std::cout << coef*vector_norm(x1,x2,n) << std::endl;

    dyn_array_destroy(B,n);
    dyn_array_destroy(c,n);
    dyn_array_destroy(x1,n);
    print_matrix(x2,n,1);
    dyn_array_destroy(x2,n);
}

int main() {
    size_t n=3;
    auto A = dyn_array(n,n);
   char filename[]="/home/maria/CLionProjects/numerical_methods/data.txt";
    std::ifstream(fin);
    fin.open(filename);
    for(size_t i=0; i<n;++i)
        for(size_t j=0; j<n; ++j)
            fin >> A[j][i];
    print_matrix(A,n,n);
   auto b= dyn_array(n,1);
    for(size_t i=0; i<n;++i)
        fin >> b[i][0];

    fix_iter_m(A,n,b,1e-3);
    dyn_array_destroy(A,n);
    dyn_array_destroy(b,n);
    fin.close();
    return 0;
}
