#include "inc/pfqmc.h"

int main() {
    
    mkl_set_num_threads(8);

    int Lx = 9;
    int Ly = 11;
    int LTau = 20;
    double dt = 0.1;
    double V = 0.7;

    int nDim = Lx * Ly * 4;
    SpinlessTvHoneycombUtils config(Lx, Ly, dt, V, LTau);
    
}