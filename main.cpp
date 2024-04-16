#include "pfqmc.h"
#include <omp.h>
int main()
{
    double start_time = omp_get_wtime();
    mkl_set_num_threads(4);
    int Lx = 2;
    int Ly = 2;
    int LTau = 100;
    double dt = 0.1;
    double V = 0.7;
    int stabilizationTime = 10;
    int thermalLength = 200;
    int evaluationLength = 1000;

    // int nDim = Lx * Ly * 4;
    SpinlessTvHoneycombUtils config(Lx, Ly, dt, V, LTau);
    rdGenerator rd(42);
    Honeycomb_tV walker(&config, &rd);
    PfQMC pfqmc(&walker, stabilizationTime);
    for (int i = 0; i < thermalLength; i++)
    {
        pfqmc.rightSweep();
        pfqmc.leftSweep();
        // std::cout << i << std::endl;
    }

    DataType energy = 0.0;
    DataType sign = 1.0;
    for (int i = 0; i < evaluationLength; i++)
    {
        pfqmc.rightSweep();
        pfqmc.leftSweep();
        // sign = pfqmc.getSign();
        energy += sign * config.energyFromGreensFunc(pfqmc.g);

        if ((i % 10) == 0)
            std::cout << "iter = " << i << " energy=" << energy / double(i + 1) << "\n";
    }
    // std::cout << pfqmc.g << "\n";

    double time = omp_get_wtime() - start_time;
    std::cout<<"total time "<<time<<std::endl;
    return 0;
}