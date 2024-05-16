#include <omp.h>
#include <string.h>

#include "inc/honeycomb.h"
#include "inc/pfqmc.h"
#include "inc/square.h"
#include "inc/singleMajoranaHoneycomb.h"

int main_honeycomb() {
    double start_time = omp_get_wtime();
    mkl_set_num_threads(8);

    int Lx = 4;
    int Ly = 4;
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
    for (int i = 0; i < thermalLength; i++) {
        pfqmc.rightSweep();
        pfqmc.leftSweep();
        // std::cout << i << std::endl;
    }

    DataType energy = 0.0;
    DataType sign = 1.0;
    for (int i = 0; i < evaluationLength; i++) {
        pfqmc.rightSweep();
        pfqmc.leftSweep();
        // check if the sign problem free condition is met
        if (std::abs(1.0 - pfqmc.sign) > 1e-2) {
            std::cout << "\n=== error in sign at round = " << i << " sign = " << pfqmc.sign << " ≠ " << 1.0 << "==== \n"; 
        }
        // sign = pfqmc.getSignRaw();
        energy += sign * config.energyFromGreensFunc(pfqmc.g);

        if ((i % 10) == 0)
            std::cout << "iter = " << i << " energy=" << energy / double(i + 1)
                      << "\n";
    }
    // std::cout << pfqmc.g << "\n";

    double time = omp_get_wtime() - start_time;
    std::cout << "total time " << time << std::endl;
    return 0;
}

int main_honeycombSingleMajorana(int Lx, int Ly, int LTau, double dt, double V, int nthreads) {
    double start_time = omp_get_wtime();
    mkl_set_num_threads(nthreads);

    // int Lx = 10;
    // int Ly = 10;
    // int LTau = 33;
    // double dt = 0.1;
    // double V = 1.65;
    int stabilizationTime = 10;
    int thermalLength = 200;
    int evaluationLength = 1000;
    std::cout << "=== Honeycomb Single Majorana Model ===\n";
    std::cout << "Lx = " << Lx << " Ly = " << Ly << " LTau = " << LTau << " dt = " << dt << " V = " << V << std::endl;
    SpinlessTvHoneycombSingleMajoranaUtils config(Lx, Ly, dt, V, LTau);
    rdGenerator rd(42);
    HoneycombSingleMajorana_tV walker(&config, &rd);
    PfQMC pfqmc(&walker, stabilizationTime);
    for (int i = 0; i < thermalLength; i++) {
	    std::cout << i << " ";
        pfqmc.rightSweep();
        pfqmc.leftSweep();
        // std::cout << i << std::endl;
    }
    std::cout << std::endl;

    DataType energy = 0.0;
    DataType structureFactorCDW = 0.0;
    DataType sign, signRaw;
    DataType srSignTot = 0.0;
    for (int i = 0; i < evaluationLength; i++) {
        pfqmc.rightSweep();
        pfqmc.leftSweep();
        sign = pfqmc.sign;

        if (i % 20 == 0) {
            signRaw = pfqmc.getSignRaw();
            if (std::abs(sign - signRaw) > 1e-2) {
                std::cout << "\n=== error in sign at round = " << i << " sign = " << sign << " ≠ " << signRaw << "==== \n"; 
            }
	        // pfqmc.sign = signRaw;
        }

        // srSignTot += sign;
        energy = config.energyFromGreensFunc(pfqmc.g);
        structureFactorCDW = config.structureFactorCDW(pfqmc.g);
        std::cout << "iter = " << i << " energy = " << energy << " structureFactorCDW = " << structureFactorCDW << " MRsign = " << pfqmc.sign << "\n";

        // if ((i % 10) == 0)
        //     std::cout << "iter = " << i << " energy=" << energy / double(i + 1) << " srSign = " << srSignTot / double(i + 1) << " structureFactorCDW = " << structureFactorCDW / double(i + 1)
        //               << "\n";
    }
    // std::cout << pfqmc.g << "\n";

    double time = omp_get_wtime() - start_time;
    std::cout << "=== End of Honeycomb Single Majorana Model ===\n";
    std::cout << "total time " << time << std::endl;
    return 0;
}

int main_square() {
    double start_time = omp_get_wtime();
    mkl_set_num_threads(8);

    int Lx = 4;
    int Ly = 4;
    int LTau = 200;
    double dt = 0.1;
    double V = 0.7;
    double delta = 0.5;
    int stabilizationTime = 10;
    int thermalLength = 200;
    int evaluationLength = 1000;

    // int nDim = Lx * Ly * 4;
    SpinlessTvSquareUtils config(Lx, Ly, dt, V, LTau, delta);
    rdGenerator rd(42);
    Square_tV walker(&config, &rd);
    PfQMC pfqmc(&walker, stabilizationTime);
    // std::cout << "here!" << std::endl;
    for (int i = 0; i < thermalLength; i++) {
        pfqmc.rightSweep();
        pfqmc.leftSweep();
    }

    DataType energy = 0.0;
    DataType sign = 1.0;
    DataType signFast;
    DataType signTot = 0.0;
    DataType rawSignTot = 0.0;
    DataType rawSign = 1.0;
    // pfqmc.sign = 1.0;
    for (int i = 0; i < evaluationLength; i++) {
        pfqmc.rightSweep();
        pfqmc.leftSweep();
        sign = pfqmc.sign;

        // rawSign = pfqmc.getSignRaw();
        rawSignTot += rawSign;
        // double deviation = std::abs(sign - rawSign);
        // if (deviation > 1e-2) {
        //     pfqmc.sign = rawSign; // update sign
        //     std::cout << "\n=== error in sign at round = " << i << " raw sign = " << rawSign << " ≠ " << sign << "==== \n"; 
        // }

        energy += sign * config.energyFromGreensFunc(pfqmc.g);
        signTot += sign;
        // std::cout << "sign = " << sign << " raw sign = " << rawSign << "\n";
        if ((i % 10) == 0)
            std::cout << "iter = " << i << " sign ave = " << signTot / double(i + 1) << " raw sign ave= " << rawSignTot / double(i + 1) 
                      << " energy=" << energy / signTot << " curSign = "  <<  sign << "\n";
    }
    // std::cout << pfqmc.g << "\n";

    double time = omp_get_wtime() - start_time;
    std::cout << "total time " << time << std::endl;
    return 0;
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cout << "Usage: " << argv[0] << " --square | --honeycomb | --SRhoneycomb\n";
        return 0;
    }

    if (std::strcmp(argv[1], "--square") == 0) {
        return main_square();
    } else if (std::strcmp(argv[1], "--honeycomb") == 0) {
        return main_honeycomb();
    } else if (std::strcmp(argv[1], "--SRhoneycomb") == 0) {
        int Lx, Ly, LTau, nthreads;
        double dt, V;
        Lx = std::stoi(argv[2]);
        Ly = std::stoi(argv[3]);
        LTau = std::stoi(argv[4]);
        dt = std::stod(argv[5]);
        V = std::stod(argv[6]);
        nthreads = std::stoi(argv[7]);
        // fstream fin()
        return main_honeycombSingleMajorana(Lx, Ly, LTau, dt, V, nthreads);
    } else {
        std::cout << argv[1] << ": invalid arguments\n";
        return 0;
    }
}
