# Pfaffian Quantum Monte Carlo (PfQMC)

usage:

```bash
# setvars for mkl first, 
# mkl dependencies should be automatically located
mkdir build && cd build
cmake .. -DEIGEN3_INCLUDE_DIR=/path/to/eigen3 & make
```

or:

```bash
mkdir obj && mkdir bin
make EIGEN3_INCLUDE_DIR=/path/to/eigen3
```

## Program Design

```mermaid
  graph LR;
    SpinlessTvUtils--inherited by-->SpinlessTvSquareUtils-.used by.->Square_tV;

    SpinlessTvUtils--> SpinlessTvHoneycombUtils-.->Honeycomb_tV;
    SpinlessTvUtils-.->SpinlessVOperator;
    Operator-.->Spinless_tV-..-> PfQMC;
    Operator-->SpinlessVOperator;
    Spinless_tV-->Square_tV;
    Spinless_tV-->Honeycomb_tV;
    SpinlessVOperator;
    Operator-->DenseOperator;
    
```