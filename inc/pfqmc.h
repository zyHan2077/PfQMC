#include "types.h"

class PfQMC
{
private:
    std::vector<MatType> Ul, Ur, Dl, Dr, Vl, Vr;


public:
    PfQMC(int L, std::vector<Operator*> &op_array, int stb);
    ~PfQMC();
};

PfQMC::PfQMC(/* args */)
{
}

PfQMC::~PfQMC()
{
}
