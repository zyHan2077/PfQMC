#ifndef CATALOG_CHECKS_H
#define CATALOG_CHECKS_H
#include "catalog/catalog.h"
#include <ostream>
namespace catalog {
void closeCheck(DataType a,DataType b,const std::string& label,double tol=2e-7);
struct Recorder {
    std::ostream& out;
    long rows=0,positive=0,zeros=0,negative=0;
    double maxError=0.,maxPhaseError=0.;
    explicit Recorder(std::ostream& stream);
    void record(const std::string& id,const std::string& kind,int index,int modes,DataType z,DataType ref,
                int expected=1,long attempts=0,long accepted=0);
};
void verifyCore(Recorder& r,rdGenerator& rng);
void verifySpecials(Recorder& r,rdGenerator& rng,int samples);
DataType independentTrace(const Hamiltonian& h,int slices,double beta,const std::vector<double>& fields);
void verifyModel(Recorder& r,rdGenerator& rng,const std::string& name,int samples,int sweeps);
}
#endif
