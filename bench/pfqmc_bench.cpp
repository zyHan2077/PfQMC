#include "pfqmc.h"
#include "singleMajoranaHoneycomb.h"
#include "honeycomb.h"
#include "square.h"
#include "kitaevChain.h"
#include <chrono>
#include <iomanip>
#include <memory>

using Clock = std::chrono::steady_clock;
double seconds(Clock::time_point t) {
    return std::chrono::duration<double>(Clock::now() - t).count();
}
int main(int argc, char** argv) {
    if (argc != 8) {
        std::cerr << "usage: pfqmc_bench MODEL L Ltau pairs threads seed check\n"
                  << "MODEL: majorana, honeycomb, square, chain0, chain1\n";
        return 2;
    }
    std::string model = argv[1];
    int L = std::stoi(argv[2]), lt = std::stoi(argv[3]);
    int pairs = std::stoi(argv[4]), threads = std::stoi(argv[5]);
    int seed = std::stoi(argv[6]);
    bool check = std::stoi(argv[7]) != 0;
    if (L < 2 || lt < 1 || pairs < 1 || threads < 1) return 2;
    mkl_set_dynamic(0);
    mkl_set_num_threads(threads);
    rdGenerator rd(seed);
    SpinlessTvHoneycombSingleMajoranaUtils cm(L, L, .1, .7, lt);
    SpinlessTvHoneycombUtils ch(L, L, .1, .7, lt);
    SpinlessTvSquareUtils cs(L, L, .1, .7, lt);
    SpinlessTvChainUtils cc(L, .1, .7, lt, 0, .3, .2, model == "chain1");
    auto totalStart = Clock::now();
    std::shared_ptr<Spinless_tV> walker;
    int sites;
    if (model == "majorana") {
        walker.reset(new HoneycombSingleMajorana_tV(&cm, &rd)); sites = cm.nsites;
    } else if (model == "honeycomb") {
        walker.reset(new Honeycomb_tV(&ch, &rd)); sites = ch.nsites;
    } else if (model == "square") {
        walker.reset(new Square_tV(&cs, &rd)); sites = cs.nsites;
    } else if (model == "chain0" || model == "chain1") {
        walker.reset(new Chain_tV(&cc, &rd)); sites = cc.nsites;
    } else return 2;
    PfQMC q(walker.get(), 10);
    double init = seconds(totalStart);
    q.rightSweep(); q.leftSweep();
    auto t = Clock::now();
    for (int i = 0; i < pairs; ++i) { q.rightSweep(); q.leftSweep(); }
    double sweeps = seconds(t);
    double total = seconds(totalStart);
    double rawError = 0, rebuildError = 0;
    if (check) {
        rawError = std::abs(q.sign - q.getSignRaw());
        MatType saved = q.g;
        q.leftInit();
        rebuildError = (q.g - saved).norm() / std::max(1.0, saved.norm());
    }
    double checksum = 0;
    for (int j = 0; j < q.nDim; ++j)
        for (int i = 0; i < q.nDim; ++i)
            checksum += q.g(i,j).real() * (1 + (i + 3*j) % 17)
                      + q.g(i,j).imag() * (1 + (3*i+j) % 19);
    long long fields = 0;
    for (auto op : walker->op_array) if (auto s = op->getAuxField())
        for (int i = 0; i < s->size(); ++i) fields += (*s)(i)*(i+1);
    std::cout << std::setprecision(17)
              << model << ',' << sites << ',' << q.nDim << ',' << lt << ','
              << pairs << ',' << threads << ',' << seed << ',' << init << ','
              << sweeps << ',' << total << ',' << q.sign.real() << ',' << q.sign.imag()
              << ',' << checksum << ',' << fields << ',' << rawError << ',' << rebuildError << '\n';
    if (!q.g.allFinite() || !std::isfinite(rawError) || !std::isfinite(rebuildError)
        || rawError > 1e-7 || rebuildError > 1e-8) return 1;
}
