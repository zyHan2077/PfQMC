#include "operator.h"
#include "honeycomb.h"

SpinlessVOperator::SpinlessVOperator(const SpinlessTvHoneycombUtils* _config, iVecType* _s, int _bondType, rdGenerator* _rd)
    :config(_config), bondType(_bondType), etaM(_config->etaM), 
        nDim(config->nsites * 2), nUnitcell(config->nUnitcell) {

    s = _s;
    rd = _rd;

    B = MatType::Identity(nDim, nDim);
    config->InteractionBGenerator(B, *s, bondType);
}

void SpinlessVOperator::reCalc() {
    B = MatType::Identity(nDim, nDim);
    config->InteractionBGenerator(B, *s, bondType);
}

bool SpinlessVOperator::singleFlip(MatType &g, int idxCell, double rand) {
    DataType r = etaM;
    auto m = config->idxCell2Coord(idxCell);
    int auxCur = (*s)(idxCell);
    int idx1, idx2;
    DataType tmp[2];
    const int inc =  1;
    DataType alpha;
    

    for (int imaj = 0; imaj < 2; imaj ++) {
        idx1 = config->majoranaCoord2Idx(m.ix, m.iy, 0, imaj);
        idx2 = config->neighborSiteIdx(m.ix, m.iy, imaj, bondType);
        // tmp = [1 + i \sigma_{12} \tanh(\lambda / 2) G_{12}]
        tmp[imaj] = ( 1.0 - ( (1.0i) * (config->thlV) * double(auxCur) * g(idx1, idx2) ) );
        r *= tmp[imaj];
    }

    bool flag = rand < std::abs(r);

    if (flag) {
        for (int imaj = 0; imaj < 2; imaj ++) {
            idx1 = config->majoranaCoord2Idx(m.ix, m.iy, 0, imaj);
            idx2 = config->neighborSiteIdx(m.ix, m.iy, imaj, bondType);

            // update aux field and B matrix
            (*s)(idxCell) = -auxCur;
            B(idx1, idx2) = -B(idx1, idx2);
            B(idx2, idx1) = -B(idx2, idx1);

            // update Green's function
            cVecType x1 = B.col(idx1);
            cVecType x2 = B.col(idx2);
            alpha = (-1.0i) * (config->thlV) / tmp[imaj];
            zgeru(&nDim, &nDim, &alpha, x1.data(), &inc, x2.data(), &inc, B.data(), &nDim);
            alpha = -alpha;
            zgeru(&nDim, &nDim, &alpha, x2.data(), &inc, x1.data(), &inc, B.data(), &nDim);
        }
    }

    return flag;
}

void SpinlessVOperator::update(MatType &g) {
    for (int i=0; i<nUnitcell; i++) {
        // random real between (0, 1)
        double rand = rd->rdUniform01();
        bool flag = singleFlip(g, i, rand);
    }
}

void SpinlessVOperator::right_multiply(const MatType &AIn, MatType &AOut) {
    AOut = AIn * B;
}