#ifndef XPAS_ALGS_BRUTE_FORCE_H
#define XPAS_ALGS_BRUTE_FORCE_H

#include "common.h"
#include "matrix.h"

class brute_force
{
public:
    brute_force(const window& window, size_t k);
    void run(score_t omega);
    bb_return bf(size_t i, size_t j, code_t prefix, score_t score, score_t eps);
    const map_t& get_map();

private:
    const window& _window;
    map_t map;
    size_t _k;
};


#endif //XPAS_ALGS_BRUTE_FORCE_H
