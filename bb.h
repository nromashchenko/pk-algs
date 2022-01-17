#ifndef XPAS_ALGS_BB_H
#define XPAS_ALGS_BB_H

#include "common.h"
#include "matrix.h"

class branch_and_bound
{
public:
    branch_and_bound(auto window, size_t k);
    void run(score_t omega);
    bb_return bb(size_t i, size_t j, code_t prefix, score_t score, score_t eps);
    const map_t& get_map();

private:
    auto window;
    map_t map;
    size_t _k;
    std::vector<score_t> _best_suffix_score;
};

class rappas
{
public:
    rappas(const matrix& matrix, size_t k);
    void run(score_t omega);
    bb_return bb(size_t i, size_t j, code_t prefix, score_t score, score_t eps);
    const map_t& get_map();

private:
    const matrix& _matrix;
    map_t map;
    size_t _k;
};


#endif //XPAS_ALGS_BB_H
