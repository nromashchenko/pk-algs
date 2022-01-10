#ifndef XPAS_ALGS_DAC_H
#define XPAS_ALGS_DAC_H

#include <algorithm>
#include <cmath>
#include "common.h"

struct phylo_kmer
{
    code_t kmer;
    score_t score;
};

class divide_and_conquer
{
public:
    divide_and_conquer(const matrix_t& matrix, size_t k);
    void run();
    std::vector<phylo_kmer> dc(size_t j, size_t h);
    const map_t& get_map();

private:
    const matrix_t& _matrix;
    map_t map;
    size_t _k;
    size_t _prefix_size;

    std::vector<phylo_kmer> _prefixes;
    std::vector<phylo_kmer> _suffixes;
};


#endif //XPAS_ALGS_DAC_H
