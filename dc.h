#ifndef XPAS_ALGS_DAC_H
#define XPAS_ALGS_DAC_H

#include <algorithm>
#include <cmath>
#include "common.h"
#include "matrix.h"

struct phylo_kmer
{
    code_t kmer;
    score_t score;
};

class divide_and_conquer
{
public:
    divide_and_conquer(map_t& map, const window& window, size_t k);
    void run(score_t omega);
    std::vector<phylo_kmer> dc(score_t omega, size_t j, size_t h);
    const map_t& get_map();

private:
    const window& _window;
    map_t& _map;
    size_t _k;
    size_t _prefix_size;

    std::vector<phylo_kmer> _prefixes;
    std::vector<phylo_kmer> _suffixes;
};


#endif //XPAS_ALGS_DAC_H
