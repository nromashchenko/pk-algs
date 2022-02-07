#ifndef XPAS_ALGS_DAC_H
#define XPAS_ALGS_DAC_H

#include <algorithm>
#include <cmath>
#include "common.h"
#include "matrix.h"

class divide_and_conquer
{
public:
    divide_and_conquer(map_t& map, const window& window, size_t k);
    void run(score_t omega);

    const map_t& get_map();

    const std::vector<phylo_kmer>& get_result() const;

    size_t get_num_kmers() const;

private:
    std::vector<phylo_kmer> dc(score_t omega, size_t j, size_t h, score_t eps);

    void preprocess();

    score_t best_score(size_t j, size_t h);

    const window& _window;
    map_t& _map;
    size_t _k;
    size_t _prefix_size;

    std::vector<phylo_kmer> _prefixes;
    std::vector<phylo_kmer> _suffixes;

    std::vector<score_t> _best_scores;

    std::vector<phylo_kmer> _result_list;
};


#endif //XPAS_ALGS_DAC_H
