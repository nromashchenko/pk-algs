#ifndef XPAS_ALGS_DAC_H
#define XPAS_ALGS_DAC_H

#include <algorithm>
#include <cmath>
#include "common.h"
#include "matrix.h"

class dccw;

class divide_and_conquer
{
    friend class dccw;
public:
    divide_and_conquer(const window& window, size_t k, score_t omega);
    void run(score_t omega);

    const map_t& get_map();

    const std::vector<phylo_kmer>& get_result() const;

    size_t get_num_kmers() const;

    void preprocess();

    std::vector<phylo_kmer> dc(score_t omega, size_t j, size_t h, score_t eps);
private:


    score_t best_score(size_t j, size_t h);

    const window& _window;
    size_t _k;
    size_t _prefix_size;

    std::vector<phylo_kmer> _prefixes;
    std::vector<phylo_kmer> _suffixes;

    std::vector<score_t> _best_scores;

    std::vector<phylo_kmer> _result_list;
};

class dccw
{
public:
    dccw(const window& window, std::vector<phylo_kmer>& prefixes, size_t k, score_t lookbehind, score_t lookahead,
         score_t omega);
    void run(score_t omega);

    const map_t& get_map();

    const std::vector<phylo_kmer>& get_result() const;

    size_t get_num_kmers() const;

    std::vector<phylo_kmer>&& get_suffixes();

    score_t get_best_suffix_score() const;

private:
    std::vector<phylo_kmer> dc(score_t omega, size_t j, size_t h, score_t eps);

    //void preprocess();

    //score_t best_score(size_t j, size_t h);

    const window& _window;
    size_t _k;
    size_t _prefix_size;

    // The second score bound for suffixes: the best suffix score of the next window
    score_t _lookahead;

    // The second score bound for prefixes: the best prefix score of the previous window
    score_t _lookbehind;

    std::vector<phylo_kmer>& _prefixes;
    std::vector<phylo_kmer> _suffixes;

    //std::vector<score_t> _best_scores;

    std::vector<phylo_kmer> _result_list;


    divide_and_conquer _dc;
};


#endif //XPAS_ALGS_DAC_H
