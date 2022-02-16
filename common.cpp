#include <cmath>
#include "common.h"


bool kmer_score_comparator(const phylo_kmer& k1, const phylo_kmer& k2)
{
    return k1.score > k2.score;
}

score_t get_threshold(score_t omega, size_t k)
{
    return std::pow((omega / 4), k);
}