#include "common.h"


bool kmer_score_comparator(const phylo_kmer& k1, const phylo_kmer& k2)
{
    return k1.score > k2.score;
}
