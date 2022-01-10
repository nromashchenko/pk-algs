#ifndef XPAS_ALGS_COMMON_H
#define XPAS_ALGS_COMMON_H

#include <vector>
#include <unordered_map>

using score_t = float;
using code_t = uint64_t;
using row_t = std::vector<score_t>;
using matrix_t = std::vector<row_t>;
using map_t = std::unordered_map<code_t, score_t>;

static const size_t sigma = 4;

std::pair<size_t, score_t> max_at(const matrix_t& matrix, size_t column);


enum class bb_return
{
    BAD_PREFIX = 0,
    GOOD_KMER = 1,
    GOOD_PRFIX = 2
};


#endif //XPAS_ALGS_COMMON_H
