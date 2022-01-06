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

bool kmer_score_comparator(const phylo_kmer& k1, const phylo_kmer& k2)
{
    return k1.score > k2.score;
}


std::vector<phylo_kmer> as_column(const matrix_t& matrix, size_t j, score_t eps)
{
    std::vector<phylo_kmer> column;
    for (size_t i = 0; i < sigma; ++i)
    {
        const auto& row = matrix[i];
        if (row[j] > eps)
        {
            column.push_back({i, row[j] });
        }
    }
    return column;
}

class divide_and_conquer
{
public:
    divide_and_conquer(const matrix_t& matrix, size_t k)
            : _matrix(matrix)
            , _k(k)
    {
        /// kmer_size can also be zero, which means the end() iterator
        const auto halfsize = size_t{ k / 2 };
        _prefix_size = (halfsize >= 1) ? halfsize : k;
    }

    void run()
    {
        const auto kmers = dc(0, _k);

        for (const auto& [kmer, score] : kmers)
        {
            map[kmer] = score;
        }
    }

    // j is the starat position of the window
    // h is the length of the window
    std::vector<phylo_kmer> dc(size_t j, size_t h)
    {
        const score_t eps = std::pow((omega / 4), _k);

        // trivial case
        if (h == 1)
        {
            return as_column(_matrix, j, eps);
        }
        else
        {
            std::vector<phylo_kmer> result;
            const auto l = dc(j, h / 2);
            auto r = dc(j + h / 2, h - h / 2);

            if (!r.empty())
            {
                std::sort(r.begin(), r.end(), kmer_score_comparator);

                for (const auto& [prefix, prefix_score] : l)
                {
                    for (const auto& [suffix, suffix_score] : r)
                    {
                        const auto kmer = (prefix << ((h - h / 2) * 2)) | suffix;
                        //const auto score = prefix_score + suffix_score;
                        const auto score = prefix_score * suffix_score;

                        if (score <= eps)
                        {
                            break;
                        }
                        result.push_back({ kmer, score });
                    }
                }
            }
            return result;
        }
    }

    const map_t& get_map()
    {
        return map;
    }

private:
    const matrix_t& _matrix;
    map_t map;
    size_t _k;
    size_t _prefix_size;

    std::vector<phylo_kmer> _prefixes;
    std::vector<phylo_kmer> _suffixes;
};


#endif //XPAS_ALGS_DAC_H
