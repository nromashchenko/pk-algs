#include <cmath>
#include "brute_force.h"

brute_force::brute_force(const window& window, size_t k)
        : _window(window)
        , _k(k)
{
}

void brute_force::run(score_t omega)
{
    const score_t eps = std::pow((omega / 4), _k);

    // generate all k-mers
    for (size_t i = 0; i < sigma; ++i)
    {
        bf(i, 0, 0, 1.0, eps);
    }

    // select ones that have the score > eps
    map_t new_map;
    for (const auto& [code, score] : map)
    {
        if (score > eps)
        {
            new_map[code] = score;
        }
    }
    map = new_map;

}

bb_return brute_force::bf(size_t i, size_t j, code_t prefix, score_t score, score_t eps)
{
    // score = score + _matrix[i][j];
    score = score * _window.get(i, j);
    prefix = (prefix << 2) | i;

    if (j == _k - 1)
    {
        map[prefix] = score;
        return bb_return::GOOD_KMER;
    }
    else
    {
        for (size_t i2 = 0; i2 < sigma; ++i2)
        {
            bf(i2, j + 1, prefix, score, eps);
        }
        return bb_return::GOOD_PRFIX;
    }
}

const map_t& brute_force::get_map()
{
    return map;
}
