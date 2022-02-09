#include <cmath>

#include "bb.h"

branch_and_bound::branch_and_bound(map_t& map, const window& window, size_t k)
        : _window(window)
        , _map(map)
        , _k(k)
        , _best_suffix_score()
{
    if (_window.empty())
    {
        throw std::runtime_error("The matrix is empty.");
    }

    if (_window.size() != k)
    {
        throw std::runtime_error("The size of the window is not k");
    }

    _result_list.reserve(std::pow(sigma, k));
}

void branch_and_bound::run(score_t omega)
{
    {
        _best_suffix_score.push_back(1.0f);

        // precalc the scores of the best suffixes
        code_t prefix = 0;
        //score_t score = 0.0;
        score_t score = 1.0;
        for (size_t i = 0; i < _k; ++i)
        {
            const auto& [index_best, score_best] = _window.max_at(_k - i - 1);

            prefix = (index_best << i * 2) | prefix;
            //score = score + score_best;
            score = score * score_best;

            //std::cout << "BEST: " << kmer_score << std::endl;
            _best_suffix_score.push_back(score);
        }
    }


    const score_t eps = std::pow((omega / 4), _k);

    /*
    for (size_t i = 0; i < sigma; ++i)
    {
        bb(i, 0, 0, 1.0, eps);
    }*/

    _stack.push_back({std::numeric_limits<code_t>::infinity(), 1.0, 0});

    while (!_stack.empty())
    {
        auto [prefix, score, j] = _stack.back();
        _stack.pop_back();

        if (j == _k)
        {
            _result_list.push_back({prefix, score});
        }
        else
        {
            for (code_t i = 0ul; i < sigma; ++i)
            {
                const auto new_score = score * _window.get(i, j);
                //const auto best_suffix = _best_suffix_score[_k - ((j + 1) + 1)];
                const auto best_suffix = _best_suffix_score[_k - (j + 1)];
                if (new_score * best_suffix > eps)
                {
                    const auto new_prefix = (prefix << 2) | i;
                    _stack.push_back({new_prefix, new_score, static_cast<unsigned short>(j + 1u)});
                }

            }
        }
    }
}
/*
bb_return branch_and_bound::bb(size_t i, size_t j, code_t prefix, score_t score, score_t eps)
{
    // score = score + _matrix[i][j];
    score = score * _window.get(i, j);
    prefix = (prefix << 2) | i;

    if (j == _k - 1)
    {
        if (score > eps)
        {
            _result_list.push_back({prefix, score});
            //_map[prefix] = score;
            return bb_return::GOOD_KMER;
        }
        else
        {
            return bb_return::BAD_PREFIX;
        }
    }

    const auto best_suffix = _best_suffix_score[_k - ((j + 1) + 1)];
    //if (score + best_suffix <= eps)
    if (score * best_suffix <= eps)
    {
        return bb_return::BAD_PREFIX;
    }
    else
    {
        for (size_t i2 = 0; i2 < sigma; ++i2)
        {
            bb(i2, j + 1, prefix, score, eps);
        }
        return bb_return::GOOD_PRFIX;
    }
}*/


const map_t& branch_and_bound::get_map()
{
    return _map;
}

std::vector<bb_return> branch_and_bound::get_returns() const
{
    return _returns;
}

const std::vector<phylo_kmer>& branch_and_bound::get_result() const
{
    return _result_list;
}

size_t branch_and_bound::get_num_kmers() const
{
    return _result_list.size();
}
/*
rappas::rappas(const matrix& matrix, size_t k)
        : _matrix(matrix)
        , _k(k)
{
    if (_matrix.empty())
    {
        throw std::runtime_error("The matrix is empty.");
    }

    if (_matrix[0].size() != k)
    {
        throw std::runtime_error("The size of the window is not k");
    }
}

void rappas::run(score_t omega)
{
    if (!_matrix.is_sorted())
    {
        throw std::runtime_error("Matrix is not sorted.");
    }

    const score_t eps = std::pow((omega / 4), _k);
    //bb(0, 0, 0, 0.0, eps);
    for (size_t i = 0; i < sigma; ++i)
    {
        bb(i, 0, 0, 1.0, eps);
    }
}

bb_return rappas::bb(size_t i, size_t j, code_t prefix, score_t score, score_t eps)
{
    // score = score + _matrix[i][j];
    score = score * _matrix[i][j];
    prefix = (prefix << 2) | _matrix.get_order()[i][j];

    if (j == _k - 1)
    {
        map[prefix] = score;
        return bb_return::GOOD_KMER;
    }
    else
    {
        for (size_t i2 = 0; i2 < sigma; ++i2)
        {
            if (score * _matrix[i2][j+1] <= eps)
            {
                break;
            }
            bb(i2, j + 1, prefix, score, eps);

        }
        return bb_return::GOOD_PRFIX;
    }
}

const map_t& rappas::get_map()
{
    return map;
}
*/