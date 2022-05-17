#include <cmath>

#include "bb.h"

branch_and_bound::branch_and_bound(const window& window, size_t k)
        : _window(window)
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

    preprocess();
}

void branch_and_bound::run(score_t omega)
{
    const score_t eps = get_threshold(omega, _k);

    // Recursive BB
    for (size_t i = 0; i < sigma; ++i)
    {
        bb(i, 0, 0, 1.0, eps);
    }

    // Iterative BB
    /*
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
    }*/
}

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

    const auto best_suffix = _best_suffix_score[_k - (j + 2)];
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
}

/*
const map_t& branch_and_bound::get_map()
{
    return _map;
}
*/
std::vector<bb_return> branch_and_bound::get_returns() const
{
    return _returns;
}

const std::vector<phylo_kmer>& branch_and_bound::get_result() const
{
    return _result_list;
}

std::vector<phylo_kmer>&& branch_and_bound::get_result()
{
    return std::move(_result_list);
}

size_t branch_and_bound::get_num_kmers() const
{
    return _result_list.size();
}

void branch_and_bound::preprocess()
{
    // For the iterative BB
    //_best_suffix_score.push_back(1.0f);

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



bbe::bbe(const window& window, std::vector<column_data> order, size_t k)
    : _window(window), _order(std::move(order)), _k(k), _best_suffix_score(k)
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

    preprocess();
}

void bbe::run(score_t omega)
{
    const score_t eps = get_threshold(omega, _k);

    // Recursive BB
    for (size_t i = 0; i < sigma; ++i)
    {
        bb(i, 0, 0, 1.0, eps);
    }
}

bb_return bbe::bb(size_t i, size_t column_id, code_t prefix, score_t score, score_t eps)
{
    const size_t j = _order[column_id].j;
    score = score * _window.get(i, j);
    //prefix = (prefix << 2) | i;
    prefix |= i << (2 * (_k - 1 - j));

    if (column_id == _k - 1)
    {
        if (score > eps)
        {
            _result_list.push_back({prefix, score});
            return bb_return::GOOD_KMER;
        }
        else
        {
            return bb_return::BAD_PREFIX;
        }
    }

    const auto best_suffix = _best_suffix_score[column_id + 1];
    if (score * best_suffix <= eps)
    {
        return bb_return::BAD_PREFIX;
    }
    else
    {
        for (size_t i2 = 0; i2 < sigma; ++i2)
        {
            bb(i2, column_id + 1, prefix, score, eps);
        }
        return bb_return::GOOD_PRFIX;
    }
}

std::vector<bb_return> bbe::get_returns() const
{
    return _returns;
}

const std::vector<phylo_kmer>& bbe::get_result() const
{
    return _result_list;
}

size_t bbe::get_num_kmers() const
{
    return _result_list.size();
}

void bbe::preprocess()
{
    // Precalculate the scores of the best suffixes, in the reverse column order given by the heap
    score_t score = 1.0;

    for (int column_id = _k - 1; column_id >= 0; --column_id)
    {
        const auto j = _order[column_id].j;
        const auto& [index_best, score_best] = _window.max_at(j);

        score = score * score_best;
        _best_suffix_score[column_id] = score;
    }
}



baseline::baseline(const window& window, size_t k, size_t num_kmers)
    : _window(window), _k(k), _num_kmers(num_kmers)
{
    _result_list.reserve(num_kmers);
}

void baseline::run(score_t omega)
{
    //size_t output_size = (rand() % (size_t)std::pow(sigma, _k)) + 1;
    const auto eps = get_threshold(omega, _k);
    for (size_t i = 0; i < _num_kmers; ++i)
    {
        _result_list.push_back({i, eps});
    }
}

const std::vector<phylo_kmer>& baseline::get_result() const
{
    return _result_list;
}

size_t baseline::get_num_kmers() const
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