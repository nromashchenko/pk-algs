#include <iostream>
#include <unordered_map>
#include <vector>
#include <cmath>
#include <algorithm>
#include <random>
#include <iomanip>
#include <cassert>

#include "common.h"
#include "dac.h"


enum class bb_return
{
    BAD_PREFIX = 0,
    GOOD_KMER = 1,
    GOOD_PRFIX = 2
};

std::pair<size_t, score_t> max_at(const matrix_t& matrix, size_t column)
{
    size_t max_index = 0;
    score_t max_score = matrix[0][column];
    for (size_t i = 1; i < matrix.size(); ++i)
    {
        if (matrix[i][column] > max_score)
        {
            max_score = matrix[i][column];
            max_index = i;
        }
    }
    return { max_index, max_score };
}

class brute_force
{
public:
    brute_force(const matrix_t& matrix, size_t k)
        : _matrix(matrix)
        , _k(k)
    {
    }

    void run()
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

    bb_return bf(size_t i, size_t j, code_t prefix, score_t score, score_t eps)
    {
        // score = score + _matrix[i][j];
        score = score * _matrix[i][j];
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

    const map_t& get_map()
    {
        return map;
    }


private:
    const matrix_t& _matrix;
    map_t map;
    size_t _k;
};

class branch_and_bound
{
public:
    branch_and_bound(const matrix_t& matrix, size_t k)
        : _matrix(matrix)
        , _k(k)
        , _best_suffix_score()
    {
        if (_matrix.empty())
        {
            throw std::runtime_error("The matrix is empty.");
        }

        if (_matrix[0].size() != k)
        {
            throw std::runtime_error("The size of the window is not k");
        }
        // precalc the scores of the best suffixes

        code_t prefix = 0;
        //score_t score = 0.0;
        score_t score = 1.0;
        for (size_t i = 0; i < _k; ++i)
        {
            const auto& [index_best, score_best] = max_at(_matrix, _k - i - 1);

            prefix = (index_best << i * 2) | prefix;
            //score = score + score_best;
            score = score * score_best;

            //std::cout << "BEST: " << kmer_score << std::endl;
            _best_suffix_score.push_back(score);
        }
    }

    void run()
    {
        const score_t eps = std::pow((omega / 4), _k);
        //bb(0, 0, 0, 0.0, eps);
        for (size_t i = 0; i < sigma; ++i)
        {
            bb(i, 0, 0, 1.0, eps);
        }
    }

    bb_return bb(size_t i, size_t j, code_t prefix, score_t score, score_t eps)
    {
        // score = score + _matrix[i][j];
        score = score * _matrix[i][j];
        prefix = (prefix << 2) | i;

        if (j == _k - 1 && score > eps)
        {
            map[prefix] = score;
            /*if (prefix == 16)
            {
                std::cout << "16" << std::endl;
            }*/
            return bb_return::GOOD_KMER;
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
    }

    const map_t& get_map()
    {
        return map;
    }

private:
    const matrix_t& _matrix;
    map_t map;
    size_t _k;
    std::vector<score_t> _best_suffix_score;
};



std::random_device rd;
std::default_random_engine eng(42);
std::uniform_real_distribution<score_t> distr(0, 1);

score_t g()
{
    return distr(eng);
}

row_t generate_row(size_t k)
{
    std::vector<score_t> row(k);
    std::generate(row.begin(), row.end(), g);
    return row;
}

matrix_t generate(size_t k)
{
    matrix_t a(sigma);
    for (size_t i = 0; i < a.size(); ++i)
    {
        a[i] = generate_row(k);
    }

    for (size_t j = 0; j < k; ++j)
    {
        score_t sum = 0.0;
        for (const auto& row : a)
        {
            sum += row[j];
        }
        for (auto& row : a)
        {
            row[j] = row[j] / sum;
        }
    }
    return a;
}

void print_matrix(const matrix_t& matrix)
{
    for (const auto& row : matrix)
    {
        for (const auto& el : row)
        {
            std::cout << std::fixed << std::setprecision(8) << el << "\t";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

void print_map(const map_t& map)
{
    for (const auto& [kmer, score] : map)
    {
        std::cout << kmer << ": " << score << std::endl;
    }
    std::cout << std::endl;
}

void assert_equal(const map_t& map1, const map_t& map2)
{
    assert(map1.size() == map2.size());
    for (const auto& [kmer, score] : map1)
    {
        assert(map2.find(kmer) != map2.end());

        assert(fabs(score - map2.at(kmer)) < 1e-6);
    }
}

void test_one(size_t k, bool print=true)
{
    const auto matrix = generate(k);
    if (print)
    {
        print_matrix(matrix);
        std::cout << "Threshold: " << std::pow((omega / 4), k) << std::endl;
    }

    branch_and_bound bb(matrix, k);
    bb.run();
    if (print)
    {
        //print_map(bb.get_map());
        std::cout << "Branch-and-bound, generated: " << bb.get_map().size() << std::endl;
    }

    divide_and_conquer dc(matrix, k);
    dc.run();
    if (print)
    {
        //print_map(dc.get_map());
        std::cout << "Divide-and-conquer, generated: " << dc.get_map().size() << std::endl;
    }

    brute_force bf(matrix, k);
    bf.run();
    if (print)
    {
        //print_map(bf.get_map());
        std::cout << "Brute force, generated: " << bf.get_map().size() << std::endl;
    }

    assert_equal(bb.get_map(), bf.get_map());
    assert_equal(dc.get_map(), bf.get_map());
}

void test_suite()
{
    const size_t num_iter = 100;
    const std::vector<size_t> k_values = { 6, 7, 8, 9 };

    for (const auto k : k_values)
    {
        for (size_t i = 0; i < num_iter; ++i)
        {
            std::cout << "\rTesting k = " << k << ". " << i << " / " << num_iter << "..." << std::flush;
            test_one(k, false);
        }
        std::cout << "\rTesting k = " << k << ". Done." << std::endl;
    }
}

int main()
{
    srand(42);

    /*{
        test_one(10);
    }*/

    test_suite();
    return 0;
}
