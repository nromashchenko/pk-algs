#include "matrix.h"

#include <utility>
#include <random>
#include <algorithm>
#include <iostream>
#include <iomanip>

std::vector<score_t> get_column(const matrix& matrix, size_t j)
{
    std::vector<score_t> column;
    for (size_t i = 0; i < sigma; ++i)
    {
        const auto& row = matrix[i];
        column.push_back(row[j]);
    }
    return column;
}

matrix::matrix(std::vector<row> d)
    : data(std::move(d))
    , sorted(false)
{
}

std::vector<score_t> matrix::operator[](size_t j) const
{
    return data[j];
}

size_t matrix::size() const
{
    return data.size();
}

bool matrix::empty() const
{
    return data.empty();
}

void matrix::sort()
{
    sorted = true;

    if (order.empty())
    {
        const size_t length = data[0].size();
        order = { sigma, std::vector<size_t>(length, sigma + 1)};
    }


    for (size_t j = 0; j < data[0].size(); ++j)
    {
        // data are stored in rows, extract the column
        auto column = get_column(data, j);

        // sort it by score
        struct score_pair
        {
            score_t score;
            size_t order;
        };
        std::vector<score_pair> pairs(sigma);
        for (size_t i = 0; i < sigma; ++i)
        {
            pairs[i] = { column[i], i };
        }

        auto compare = [](const score_pair& p1, const score_pair& p2) { return p1.score > p2.score; };
        std::sort(begin(pairs), end(pairs), compare);

        // save the sorting
        for (size_t i = 0; i < sigma; ++i)
        {
            data[i][j] = pairs[i].score;
            order[i][j] = pairs[i].order;
        }
    }
}

bool matrix::is_sorted() const
{
    return sorted;
}

std::pair<size_t, score_t> matrix::max_at(size_t column) const
{
    size_t max_index = 0;
    score_t max_score = data[0][column];
    for (size_t i = 1; i < data.size(); ++i)
    {
        if (data[i][column] > max_score)
        {
            max_score = data[i][column];
            max_index = i;
        }
    }
    return { max_index, max_score };
}

std::vector<matrix::row> matrix::get_data() const
{
    return data;
}

std::vector<std::vector<size_t>> matrix::get_order() const
{
    return order;
}

score_t g()
{
    return distr(eng);
}

matrix::row generate_row(size_t k)
{
    std::vector<score_t> row(k);
    std::generate(row.begin(), row.end(), g);
    return row;
}

matrix generate(size_t length)
{
    // generate rows
    std::vector<matrix::row> a(sigma);
    for (auto& row : a)
    {
        row = generate_row(length);
    }

    // normalize columns
    for (size_t j = 0; j < length; ++j)
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

void print_matrix(const matrix& matrix)
{
    for (const auto& row : matrix.get_data())
    {
        for (const auto& el : row)
        {
            std::cout << std::fixed << std::setprecision(8) << el << "\t";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}
