#include "matrix.h"

#include <utility>
#include <random>
#include <algorithm>
#include <iostream>
#include <iomanip>

matrix::matrix(std::vector<row>  d)
    : data(std::move(d))
{}

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
