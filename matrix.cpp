#include "matrix.h"
#include <range/v3/numeric/accumulate.hpp>
#include <iostream>

#include <iomanip>
#include <range/v3/algorithm/for_each.hpp>
#include "matrix.h"

auto print_elem = [](auto elem){
    std::cout << std::fixed << std::setprecision(8) << elem << "\t";
    //std::cout << elem << ' ';
};

auto print = [](auto rng) {
    rs::for_each(rng, print_elem);
    std::cout << std::endl;
};

auto print2D = [](auto rng) {
    for (auto r: rng)
    {
        print(r);
    }
};

matrix::matrix(std::vector<score_t> data, size_t num_cols)
    : data(std::move(data))
    , num_cols(num_cols)
{
    normalize();
}

void matrix::normalize()
{
    // normalizes a column
    auto normalize_column = [](auto column)
    {
        score_t sum = rs::accumulate(column, 0.0f);
        print(column);
        std::cout << "SUM: " << sum << std::endl;
        return column | vs::transform([sum](score_t v) { return v / sum; });
    };

    // transposes a matrix
    auto transpose = [](auto rng) {
        auto flat = rng | vs::join;
        int height = rs::distance(rng);
        int width = rs::distance(flat) / height;
        auto inner = [=](int i) {
            return flat | vs::drop(i) | vs::stride(width);
        };
        return vs::ints(0,width) | vs::transform(inner);
    };

    // normalizes a matrix
    auto normalized =
        get_columns(0, num_cols)
        | vs::transform(normalize_column)
        | vs::transform(transpose)
        | vs::join
        | rs::to_vector;


    std::vector<score_t> copy;
    std::copy(normalized.begin(), normalized.end(), std::back_inserter(copy));
    data = std::move(copy);

    //data = normalized;
}

/*
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


window::window(const matrix& m, size_t start_pos, size_t k)
    : data(m)
    , start_pos(start_pos)
    , k(k)
{
    if (start_pos < m.size())
    {
        throw std::runtime_error("Wrong start position of the window: " + std::to_string(start_pos));
    }

    if (start_pos + k > m.size())
    {
        throw std::runtime_error("The window of size " + std::to_string(k) + " does not fit");
    }
}


std::pair<size_t, score_t> window::max_at(size_t column) const
{
    size_t max_index = 0;

    size_t pos = start_pos;
    score_t max_score = data.get_data()[pos][column];
    for (size_t i = pos + 1; i < start_pos + k; ++i)
    {
        if (data.get_data()[i][column] > max_score)
        {
            max_score = data.get_data()[i][column];
            max_index = i;
        }
    }
    return { max_index, max_score };
}

score_t window::get_value(size_t row, size_t column) const
{
    return data.get_data()[row][start_pos + column];
}


size_t window::get_order(size_t row, size_t column) const
{
    return data.get_order()[row][start_pos + column];
}

std::vector<score_t> window::operator[](size_t row) const
{
    return data.get_data()[row];
}

size_t window::size() const
{
    return k;
}

bool window::empty() const
{
    return k == 0;
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
    return matrix(a);
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
*/


score_t g()
{
    return distr(eng);
}

matrix generate(size_t num_cols)
{
    // generate rows
    std::vector<score_t> a(sigma * num_cols);
    std::generate(a.begin(), a.end(), g);

    // normalize columns
    matrix m(a, num_cols);


    /*for (size_t j = 0; j < length; ++j)
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
    }*/

    return m;
}
