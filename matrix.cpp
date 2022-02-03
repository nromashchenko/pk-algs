#include "matrix.h"

#include <utility>
#include <algorithm>
#include <iostream>
#include <iomanip>
/*
std::vector<score_t> get_column(const matrix& matrix, size_t j)
{
    std::vector<score_t> column;
    for (size_t i = 0; i < sigma; ++i)
    {
        const auto& row = matrix[i];
        column.push_back(row[j]);
    }
    return column;
}*/

matrix::matrix()
    //: sorted(false)
{

}

matrix::matrix(std::vector<column> d)
    : data(std::move(d))
    //  , sorted(false)
{
}

score_t matrix::get(size_t i, size_t j) const
{
    return data[j][i];
}

size_t matrix::width() const
{
    return data.size();
}

bool matrix::empty() const
{
    return data.empty();
}

void matrix::sort()
{
    /*
    sorted = true;

    if (order.empty())
    {
        const size_t length = data.size();
        order = { sigma, std::vector<size_t>(length, sigma + 1)};
    }*/


    for (size_t j = 0; j < data.size(); ++j)
    {
        // data are stored in rows, extract the column
        //auto column = get_column(data, j);
        const auto& column = data[j];

        // sort it by score
        /*struct score_pair
        {
            score_t score;
            size_t order;
        };

        std::vector<score_pair> pairs(sigma);
        for (size_t i = 0; i < sigma; ++i)
        {
            pairs[i] = { get(i, j), i };
        }

        auto compare = [](const score_pair& p1, const score_pair& p2) { return p1.score > p2.score; };
        std::sort(begin(pairs), end(pairs), compare);


        // save the sorting
        for (size_t i = 0; i < sigma; ++i)
        {
            data[i][j] = pairs[i].score;
            order[i][j] = pairs[i].order;
        }
        */
    }
}
/*
bool matrix::is_sorted() const
{
    return sorted;
}*/

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

const std::vector<matrix::column>& matrix::get_data() const
{
    return data;
}

std::vector<matrix::column>& matrix::get_data()
{
    return data;
}

/*
matrix::column& matrix::get_row(size_t i)
{
    return data[i];
}

std::vector<std::vector<size_t>> matrix::get_order() const
{
    return order;
}*/


window::window(matrix& m, size_t start_pos, size_t size)
    : _matrix(m), _start_pos(start_pos), _size(size)
{

}

window& window::operator=(const window& other)
{
    if (*this != other)
    {
        _matrix = other._matrix;
        _start_pos = other._start_pos;
        _size = other._size;
    }
    return *this;
}

bool window::operator==(const window& other) const
{
    /// FIXME:
    /// Let us imagine those windows are over the same matrix
    return _start_pos == other._start_pos && _size == other._size;
}

bool window::operator!=(const window& other) const
{
    return !(*this == other);
}

score_t window::get(size_t i, size_t j) const
{
    return _matrix.get(i, _start_pos + j);
}

size_t window::size() const
{
    return _size;
}

bool window::empty() const
{
    return _size == 0;
}

size_t window::get_position() const
{
    return _start_pos;
}

std::pair<size_t, score_t> window::max_at(size_t column) const
{
    size_t max_index = 0;
    score_t max_score = get(0, column);
    for (size_t i = 1; i < sigma; ++i)
    {
        if (get(i, column) > max_score)
        {
            max_score = get(i, column);
            max_index = i;
        }
    }
    return { max_index, max_score };
}


impl::window_iterator::window_iterator(matrix& matrix, size_t kmer_size) noexcept
    : _matrix(matrix), _window(matrix, 0, kmer_size), _kmer_size(kmer_size), _current_pos(0)
{
}

impl::window_iterator& impl::window_iterator::operator++()
{
    _current_pos++;
    if (_current_pos + _kmer_size < _matrix.width())
    {
        _window = window(_matrix, _current_pos, _kmer_size);
    }
    else
    {
        // end iterator
        _window = window(_matrix, 0, 0);
        _kmer_size = 0;
    }
    return *this;
}

bool impl::window_iterator::operator==(const window_iterator& rhs) const noexcept
{
    return _window == rhs._window;
}

bool impl::window_iterator::operator!=(const window_iterator& rhs) const noexcept
{
    return !(*this == rhs);
}

impl::window_iterator::reference impl::window_iterator::operator*() noexcept
{
    return _window;
}

to_windows::to_windows(matrix& matrix, size_t kmer_size)
    : _matrix{ matrix }, _kmer_size{ kmer_size }//, _start_pos{ 0 }
{}

to_windows::const_iterator to_windows::begin() const
{
    return { _matrix, _kmer_size };
}

to_windows::const_iterator to_windows::end() const noexcept
{
    return { _matrix, 0 };
}


score_t g()
{
    return distr(eng);
}

matrix::column generate_column(size_t sigma)
{
    std::vector<score_t> column(sigma);
    std::generate(column.begin(), column.end(), g);
    return column;
}

matrix generate(size_t length)
{
    // generate columns
    std::vector<matrix::column> a(length);
    for (auto& column : a)
    {
        column = generate_column(sigma);
    }

    // normalize columns
    for (auto& column : a)
    {
        score_t sum = std::accumulate(column.begin(), column.end(), 0.0f);
        for (auto& elem : column)
        {
            elem /= sum;
        }
    }
    return a;
}

void print_matrix(const matrix& matrix)
{
    for (const auto& column : matrix.get_data())
    {
        for (const auto& el : column)
        {
            std::cout << std::fixed << std::setprecision(8) << el << "\t";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}