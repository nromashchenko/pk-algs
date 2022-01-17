#ifndef XPAS_ALGS_MATRIX_H
#define XPAS_ALGS_MATRIX_H

#include "common.h"
#include <random>
#include <range/v3/view.hpp>
#include <type_traits>


namespace rs = ranges;
namespace vs = ranges::views;


class matrix
{
private:
    std::vector<score_t> data;
    size_t num_cols;

    auto to_row()
    {
        return [this](size_t i){
            return data | vs::drop(i * num_cols) | vs::take(num_cols);
        };
    };

    auto to_column()
    {
        return [this](size_t i){
            return data | vs::drop(i) | vs::stride(num_cols); };
    };


public:
    matrix(std::vector<score_t> data, size_t num_cols);

    void normalize();

    auto get_rows(size_t start, size_t end)
    {
        return vs::ints(start, end) | vs::transform(to_row());
    }

    auto get_columns(size_t start, size_t end)
    {
        return vs::ints(start, end) | vs::transform(to_column());
    }
};


// window is a range of columns
using window = std::invoke_result<decltype(&matrix::get_columns)>;

/*
class matrix {
public:
    using row = std::vector<score_t>;

    explicit matrix(std::vector<row> d);

    //std::vector<score_t> operator[](size_t j) const;

    [[nodiscard]]
    size_t size() const;

    [[nodiscard]]
    bool empty() const;

    void sort();

    [[nodiscard]]
    bool is_sorted() const;


    [[nodiscard]]
    std::vector<row> get_data() const;

    [[nodiscard]]
    std::vector<std::vector<size_t>> get_order() const;

private:
    std::vector<row> data;

    bool sorted;
    std::vector<std::vector<size_t>> order;
};

class window
{
public:
    window(const matrix& m, size_t start_pos, size_t k);

    [[nodiscard]]
    std::pair<size_t, score_t> max_at(size_t column) const;

    [[nodiscard]]
    score_t get_value(size_t row, size_t column) const;

    [[nodiscard]]
    size_t get_order(size_t row, size_t column) const;

    [[nodiscard]]
    std::vector<score_t> operator[](size_t row) const;

    [[nodiscard]]
    size_t size() const;

    [[nodiscard]]
    bool empty() const;

private:
    const matrix& data;
    size_t start_pos;
    size_t k;
};*/

matrix generate(size_t length);

static std::random_device rd;
static std::default_random_engine eng(42);
static std::uniform_real_distribution<score_t> distr(0, 1);

#endif //XPAS_ALGS_MATRIX_H
