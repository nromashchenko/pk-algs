#ifndef XPAS_ALGS_MATRIX_H
#define XPAS_ALGS_MATRIX_H

#include "common.h"
#include <random>

class matrix {
public:
    using row = std::vector<score_t>;

    matrix(std::vector<row> d);

    //std::vector<score_t> operator[](size_t j) const;

    score_t get(size_t i, size_t j) const;

    size_t width() const;

    bool empty() const;

    void sort();

    bool is_sorted() const;

    [[nodiscard]]
    std::pair<size_t, score_t> max_at(size_t column) const;

    std::vector<row> get_data() const;

    std::vector<std::vector<size_t>> get_order() const;

private:
    std::vector<row> data;

    bool sorted;
    std::vector<std::vector<size_t>> order;
};

class window
{
public:
    window(matrix& m, size_t start_pos, size_t size);

    window(window&&) noexcept = default;
    window& operator=(const window& other);
    window& operator=(window&&) noexcept = default;


    bool operator==(const window& other) const;
    bool operator!=(const window& other) const;

    score_t get(size_t i, size_t j) const;

    size_t size() const;

    bool empty() const;

    size_t get_position() const;

    void sort();

    bool is_sorted() const;

    [[nodiscard]]
    std::pair<size_t, score_t> max_at(size_t column) const;

    std::vector<matrix::row> get_data() const;

    std::vector<std::vector<size_t>> get_order() const;

private:
    matrix& _matrix;
    size_t _start_pos;
    size_t _size;
};

namespace impl
{
    class window_iterator
    {
    public:
        using iterator_category = std::forward_iterator_tag;
        using reference = window&;

        window_iterator(matrix& matrix, size_t kmer_size) noexcept;
        window_iterator(const window_iterator&) = delete;
        window_iterator(window_iterator&&) = delete;
        window_iterator& operator=(const window_iterator&) = delete;
        window_iterator& operator=(window_iterator&&) = delete;
        ~window_iterator() = default;

        window_iterator& operator++();

        bool operator==(const window_iterator& rhs) const noexcept;
        bool operator!=(const window_iterator& rhs) const noexcept;

        reference operator*() noexcept;
    private:
        matrix& _matrix;

        window _window;

        size_t _kmer_size;

        size_t _current_pos;
    };
}

class to_windows
{
public:
    using iterator_category = std::forward_iterator_tag;
    using const_iterator = impl::window_iterator;

    using reference = window&;

    to_windows(matrix& matrix, size_t kmer_size);
    to_windows(const to_windows&) = delete;
    to_windows(to_windows&&) = delete;
    to_windows& operator=(const to_windows&) = delete;
    to_windows& operator=(to_windows&&) = delete;
    ~to_windows() noexcept = default;

    [[nodiscard]]
    const_iterator begin() const;

    [[nodiscard]]
    const_iterator end() const noexcept;

private:
    matrix& _matrix;
    size_t _kmer_size;
};


matrix generate(size_t length);
void print_matrix(const matrix& matrix);

static std::random_device rd;
static std::default_random_engine eng(42);
static std::uniform_real_distribution<score_t> distr(0, 1);

#endif //XPAS_ALGS_MATRIX_H