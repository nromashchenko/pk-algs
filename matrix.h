#ifndef XPAS_ALGS_MATRIX_H
#define XPAS_ALGS_MATRIX_H

#include "common.h"
#include <random>

class matrix {
public:
    using row = std::vector<score_t>;

    matrix(std::vector<row>  d);

    std::vector<score_t> operator[](size_t j) const;

    size_t size() const;

    bool empty() const;

    [[nodiscard]]
    std::pair<size_t, score_t> max_at(size_t column) const;

    std::vector<row> get_data() const;
private:
    std::vector<row> data;
};

matrix generate(size_t length);
void print_matrix(const matrix& matrix);

static std::random_device rd;
static std::default_random_engine eng(42);
static std::uniform_real_distribution<score_t> distr(0, 1);

#endif //XPAS_ALGS_MATRIX_H
