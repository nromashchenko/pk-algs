#include "common.h"

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