#include "dc.h"


std::vector<phylo_kmer> as_column(const window& window, size_t j, score_t eps)
{
    std::vector<phylo_kmer> column;
    for (size_t i = 0; i < sigma; ++i)
    {
        const auto& element = window.get(i, j);
        if (element > eps)
        {
            column.push_back({i, element });
        }
    }
    return column;
}



divide_and_conquer::divide_and_conquer(map_t& map, const window& window, size_t k)
        : _window(window)
        , _map(map)
        , _k(k)
{
    /// kmer_size can also be zero, which means the end() iterator
    const auto halfsize = size_t{ k / 2 };
    _prefix_size = (halfsize >= 1) ? halfsize : k;

    _result_list.reserve(std::pow(sigma, k));
}

void divide_and_conquer::run(score_t omega)
{
    // preprocessing O(k): range product query
    preprocess();

    const score_t eps = std::pow((omega / 4), _k);
    _result_list = dc(omega, 0, _k, eps);

    //const auto kmers = dc(omega, 0, _k, eps);
    /*
    for (const auto& [kmer, score] : kmers)
    {
        _map[kmer] = score;
    }*/
}

// j is the starat position of the window
// h is the length of the window
std::vector<phylo_kmer> divide_and_conquer::dc(score_t omega, size_t j, size_t h, score_t eps)
{
    // trivial case
    if (h == 1)
    {
        return as_column(_window, j, eps);
    }
    else
    {
        std::vector<phylo_kmer> result_vector;
        std::vector<phylo_kmer>& result = (h == _k) ? _result_list : result_vector;

        score_t eps_l = eps / best_score(j + h / 2, h - h / 2);
        score_t eps_r = eps / best_score(j, h / 2);

        auto l = dc(omega, j, h / 2, eps_l);
        auto r = dc(omega, j + h / 2, h - h / 2, eps_r);

        // let's sort not suffixes, but whichever is less to sort, suffixes or prefixes
        auto min = (l.size() < r.size()) ? l : r;
        auto max = (l.size() < r.size()) ? r : l;
        bool prefix_sort = l.size() < r.size();
        if (!min.empty())
        {
            std::sort(min.begin(), min.end(), kmer_score_comparator);

            for (const auto& [a, a_score] : max)
            {
                for (const auto& [b, b_score] : min)
                {
                    //const auto score = prefix_score + suffix_score;
                    const auto score = a_score * b_score;
                    if (score <= eps)
                    {
                        break;
                    }

                    code_t kmer;
                    if (prefix_sort)
                    {
                        kmer = (b << ((h - h / 2) * 2)) | a;
                    }
                    else
                    {
                        kmer = (a << ((h - h / 2) * 2)) | b;
                    }
                    result.push_back({ kmer, score });
                }
            }
        }
        return result;
    }
}

void divide_and_conquer::preprocess()
{
    _best_scores = std::vector<score_t>(_k, 0.0f);
    score_t product = 1.0f;
    for (size_t j = 0; j < _k; ++j)
    {
        const auto& [index_best, score_best] = _window.max_at(j);
        product *= score_best;
        _best_scores[j] = product;
    }

    /// We ignore 0s because they can not be highest values in the column
}

score_t divide_and_conquer::best_score(size_t start_pos, size_t h)
{
    /*score_t result = 1.0f;
    for (size_t j = start_pos; j < start_pos + h; ++j)
    {
        const auto& [index_best, score_best] = _window.max_at(j);
        result *= score_best;
    }
    return result;*/
    return _best_scores[start_pos + h - 1] / _best_scores[start_pos];
}

const map_t& divide_and_conquer::get_map()
{
    return _map;
}

const std::vector<phylo_kmer>& divide_and_conquer::get_result() const
{
    return _result_list;
}


size_t divide_and_conquer::get_num_kmers() const
{
    return _result_list.size();
}