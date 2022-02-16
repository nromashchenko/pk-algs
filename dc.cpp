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



divide_and_conquer::divide_and_conquer(const window& window, size_t k)
        : _window(window)
        , _k(k)
{
    /// kmer_size can also be zero, which means the end() iterator
    const auto halfsize = size_t{ k / 2 };
    _prefix_size = (halfsize >= 1) ? halfsize : k;

    _result_list.reserve(std::pow(sigma, k));

    // preprocessing O(k): range product query
    preprocess();
}

void divide_and_conquer::run(score_t omega)
{
    const auto eps = get_threshold(omega, _k);

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
        bool prefix_sort = l.size() < r.size();
        auto min = prefix_sort ? l : r;
        auto max = prefix_sort ? r : l;

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
/*
const map_t& divide_and_conquer::get_map()
{
    return _map;
}*/

const std::vector<phylo_kmer>& divide_and_conquer::get_result() const
{
    return _result_list;
}


size_t divide_and_conquer::get_num_kmers() const
{
    return _result_list.size();
}

dccw::dccw(const window& window, std::vector<phylo_kmer>& prefixes, size_t k, score_t lookahead, score_t lookbehind)
    : _window(window)
    , _prefixes(prefixes)
    , _k(k)
    , _lookahead(lookahead)
    , _lookbehind(lookbehind)
    , _dc(window, k)
{
    const auto halfsize = size_t{ k / 2 };
    _prefix_size = (halfsize >= 1) ? halfsize : k;

    _result_list.reserve(std::pow(sigma, k));

    _dc.preprocess();
}

void dccw::run(score_t omega)
{
    const score_t eps = std::pow((omega / 4), _k);
    score_t eps_r = eps / _dc.best_score(0, _k / 2);
    score_t eps_l = eps / _dc.best_score(_k / 2, _k - _k / 2);

    auto& L = _prefixes;
    if (L.empty())
    {
        L = _dc.dc(omega, 0, _k / 2, eps_l);
    }

    _suffixes = std::move(_dc.dc(omega, _k / 2, _k - _k / 2, std::max(eps_r, _lookahead)));
    auto& R = _suffixes;

    // Let's sort not suffixes, but whichever is less to sort, suffixes or prefixes
    // The trick is that L can contain more dead prefixes (which were alive suffixes in the previous window),
    // and R can contain dead suffixes (which will be alive prefixes in the next window).
    // Let's first find how many alive prefixes and suffixes we have in the current window.
    auto last_prefix = L.end();
    auto last_suffix = R.end();

    // The best prefix score of the previous window was better than the best suffix score of
    // the current window. Then, more strings of L are alive in W_prev than in W.
    // => Need to partition L to find only the part of strings that are alive in W.
    if (_lookbehind > eps_l)
    {
        auto is_dead_prefix = [eps_l](const auto& pk) { return pk.score < eps_l; };
        last_prefix = std::partition(L.begin(), L.end(), is_dead_prefix);
    }

    // The same for strings of R and the lookahead score for the next window
    if (_lookahead > eps_r)
    {
        auto is_dead_suffix = [eps_r](const auto& pk) { return pk.score < eps_r; };
        last_suffix = std::partition(R.begin(), R.end(), is_dead_suffix);
    }

    size_t num_alive_prefixes = std::distance(L.begin(), last_prefix);
    size_t num_alive_suffixes = std::distance(R.begin(), last_suffix);

    bool prefix_sort = num_alive_prefixes < num_alive_suffixes;
    auto min = prefix_sort ? L : R;
    auto max = prefix_sort ? R : L;
    auto last_min = prefix_sort ? last_prefix : last_suffix;
    //auto last_max = prefix_sort ? last_suffix : last_prefix;

    if (!min.empty())
    {
        auto eps_min = prefix_sort ? eps_l : eps_r;
        auto eps_max = prefix_sort ? eps_r : eps_l;

        std::sort(min.begin(), last_min, kmer_score_comparator);

        //for (const auto& [a, a_score] : max)
        size_t i = 0;
        while (i < max.size())
        {
            const auto& [a, a_score] = max[i];
            if (a_score < eps_max)
            {
                break;
            }

            //    for (const auto& [b, b_score] : min)
            size_t j = 0;
            while (j < min.size())
            {
                const auto& [b, b_score] = min[j];
                if (b_score < eps_min)
                {
                    break;
                }

                //const auto score = prefix_score + suffix_score;
                const auto score = a_score * b_score;
                if (score <= eps)
                {
                    break;
                }

                code_t kmer;
                if (prefix_sort)
                {
                    kmer = (b << ((_k - _k / 2) * 2)) | a;
                }
                else
                {
                    kmer = (a << ((_k - _k / 2) * 2)) | b;
                }
                _result_list.push_back({ kmer, score });

                j++;
            }

            i++;
        }
    }
}


const std::vector<phylo_kmer>& dccw::get_result() const
{
    return _result_list;
}


size_t dccw::get_num_kmers() const
{
    return _result_list.size();
}

std::vector<phylo_kmer>&& dccw::get_suffixes()
{
    return std::move(_suffixes);
}