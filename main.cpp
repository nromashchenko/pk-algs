#include <iostream>
#include <unordered_map>
#include <vector>
#include <cmath>
#include <cassert>
#include <chrono>
#include <fstream>
#include <iterator>
#include <filesystem>

#include "common.h"
#include "dc.h"
#include "bb.h"
#include "brute_force.h"
#include "ar.h"

struct run_params
{
    size_t k;
    score_t omega;
};

struct flags
{
    bool run_bb;
    bool run_dc;
    bool run_dccw;
};

const std::vector<run_params> params =
    {
        { 6, 1.0},
        { 7, 1.0},
        { 8, 1.0},
        { 9, 1.0},
        { 10, 1.0},
        { 11, 1.0},
        { 12, 1.0},

        { 6, 1.25},
        { 7, 1.25},
        { 8, 1.25},
        { 9, 1.25},
        { 10, 1.25},
        { 11, 1.25},
        { 12, 1.25},
        //{ 13, 1.25},
        //{ 14, 1.25},

        { 6, 1.5},
        { 7, 1.5},
        { 8, 1.5},
        { 9, 1.5},
        { 10, 1.5},
        { 11, 1.5},
        { 12, 1.5},
        { 13, 1.5},
        { 14, 1.5},

        { 6, 1.75},
        { 7, 1.75},
        { 8, 1.75},
        { 9, 1.75},
        { 10, 1.75},
        { 11, 1.75},
        { 12, 1.75},
        { 13, 1.75},
        { 14, 1.75},

        { 6, 2.0},
        { 7, 2.0},
        { 8, 2.0},
        { 9, 2.0},
        { 10, 2.0},
        { 11, 2.0},
        { 12, 2.0},
        { 13, 2.0},
        { 14, 2.0},

    };

const size_t const_k = 10;

const std::vector<run_params> params_default =
    {
        { 10, 1.5 },
    };

const std::vector<run_params> params_k12 =
    {
        { 12, 1.5 },
    };


//const size_t const_k = 5;
const std::vector<run_params> params_one_k =
{
    { 6, 1.0},
    { 7, 1.0},
    { 8, 1.0},
    { 9, 1.0},
    { 10, 1.0},
    { 11, 1.0},
    { 12, 1.0},
};

const std::vector<run_params> params_omega_1 =
    {
        { 6, 1.0},
        //{ 7, 1.0},
        { 8, 1.0},
        //{ 9, 1.0},
        { 10, 1.0},
        //{ 11, 1.0},
        //{ 12, 1.0},
    };

const std::vector<run_params> params_omega_1_even_k =
    {
        { 6, 1.0},
        { 8, 1.0},
        { 10, 1.0},
        { 12, 1.0},
    };

const std::vector<run_params> params_omega_1_5 =
    {
        { 6, 1.5},
        { 7, 1.5},
        { 8, 1.5},
        { 9, 1.5},
        { 10, 1.5},
        { 11, 1.5},
        { 12, 1.5},
    };

const std::vector<run_params> params_omega_1_5_even_k =
    {
        { 6, 1.5},
        { 8, 1.5},
        { 10, 1.5},
        { 12, 1.5},
        //{ 14, 1.5},
    };

const std::vector<run_params> params_omega_0 =
    {
        { 6, 0},
        { 8, 0},
    };

const std::vector<run_params> params_test =
    {
        //{ 4, 1.0},
        { 6, 1.5},
        //{ 8, 1.0},
        //{ 4, 1.5},
        //{ 6, 1.5},
        //{ 8, 1.5},
    };

const std::vector<run_params> params_omega_2_even_k =
    {
        { 6, 2.0 },
        { 8, 2.0 },
        { 10, 2.0 },
        { 12, 2.0 },
        //{ 14, 1.5},
    };

void print_map(const map_t& map)
{
    for (const auto& [kmer, score] : map)
    {
        std::cout << kmer << ": " << score << std::endl;
    }
    std::cout << std::endl;
}

void assert_equal_map(const map_t& map1, const map_t& map2)
{
    if (map1.size() != map2.size())
    {
        std::cout << "Different sizes: " << map1.size() << ", " << map2.size() << std::endl;
        return;
    }

    // allows one mismatch
    //assert(fabs(map1.size() - map2.size()) < 2);
    assert(map1.size() == map2.size());

    for (const auto& [kmer, score] : map1)
    {

        assert(map2.find(kmer) != map2.end());

        assert(fabs(score - map2.at(kmer)) < 1e-6);
    }
}

void assert_equal(const std::vector<phylo_kmer>& a, const std::vector<phylo_kmer>& b)
{
    assert(a.size() == b.size());

    std::unordered_map<code_t, score_t> map_a;
    for (const auto& [kmer, score] : a)
    {
        map_a[kmer] = score;
    }

    std::unordered_map<code_t, score_t> map_b;
    for (const auto& [kmer, score] : b)
    {
        map_b[kmer] = score;
    }
    assert_equal_map(map_a, map_b);
}

void check_size(const std::vector<phylo_kmer>& a, const std::vector<phylo_kmer>& b)
{
    if (a.size() != b.size())
    {
        std::cout << "Different sizes: " << a.size() << ", " << b.size() << std::endl;
    }
}

score_t shannon(const std::vector<score_t>& values)
{
    score_t result = 0.0;
    for (const auto& v : values)
    {
        result += v * static_cast<score_t>(log2(v));
    }
    return - result;
}

// Get the order of columns by entropy
std::vector<column_data> get_order(const window& window)
{
    std::vector<column_data> order;
    const auto k = window.size();
    for (size_t j = 0; j < k; ++j)
    {
        order.push_back({j,  shannon(window.get_column(j))});
    }

    auto entropy_compare = [](const auto& a, const auto& b) {
        return a.entropy < b.entropy;
    };
    std::sort(order.begin(), order.end(), entropy_compare);
    return order;
}

void test_one(size_t k, bool print=true)
{
    const score_t omega = 1.0;

    auto matrix = generate(2 * k);
    if (print)
    {
        print_matrix(matrix);
        std::cout << "Threshold: " << std::pow((omega / 4), k) << std::endl;
    }

    for (const auto& window : to_windows(matrix, k))
    {
        //std::cout << "WINDOW: " << window.get_position() << std::endl;
        branch_and_bound bb(window, k, omega);
        bb.run(omega);
        if (print)
        {
            //print_map(bb.get_map());
            std::cout << "Branch-and-bound, generated: " << bb.get_result().size() << std::endl;
        }


        /*
        std::cout << "Order:" << std::endl;
        for (const auto &[j, e] : order)
        {
            std::cout << j << ":" << e << std::endl;
        }*/
/*
        bbe bbe(window, std::move(get_order(window)), k);
        bbe.run(omega);
        if (print)
        {
            //print_map(bb.get_map());
            std::cout << "Branch-and-bound Sorted, generated: " << bbe.get_result().size() << std::endl;
        }

        assert_equal(bb.get_result(), bbe.get_result());*/

        divide_and_conquer dc(window, k, omega);
        dc.run(omega);
        if (print)
        {
            //print_map(dc.get_map());
            std::cout << "Divide-and-conquer, generated: " << dc.get_result().size() << std::endl;
            std::cout << std::endl;
        }

/*
        brute_force bf(window, k);
        bf.run(omega);
        if (print)
        {
            //print_map(bf.get_map());
            std::cout << "Brute force, generated: " << bf.get_map().size() << std::endl;
        }

        matrix.sort();*/


        assert_equal(bb.get_result(), dc.get_result());
        //assert_equal(dc.get_result(), dccw.get_result());

        //assert_equal(bb.get_map(), bf.get_map());
        //assert_equal(rap.get_map(), bf.get_map());
    }

}

void test_suite()
{
    const size_t num_iter = 100;
    const std::vector<size_t> k_values = { 6, 7, 8, 9 };

    for (const auto k : k_values)
    {
        for (size_t i = 0; i < num_iter; ++i)
        {
            std::cout << "\rTesting k = " << k << ". " << i << " / " << num_iter << "..." << std::flush;
            test_one(k, false);
        }
        std::cout << "\rTesting k = " << k << ". Done." << std::endl;
    }
}


enum class algorithm
{
    bb = 0,
    dc = 1,
    rappas = 2,
    dccw = 3,
    baseline = 4,
    bbe = 5
};

struct run_stats
{
    algorithm alg;
    size_t num_kmers;
    unsigned long time;
    size_t k;
    score_t omega;
    std::string node;
    size_t window_pos;
};

void print_as_csv(const std::vector<run_stats>& stats, const std::string& filename)
{
    std::cout << "Writing results: " << filename << "...";

    std::ofstream file(filename);
    file << "alg,num_kmers,time,k,omega,node,position" << std::endl;
    for (const auto& stat: stats)
    {
        const auto& [alg, num_kmers, time, k, omega, node, position] = stat;

        switch (alg)
        {
            case algorithm::bb:
                file << "bb";
                break;
            case algorithm::dc:
                file << "dc";
                break;
            case algorithm::dccw:
                file << "dccw";
                break;
            case algorithm::rappas:
                file << "rappas";
                break;
            case algorithm::baseline:
                file << "bl";
                break;
            case algorithm::bbe:
                file << "bbe";
                break;
        }
        file << "," << num_kmers << "," << time << "," << k << "," << omega << "," << node << "," << position << std::endl;
    }
    file.close();

    std::cout << std::endl;
}

struct bb_stats
{
    size_t k;
    score_t omega;
    std::vector<bb_return> returns;
};


void print_returns(const std::vector<bb_stats>& stats, const std::string& filename)
{
    std::cout << "Writing results: " << filename << "...";

    std::ofstream file(filename);
    for (const auto& [k, omega, stat]: stats)
    {
        file << k << "," << omega << ",";
        for (const auto& ret : stat)
        {
            switch (ret)
            {
                case bb_return::BAD_PREFIX:
                    file << 0;
                    break;
                case bb_return::GOOD_KMER:
                    file << 1;
                    break;
                case bb_return::GOOD_PRFIX:
                    file << 2;
                    break;
            }
            //file << ",";
        }
        file << std::endl;
    }
    file.close();

    std::cout << std::endl;
}

std::tuple<std::vector<phylo_kmer>, run_stats> run_bb(const window& window, size_t k, float omega,
                                                            const std::string& node_name)
{
    branch_and_bound bb(window, k, omega);
    auto begin = std::chrono::steady_clock::now();
    bb.run(omega);
    auto end = std::chrono::steady_clock::now();
    unsigned long time = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();
    const auto stats = run_stats{
                             algorithm::bb,
                             bb.get_num_kmers(),
                             time,
                             k, omega,
                             node_name,
                             window.get_position()
                         };
    return { bb.get_result(), stats };
}

std::tuple<std::vector<phylo_kmer>, run_stats> run_dc(const window& window, size_t k, float omega,
                                                      const std::string& node_name)
{
    divide_and_conquer dc(window, k, omega);
    auto begin = std::chrono::steady_clock::now();
    dc.run(omega);
    auto end = std::chrono::steady_clock::now();
    unsigned long time = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();
    const auto stats = run_stats{
        algorithm::dc,
        dc.get_num_kmers(),
        time,
        k, omega,
        node_name,
        window.get_position()
    };
    return { dc.get_result(), stats };
}

std::tuple<std::vector<phylo_kmer>, run_stats> run_dccw(std::vector<phylo_kmer>& prefixes,
                                                        const window& prev, const window& current, const window& next,
                                                        size_t k, float omega,
                                                        const std::string& node_name)
{
    score_t lookbehind = get_threshold(omega, k);
    if (prev.get_position() < current.get_position())
    {
        lookbehind = prev.range_product(0, k / 2);;
    }
    else
    {
        prefixes.clear();
    }
    score_t lookahead = get_threshold(omega, k);
    if (next.get_position() > current.get_position())
    {
        lookahead = next.range_product(k / 2, k - k / 2);
    }

    dccw dccw(current, prefixes, k, lookbehind, lookahead, omega);
    auto begin = std::chrono::steady_clock::now();
    dccw.run(omega);
    auto end = std::chrono::steady_clock::now();
    unsigned long time = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();
    prefixes = std::move(dccw.get_suffixes());

    const auto stats = run_stats{
        algorithm::dccw,
        dccw.get_num_kmers(),
        time,
        k, omega,
        node_name,
        current.get_position()
    };
    return { dccw.get_result(), stats };
}

void test_random(const flags& flags,
                 const std::vector<run_params>& parameters,
                 size_t num_iter, const std::string& filename)
{
    const auto node_name = "Random" + std::to_string(num_iter);

    std::vector<run_stats> stats;

    /// for debugging
    std::vector<phylo_kmer> bb_result;
    std::vector<phylo_kmer> dc_result;
    std::vector<phylo_kmer> dccw_result;
    (void)bb_result; (void)dc_result; (void)dccw_result;


    for (const auto& [k, omega] : parameters)
    {
        for (size_t i = 0; i < num_iter; ++i)
        {
            std::cout << "\r\tRunning for k = " << k << ", omega = " << omega << ". " << i << " / " << num_iter << "..." << std::flush;

            auto matrix = generate(1000);

            std::vector<phylo_kmer> prefixes;

            for (const auto& [prev, window, next] : chain_windows(matrix, k))
            {
                if (flags.run_bb)
                {
                    const auto& [result, stat] = run_bb(window, k, omega, node_name);
                    stats.push_back(stat);
                    bb_result = result;
                }

                if (flags.run_dc)
                {
                    const auto& [result, stat] = run_dc(window, k, omega, node_name);
                    stats.push_back(stat);
                    dc_result = result;
                }

                if (flags.run_dccw)
                {
                    const auto& [result, stat] = run_dccw(prefixes, prev, window, next, k, omega, node_name);
                    stats.push_back(stat);
                    dccw_result = result;
                }

                //assert_equal(bb_result, dc_result);
                //assert_equal(dc_result, dccw_result);
            }
        }
        std::cout << "\r\tRunning for k = " << k << ", omega = " << omega << ". Done." << std::endl;
    }

    print_as_csv(stats, filename);
    //print_returns(bb_stats, "returns.txt");
}


std::vector<std::string> get_ghost_ids(const std::string& filename)
{
    std::vector<std::string> result;
    std::ifstream file(filename);

    std::copy(std::istream_iterator<std::string>(file),
              std::istream_iterator<std::string>(),
              back_inserter(result));
    return result;
}



void test_data(const flags& flags,
               const std::vector<run_params>& parameters, const std::string& input,
               const std::string& ghost_ids_file, const std::string& output)
{
    const auto ghost_ids = get_ghost_ids(ghost_ids_file);

    raxmlng_reader reader(input);
    auto matrices = reader.read();

    std::unordered_map<std::string, matrix> sample;
    for (const auto& [k, v] : matrices)
    {
        if (const auto& it = std::find(ghost_ids.begin(), ghost_ids.end(), k); it != ghost_ids.end())
        {
            sample[k] = v;
        }
    }

    /// for debugging
    std::vector<phylo_kmer> bb_result;
    std::vector<phylo_kmer> dc_result;
    std::vector<phylo_kmer> dccw_result;
    (void)bb_result; (void)dc_result; (void)dccw_result;

    std::cout << "Num matrices: " << sample.size() << std::endl;

    std::vector<run_stats> stats;
    size_t node_i = 0;
    for (auto& [node_name, matrix] : sample)
    {
        if (node_i % 1 == 0)
        {
            std::cout << "\r\tRunning for node " << node_name << ", " << node_i << " / " << sample.size() << "..." << std::flush;
        }

        for (const auto& [k, omega] : parameters)
        {

            std::vector<phylo_kmer> prefixes;
            for (const auto& [prev, window, next] : chain_windows(matrix, k))
            //for (const auto& window : to_windows(matrix, k))
            {
                if (flags.run_bb)
                {
                    const auto& [result, stat] = run_bb(window, k, omega, node_name);
                    stats.push_back(stat);
                    bb_result = result;
                }

                if (flags.run_dc)
                {
                    const auto& [result, stat] = run_dc(window, k, omega, node_name);
                    stats.push_back(stat);
                    dc_result = result;
                }

                if (flags.run_dccw)
                {
                    const auto& [result, stat] = run_dccw(prefixes, prev, window, next, k, omega, node_name);
                    stats.push_back(stat);
                    dccw_result = result;
                }

                //assert_equal(bb_result, dc_result);
                //assert_equal(dc_result, dccw_result);
            }
        }

        if (node_i % 1 == 0)
        {
            std::cout << "\r\tRunning for node " << node_name << ", " << node_i << " / " << sample.size() << ". Done.\n"
                      << std::flush;
        }
        node_i++;
    }

    print_as_csv(stats, output);
}

int main(int argc, char** argv)
{
    srand(42);

    //test_one(10);
    //test_suite();

    //const auto parameters = params;
    //const auto parameters = params_test;
    const auto parameters = params_default;
    //const auto parameters = params_k12;

    //const auto parameters = params_one_k;
    //const auto parameters = params_omega_1;
    //const auto parameters = params_omega_1_even_k;
    //const auto parameters = params_omega_1_5;
    //const auto parameters = params_omega_1_5_even_k;
    //const auto parameters = params_omega_0;
    //const auto parameters = params_omega_2_even_k;

    flags alg_flags = { true, true, true };

    if (argc > 2)
    {
        if (argc != 7)
        {
            std::cout << "Usage:\n\t"
                << argv[0] << "\n\n or \n\n\t"
                << argv[0] << " <RAxML-NG output file> <Ghost ID file> 0/1[run BB] 0/1[run DC] 0/1[run DCCW] OUTPUT_FILE" << std::endl;
            return 1;
        }
        std::string filename = argv[1];
        std::string ghost_ids_file = argv[2];
        bool run_bb = static_cast<bool>(std::stoi(argv[3]));
        bool run_dc = static_cast<bool>(std::stoi(argv[4]));
        bool run_dccw = static_cast<bool>(std::stoi(argv[5]));
        std::string output_file = argv[6];

        if (std::filesystem::exists(output_file))
        {
            std::cerr << "File exists: " << output_file << std::endl;
            return 1;
        }

        test_data({ run_bb, run_dc, run_dccw }, parameters, filename, ghost_ids_file, output_file);
    }
    else
    {
        test_random(alg_flags, parameters, 1000, std::string(std::tmpnam(nullptr)) + ".csv");
    }

    return 0;
}
