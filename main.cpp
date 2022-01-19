#include <iostream>
#include <unordered_map>
#include <vector>
#include <cmath>
#include <cassert>
#include <chrono>
#include <fstream>

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

const std::vector<run_params> params =
    {
        { 6, 1.0},
        { 7, 1.0},
        { 8, 1.0},
        { 9, 1.0},
        { 10, 1.0},
        //{ 11, 1.0},

        { 6, 1.25},
        { 7, 1.25},
        { 8, 1.25},
        { 9, 1.25},
        { 10, 1.25},
        { 11, 1.25},
        { 12, 1.25},
        /*{ 13, 1.25},
        { 14, 1.25},*/

        { 6, 1.5},
        { 7, 1.5},
        { 8, 1.5},
        { 9, 1.5},
        { 10, 1.5},
        { 11, 1.5},
        { 12, 1.5},
        //{ 13, 1.5},
        //{ 14, 1.5},

        { 6, 1.75},
        { 7, 1.75},
        { 8, 1.75},
        { 9, 1.75},
        { 10, 1.75},
        { 11, 1.75},
        { 12, 1.75},
        //{ 13, 1.75},
        //{ 14, 1.75},

        { 6, 2.0},
        { 7, 2.0},
        { 8, 2.0},
        { 9, 2.0},
        { 10, 2.0},
        { 11, 2.0},
        { 12, 2.0},
        /*{ 13, 2.0},
        { 14, 2.0},*/
    };

void print_map(const map_t& map)
{
    for (const auto& [kmer, score] : map)
    {
        std::cout << kmer << ": " << score << std::endl;
    }
    std::cout << std::endl;
}

void assert_equal(const map_t& map1, const map_t& map2)
{
    /*if (map1.size() != map2.size())
    {
        std::cout << "Different sizes: " << map1.size() << ", " << map2.size() << std::endl;
    }*/

    // allows one mismatch
    //assert(fabs(map1.size() - map2.size()) < 2);
    assert(map1.size() == map2.size());

    for (const auto& [kmer, score] : map1)
    {

        assert(map2.find(kmer) != map2.end());

        assert(fabs(score - map2.at(kmer)) < 1e-6);
    }
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
        branch_and_bound bb(window, k);
        bb.run(omega);
        if (print)
        {
            //print_map(bb.get_map());
            std::cout << "Branch-and-bound, generated: " << bb.get_map().size() << std::endl;
        }

        divide_and_conquer dc(window, k);
        dc.run(omega);
        if (print)
        {
            //print_map(dc.get_map());
            std::cout << "Divide-and-conquer, generated: " << dc.get_map().size() << std::endl;
        }

        brute_force bf(window, k);
        bf.run(omega);
        if (print)
        {
            //print_map(bf.get_map());
            std::cout << "Brute force, generated: " << bf.get_map().size() << std::endl;
        }

        matrix.sort();
        /*rappas rap(matrix, k);
        rap.run(omega);
        if (print)
        {
            //print_map(bb.get_map());
            std::cout << "Rappas, generated: " << rap.get_map().size() << std::endl;
        }*/

        assert_equal(bb.get_map(), bf.get_map());
        assert_equal(dc.get_map(), bf.get_map());
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
    rappas = 2
};

struct run_stats
{
    algorithm alg;
    size_t num_kmers;
    long time;
    size_t k;
    score_t omega;
};

void print_as_csv(const std::vector<run_stats>& stats, const std::string& filename)
{
    std::cout << "Writing results: " << filename << "...";

    std::ofstream file(filename);
    file << "alg,num_kmers,time,k,omega" << std::endl;
    for (const auto& stat: stats)
    {
        const auto& [alg, num_kmers, time, k, omega] = stat;

        switch (alg)
        {
            case algorithm::bb:
                file << "bb";
                break;
            case algorithm::dc:
                file << "dc";
                break;
            case algorithm::rappas:
                file << "rappas";
                break;
        }
        file << "," << num_kmers << "," << time << "," << k << "," << omega << std::endl;
    }
    file.close();

    std::cout << std::endl;
}

void test_random(size_t num_iter, const std::string& filename)
{
    std::vector<run_stats> stats;

    for (const auto& [k, omega] : params)
    {
        for (size_t i = 0; i < num_iter; ++i)
        {
            std::cout << "\r\tRunning for k = " << k << ", omega = " << omega << ". " << i << " / " << num_iter << "..." << std::flush;

            auto matrix = generate(5 * k);


            for (const auto& window : to_windows(matrix, k))
            {
                auto begin = std::chrono::steady_clock::now();
                branch_and_bound bb(window, k);
                bb.run(omega);
                auto end = std::chrono::steady_clock::now();
                long bb_time = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();
                stats.push_back({
                                    algorithm::bb,
                                    bb.get_map().size(), bb_time,
                                    k, omega
                                });

                begin = std::chrono::steady_clock::now();
                divide_and_conquer dc(window, k);
                dc.run(omega);
                end = std::chrono::steady_clock::now();
                long dc_time = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();
                stats.push_back({
                                    algorithm::dc,
                                    dc.get_map().size(), dc_time,
                                    k, omega
                                });

                /*
                matrix.sort();

                begin = std::chrono::steady_clock::now();
                rappas rap(matrix, k);
                rap.run(omega);
                end = std::chrono::steady_clock::now();
                long rap_time = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();
                stats.push_back({
                                        algorithm::rappas,
                                        rap.get_map().size(), rap_time,
                                        k, omega
                                });
                */

                assert_equal(bb.get_map(), dc.get_map());
                //assert_equal(bb.get_map(), rap.get_map());
            }
        }
        std::cout << "\r\tRunning for k = " << k << ", omega = " << omega << ". Done." << std::endl;
    }

    print_as_csv(stats, filename);
}

void test_data(const std::string& input, const std::string& output)
{
    raxmlng_reader reader(input);
    auto matrices = reader.read();

    const size_t sample_size = 10;
    std::unordered_map<std::string, matrix> sample;
    std::sample(matrices.begin(), matrices.end(), std::inserter(sample, sample.begin()),
                sample_size, std::mt19937{std::random_device{}()});

    std::vector<run_stats> stats;
    size_t node_i = 0;
    for (auto& [node, matrix] : sample)
    {
        if (node_i % 1 == 0)
        {
            std::cout << "\r\tRunning for node " << node << ", " << node_i << " / " << matrices.size() << "..." << std::flush;
        }

        for (const auto& [k, omega] : params)
        {
            for (const auto& window : to_windows(matrix, k))
            {
                auto begin = std::chrono::steady_clock::now();
                branch_and_bound bb(window, k);
                bb.run(omega);
                auto end = std::chrono::steady_clock::now();
                long bb_time = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();
                stats.push_back({
                                    algorithm::bb,
                                    bb.get_map().size(), bb_time,
                                    k, omega
                                });

                begin = std::chrono::steady_clock::now();
                divide_and_conquer dc(window, k);
                dc.run(omega);
                end = std::chrono::steady_clock::now();
                long dc_time = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();
                stats.push_back({
                                    algorithm::dc,
                                    dc.get_map().size(), dc_time,
                                    k, omega
                                });

                /*
                matrix.sort();

                begin = std::chrono::steady_clock::now();
                rappas rap(matrix, k);
                rap.run(omega);
                end = std::chrono::steady_clock::now();
                long rap_time = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();
                stats.push_back({
                                        algorithm::rappas,
                                        rap.get_map().size(), rap_time,
                                        k, omega
                                });
                */

                assert_equal(bb.get_map(), dc.get_map());
                //assert_equal(bb.get_map(), rap.get_map());
            }
        }

        if (node_i % 1 == 0)
        {
            std::cout << "\r\tRunning for node " << node << ", " << node_i << " / " << matrices.size() << ". Done.\n"
                      << std::flush;
        }
        node_i++;
    }

    print_as_csv(stats, output);
}


int main(int argc, char** argv)
{
    srand(42);

    //test_one(5);
    //test_suite();

    test_random(100, std::string(std::tmpnam(nullptr)) + ".csv");

    /*
    if (argc > 1)
    {
        std::string filename = argv[1];
        test_data(filename, std::string(std::tmpnam(nullptr)) + ".csv");
    }
    else
    {
        std::cout << "Usage:\n\t" << argv[0] << " FILENAME" << std::endl;
        std::cout << "The filename should be the AR result of RAxML-ng." << std::endl;
    }*/

    return 0;
}
