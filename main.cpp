#include <iostream>
#include <unordered_map>
#include <vector>
#include <cmath>
#include <algorithm>
#include <random>
#include <iomanip>
#include <cassert>
#include <chrono>
#include <fstream>

#include "common.h"
#include "dc.h"
#include "bb.h"
#include "brute_force.h"


std::random_device rd;
std::default_random_engine eng(42);
std::uniform_real_distribution<score_t> distr(0, 1);

score_t g()
{
    return distr(eng);
}

row_t generate_row(size_t k)
{
    std::vector<score_t> row(k);
    std::generate(row.begin(), row.end(), g);
    return row;
}

matrix_t generate(size_t k)
{
    matrix_t a(sigma);
    for (size_t i = 0; i < a.size(); ++i)
    {
        a[i] = generate_row(k);
    }

    for (size_t j = 0; j < k; ++j)
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
    return a;
}

void print_matrix(const matrix_t& matrix)
{
    for (const auto& row : matrix)
    {
        for (const auto& el : row)
        {
            std::cout << std::fixed << std::setprecision(8) << el << "\t";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

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

    assert(map1.size() == map2.size());
    for (const auto& [kmer, score] : map1)
    {

        assert(map2.find(kmer) != map2.end());

        assert(fabs(score - map2.at(kmer)) < 1e-6);
    }
}

void test_one(size_t k, bool print=true)
{
    const auto omega = 1.0;

    const auto matrix = generate(k);
    if (print)
    {
        print_matrix(matrix);
        std::cout << "Threshold: " << std::pow((omega / 4), k) << std::endl;
    }

    branch_and_bound bb(matrix, k);
    bb.run();
    if (print)
    {
        //print_map(bb.get_map());
        std::cout << "Branch-and-bound, generated: " << bb.get_map().size() << std::endl;
    }

    divide_and_conquer dc(matrix, k);
    dc.run();
    if (print)
    {
        //print_map(dc.get_map());
        std::cout << "Divide-and-conquer, generated: " << dc.get_map().size() << std::endl;
    }

    brute_force bf(matrix, k);
    bf.run();
    if (print)
    {
        //print_map(bf.get_map());
        std::cout << "Brute force, generated: " << bf.get_map().size() << std::endl;
    }

    assert_equal(bb.get_map(), bf.get_map());
    assert_equal(dc.get_map(), bf.get_map());
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

void benchmark(size_t num_iter, const std::string& filename)
{
    //const std::vector<size_t> k_values = { 6, 7, 8, 9, 10, 11 };
    const std::vector<size_t> k_values = { 6, 7, 8, 9, 10, 11 };
    const std::vector<float> omega_values = { 1.0, 1.5 };

    std::vector<run_stats> stats;

    for (const auto omega : omega_values)
    {
        std::cout << "Omega: " << omega << std::endl;
        for (const auto k : k_values)
        {
            for (size_t i = 0; i < num_iter; ++i)
            {
                std::cout << "\r\tRunning for k = " << k << ". " << i << " / " << num_iter << "..." << std::flush;

                const auto matrix = generate(k);

                auto begin = std::chrono::steady_clock::now();
                branch_and_bound bb(matrix, k);
                bb.run();
                auto end = std::chrono::steady_clock::now();
                long bb_time = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();
                stats.push_back({
                    algorithm::bb,
                    bb.get_map().size(), bb_time,
                    k, omega
                });

                begin = std::chrono::steady_clock::now();
                divide_and_conquer dc(matrix, k);
                dc.run();
                end = std::chrono::steady_clock::now();
                long dc_time = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();
                stats.push_back({
                    algorithm::dc,
                    dc.get_map().size(), dc_time,
                    k, omega
                });

                begin = std::chrono::steady_clock::now();
                rappas rap(matrix, k);
                rap.run();
                end = std::chrono::steady_clock::now();
                long rap_time = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();
                stats.push_back({
                                        algorithm::rappas,
                                        rap.get_map().size(), rap_time,
                                        k, omega
                                });

                assert_equal(bb.get_map(), dc.get_map());
                assert_equal(bb.get_map(), rap.get_map());
            }
            std::cout << "\r\tRunning for k = " << k << ". Done." << std::endl;
        }
    }

    print_as_csv(stats, filename);
}

int main()
{
    srand(42);

    //test_one(10);
    //test_suite();

    benchmark(500, std::string(std::tmpnam(nullptr)) + ".csv");

    return 0;
}
