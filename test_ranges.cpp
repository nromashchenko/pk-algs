#include <iostream>
#include <iomanip>
#include <range/v3/algorithm/for_each.hpp>
#include "matrix.h"

auto print_elem = [](auto elem){
    std::cout << std::fixed << std::setprecision(8) << elem << "\t";
    //std::cout << elem << ' ';
};

auto print = [](auto rng) {
    rs::for_each(rng, print_elem);
    std::cout << std::endl;
};

auto print2D = [](auto rng) {
    for (auto r: rng)
    {
        print(r);
    }
};

int main()
{
    const size_t num_cols = 10;

    auto data = generate(num_cols);
    //print2D(data.get_rows(0, 4));

    print2D(data.get_columns(0, num_cols));

    return 0;
}