#ifndef XPAS_ALGS_AR_H
#define XPAS_ALGS_AR_H

#include <string>
#include <vector>
#include <unordered_map>
#include "matrix.h"

using ar_result = std::unordered_map<std::string, matrix>;

class raxmlng_reader
{
public:
    raxmlng_reader(const std::string& file_name) noexcept;
    raxmlng_reader(const raxmlng_reader&) = delete;
    raxmlng_reader(raxmlng_reader&&) = delete;
    raxmlng_reader& operator=(const raxmlng_reader&) = delete;
    raxmlng_reader& operator=(raxmlng_reader&&) = delete;
    ~raxmlng_reader() noexcept = default;

    ar_result read();

private:
    ar_result read_matrix();

    std::string _file_name;
};

#endif //XPAS_ALGS_AR_H
