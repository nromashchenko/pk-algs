#include <iostream>
#include <fast-cpp-csv-parser/csv.h>
#include "ar.h"
#include "matrix.h"


raxmlng_reader::raxmlng_reader(const std::string& file_name) noexcept
    : _file_name{ file_name }
{}

ar_result raxmlng_reader::read()
{
    try
    {
        std::cout << "Loading RAXML-NG results: " + _file_name << "..." << std::endl;
        auto matrices = read_matrix();

        std::cout << "Loaded " << matrices.size() << " matrices" << std::endl;
        return matrices;
    }
    catch (::io::error::integer_overflow& error)
    {
        throw std::runtime_error("RAXML-NG result parsing error: " + std::string(error.what()));
    }
}

ar_result raxmlng_reader::read_matrix()
{
    // column-based
    ar_result result;

    ::io::CSVReader<5,
            ::io::trim_chars<' '>,
            ::io::no_quote_escape<'\t'>,
            ::io::throw_on_overflow,
            ::io::single_and_empty_line_comment<'.'>> _in(_file_name);

        _in.read_header(::io::ignore_extra_column, "Node", "p_A", "p_C", "p_G", "p_T");

        std::string node_label;
        score_t a, c, g, t;
        while (_in.read_row(node_label, a, c, g, t))
        {
            auto new_column = std::vector<score_t>{ a, c, g, t };
            result[node_label].get_data().push_back(new_column);
        }

        return result;
}