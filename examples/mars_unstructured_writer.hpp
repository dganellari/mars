#ifndef MARS_UNSTRUCTURED_WRITER_HPP
#define MARS_UNSTRUCTURED_WRITER_HPP

#include "adios2.h"

class UnstructuredWriter {
public:
    UnstructuredWriter(adios2::IO io);
    void open(const std::string &fname);
    void write(int step);
    void close();

protected:
    adios2::IO io;
    adios2::Engine writer;
    adios2::Variable<double> var_data;
    std::vector<double> data;
};
