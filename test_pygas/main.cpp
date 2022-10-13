#include "libthermo/gas.h"

#include <pybind11/pybind11.h>

#include <iostream>


namespace py = pybind11;
using namespace thermo;

struct Row
{
    Row(ThermoInterface* gas_)
    {
        gas = gas_;
        std::cout << gas->H(288);
    }

    ThermoInterface* gas;

    void print_h()
    {
        std::cout << gas->H(288.15);
    }
};


PYBIND11_MODULE(test_pygas, m)
{
    py::class_<Row>(m, "Row").def(py::init<ThermoInterface*>()).def("print_h", &Row::print_h);
}