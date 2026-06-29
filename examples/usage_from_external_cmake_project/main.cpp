// Downstream-consumer smoke test for an installed MARS.
// Purpose: prove that an external CMake project can find_package(Mars), compile
// against the installed headers, and link Mars::mars. It deliberately does no
// real computation (no GPU, no MPI, no mesh) -- if this builds and links, the
// installed package is consumable, which is all the install gate needs to check.
#include <mars_config.hpp>

#include <iostream>

int main()
{
    std::cout << "MARS install smoke: find_package(Mars) + Mars::mars link OK\n";
    return 0;
}
