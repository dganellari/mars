#ifndef MARS_ENV_HPP
#define MARS_ENV_HPP

#include <memory>
#include "mars_base.hpp"

namespace mars {
    class Env {
    public:
        class Impl;

        Env(int argc, char *argv[]);
        ~Env();
        Integer exit_code();

    private:
        std::unique_ptr<Impl> impl_;
    };
}  // namespace mars

#endif  // MARS_ENV_HPP
