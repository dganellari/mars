#ifndef MARS_ADIOS2_IO_HPP
#define MARS_ADIOS2_IO_HPP

#include <memory>

namespace mars {
    namespace adios2 {

        template <class Mesh>
        class IO {
        public:
            IO(Mesh &mesh);
            ~IO();

        private:
            class Impl;
            std::unique_ptr<Impl> impl_;
        };
    }  // namespace adios2
}  // namespace mars

#endif  // MARS_ADIOS2_IO_HPP