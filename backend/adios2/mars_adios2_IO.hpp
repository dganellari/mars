#ifndef MARS_ADIOS2_IO_HPP
#define MARS_ADIOS2_IO_HPP

#include "mars_utils_kokkos.hpp"

#include <memory>

namespace adios2 {
    class IO;
    class Engine;
}  // namespace adios2

namespace mars {
    namespace adios2 {

        template <class Mesh>
        class IO {
        public:
            IO(Mesh &mesh);
            ~IO();

            class Field {
            public:
                virtual ~Field() = default;

                Field(std::string name, int n_components) : name(name), n_components(n_components) {}

                virtual void define(::adios2::IO &io) = 0;
                virtual void set(::adios2::Engine &engine, ::adios2::IO &io) = 0;

                std::string name;
                int n_components{1};
            };

            class RealField : public Field {
            public:
                RealField(std::string name, int n_components, const ViewVectorType<Real> &data)
                    : Field(std::move(name), n_components), data(data) {}

                void define(::adios2::IO &io) override;
                void set(::adios2::Engine &engine, ::adios2::IO &io) override;

                ViewVectorType<Real> data;
            };

            void add_field(const std::string &name, const int n_components, const ViewVectorType<Real> &data);
            void set_output_path(const std::string &path);
            bool open_output();
            void close();
            bool write_step();

        private:
            class Impl;
            std::unique_ptr<Impl> impl_;
        };
    }  // namespace adios2
}  // namespace mars

#endif  // MARS_ADIOS2_IO_HPP