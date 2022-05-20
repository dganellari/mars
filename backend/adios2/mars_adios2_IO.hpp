#ifndef MARS_ADIOS2_IO_HPP
#define MARS_ADIOS2_IO_HPP

#ifdef WITH_KOKKOS
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

                Field(std::string name, const Integer n_nodes, int n_components)
                    : name(name), n_nodes(n_nodes), n_components(n_components) {}

                virtual void define(::adios2::IO &io) = 0;
                virtual void set(const Mesh &mesh, ::adios2::Engine &engine, ::adios2::IO &io) = 0;

                std::string name;
                int n_nodes;
                int n_components{1};
            };

            class RealField : public Field {
            public:
                RealField(std::string name, const Integer n_nodes, int n_components, const ViewVectorType<Real> &data)
                    : Field(std::move(name), n_nodes, n_components), data(data) {}

                void define(::adios2::IO &io) override;
                void set(const Mesh &mesh, ::adios2::Engine &engine, ::adios2::IO &io) override;

                ViewVectorType<Real> data;
            };

            void add_field(const std::string &name, const int n_components, const ViewVectorType<Real> &data);
            void set_output_path(const std::string &path);
            bool open_output();
            void close();
            bool write_step();

            template <class View>
            bool write_tpetra(const std::string &path, const View &data) {
                ::mars::ViewVectorType<Real> x_view_rank1(::Kokkos::subview(data, ::Kokkos::ALL, 0));
                return aux_write(path, x_view_rank1);
            }

            template <class View>
            bool write(const std::string &path, const View &data) {
                return aux_write(path, data);
            }

        private:
            class Impl;
            std::unique_ptr<Impl> impl_;

            bool aux_write(const std::string &path, const ViewVectorType<Real> &data);
        };
    }  // namespace adios2
}  // namespace mars
#endif
#endif  // MARS_ADIOS2_IO_HPP