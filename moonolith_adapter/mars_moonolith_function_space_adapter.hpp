#ifndef MOONOLITH_MARS_FUNCTION_SPACE_ADAPTER_HPP
#define MOONOLITH_MARS_FUNCTION_SPACE_ADAPTER_HPP

#include "moonolith_mars_mesh_adapter.hpp"
#include "moonolith_function_space.hpp"
#include "mars_par_mesh.hpp"

namespace mars {

    template<class FunctionSpace>
    class MarsFunctionSpaceManager {
    public:
        using MeshManager = mars::MarsCollectionManager<typename FunctionSpace::Mesh>;
        using Elem = typename FunctionSpace::Mesh::Elem;

        template<class Adapter>
        static void data_added(const Integer local_idx, const Adapter &a)
        {
            MeshManager::data_added(local_idx, a);
        }

        static const Elem &elem(const FunctionSpace &space, const Integer handle) 
        {
            return MeshManager::elem(space.mesh(), handle);
        }

        static Integer tag(const FunctionSpace &space, const Integer handle)
        {
            return space.dof_map().block(handle);
        }

        static Integer n_elements(const FunctionSpace &space)
        {
            return MeshManager::n_elements(space.mesh());
        }

        static Integer elements_begin(const FunctionSpace &space)
        {
            return MeshManager::elements_begin(space.mesh());
        }

        static Integer elements_end(const FunctionSpace &space)
        {
            return MeshManager::elements_end(space.mesh());
        }

        static Integer handle(const FunctionSpace &space, const Integer element_index)
        {
            return MeshManager::handle(space.mesh(), element_index);
        }

        static bool skip(const FunctionSpace &space, const Integer &element_index)
        {
            return MeshManager::skip(space.mesh(), element_index);
        }

        template<class Bound>
        static void fill_bound(const FunctionSpace &space, const Integer handle, Bound &bound, const double blow_up)
        {
           MeshManager::fill_bound(space.mesh(), handle, bound, blow_up);
        }

        template<class Iterator>
        static void serialize(const FunctionSpace &space, const Iterator &begin, const Iterator &end, moonolith::OutputStream &os)
        {
            CHECK_STREAM_WRITE_BEGIN("fs_serialize", os);
            MeshManager::serialize(space.mesh(), begin, end, os);
            space.dof_map().serialize(begin, end, os);
            CHECK_STREAM_WRITE_END("fs_serialize", os);
        }

        static std::unique_ptr<FunctionSpace> build(moonolith::InputStream &is)
        {
            CHECK_STREAM_READ_BEGIN("fs_serialize", is);

            auto space = moonolith::make_unique<FunctionSpace>(  MeshManager::build(is) );
            space->dof_map().build(is);

            CHECK_STREAM_READ_END("fs_serialize", is);
            return space;
        }

    private:
        MeshManager mesh_manager_;
    };

    template<Integer Dim, Integer ManifoldDim, class Out>
    inline void make(const moonolith::FunctionSpace<mars::Mesh<Dim, ManifoldDim>> &space, const typename mars::Mesh<Dim, ManifoldDim>::Elem &elem, Out &out)
    {
        make(space.mesh(), elem, out);
    }

    template<Integer Dim, Integer ManifoldDim, class Out>
    inline void make(const moonolith::FunctionSpace<mars::IMesh<Dim>> &space, const typename mars::IMesh<Dim>::Elem &elem, Out &out)
    {
        make(space.mesh(), elem, out);
    }

    }

    namespace moonolith {
        template<Integer Dim, Integer ManifoldDim>
        class CollectionManager< FunctionSpace<mars::Mesh<Dim, ManifoldDim>> > : public mars::MarsFunctionSpaceManager< FunctionSpace<mars::Mesh<Dim, ManifoldDim>> > {};

        template<Integer Dim>
        class CollectionManager< FunctionSpace<mars::IMesh<Dim>> > : public mars::MarsFunctionSpaceManager< FunctionSpace<mars::IMesh<Dim>> > {};
    }

    namespace mars {
    
    template<Integer Dim, Integer ManifoldDim>
    std::unique_ptr< moonolith::FunctionSpace<mars::Mesh<Dim, ManifoldDim>> > make_function_space(
        mars::ParMesh<Dim, ManifoldDim> &par_mesh)
    {   
        mars::Mesh<Dim, ManifoldDim> &serial_mesh = par_mesh.get_serial_mesh();
        auto space = moonolith::make_unique< moonolith::FunctionSpace<mars::Mesh<Dim, ManifoldDim>> >(serial_mesh);
        auto &dof_map = space->dof_map();

        dof_map.resize(serial_mesh.n_elements());

        const moonolith::FEType type = ManifoldDim == 1? moonolith::EDGE2 :(ManifoldDim == 2? moonolith::TRI3 : (ManifoldDim == 3? moonolith::TET4 : moonolith::INVALID));

        for(Integer i = 0; i < serial_mesh.n_elements(); ++i) {
            const auto &e = serial_mesh.elem(i);
            auto &dof = dof_map.dof_object(i);
            dof.global_idx = par_mesh.elem_map().global(i);
            dof.block = e.block;
            dof.dofs.resize(mars::n_nodes(e));

            for(Integer k = 0; k < mars::n_nodes(e); ++k) {
                dof.dofs[k] = par_mesh.node_map().global(e.node(k));
            }

            dof.type = type;
        }

        return space;
    }
}

#endif //MOONOLITH_MARS_FUNCTION_SPACE_ADAPTER_HPP
