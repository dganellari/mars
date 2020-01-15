#ifndef MARS_SIGNED_NORMAL_HPP
#define MARS_SIGNED_NORMAL_HPP


#include "mars_base.hpp"
#include "mars_simplex.hpp"
#include "mars_array.hpp"


namespace mars{





///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////// class SignedNormal: for each face f of the element e                             ///////////////////
////////////// normal[e][f] returns the normal, outward[e][f] returns the sign of the normal    ///////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename Elem>
class SignedNormal;


template<Integer Dim, Integer ManifoldDim>
class SignedNormal<Simplex<Dim,ManifoldDim>>{
public:
    using NormalType=Array< Vector< Real, Dim> , ManifoldDim + 1 >;
    using SignType=Array<Real, ManifoldDim + 1 >;
    
    auto& operator () () const { return normal_; };
     auto& normals() const { return normal_; };
    // const Vector< Vector< Real, Dim> , ManifoldDim + 1 >& normal(const Integer elem_id) const {return normal_[elem_id];};
    auto& sign() const {return outward_;};
    const auto& sign(const Integer elem_id) const {return outward_[elem_id];};
    
    template< typename MeshT>
    void init(MeshT& mesh)
    {
        using Elem=typename MeshT::Elem;
        static_assert(Dim==ManifoldDim, "SignedNormal: computing internal normals of a mesh makes sense only if Dim==ManifoldDim");
        Integer n_elements=mesh.n_elements();
        Integer facenodes_carray[ManifoldDim];
        std::array<Integer,ManifoldDim+1> elemnodes_local;
        std::array<Integer,ManifoldDim+1> elemnodes_global;
        std::array<Integer,ManifoldDim+1> elemnodes_local_sort;
        
        std::array< Integer, ManifoldDim> facenodes_local;
        std::array< Integer, ManifoldDim> facenodes_local_sort;
        std::array<Integer,ManifoldDim> facenodes_global;
        std::array<Integer,ManifoldDim> facenodes_global_sort;
        
        std::array<Vector<Real,Dim>, ManifoldDim+1> points;
        Vector<Vector<Real,Dim>, ManifoldDim> facepoints;
        Vector<Vector<Real,Dim>, ManifoldDim+1> elempoints_sort;
        
        Vector<Real,Dim> facepoint_mean;
        Vector<Real,Dim> elempoint_mean;
        
        // we initialize the simplex
        Simplex<Dim,ManifoldDim> simplex_elem;
        Simplex<Dim, ManifoldDim-1> simplex_side;

        mesh.update_dual_graph();

        for(Integer nn=0;nn<ManifoldDim+1;nn++)
            simplex_elem.nodes[nn]=nn;
        
        normal_.resize(n_elements);
        outward_.resize(n_elements);
        
        
        // elemnodes_local is the array{ManifoldDim,ManifoldDim-1,...,1,0}
        for(Integer nn=0;nn<ManifoldDim+1 ;nn++)
            elemnodes_local[nn]=nn;
        std::reverse(std::begin(elemnodes_local), std::end(elemnodes_local));
        
        // loop on all the elements
        for(Integer ee=n_elements_;ee<mesh.n_elements();ee++)
        {
            auto &e = mesh.elem(ee);
      
            auto &adj = mesh.dual_graph().adj(ee);
            
            for(Integer nn=0;nn<ManifoldDim+1;nn++)
                elempoints_sort[nn]=mesh.points()[e.nodes[nn]];

            elempoint_mean=elempoints_sort.Tmean();


            for(Integer k = 0; k < n_sides(e); ++k) {
            // we use side and not side_sorted because we need normals
                e.side(k, simplex_side);

                auto n = normal(simplex_side, mesh.points());

                if(adj[k]>ee && adj[k]!=INVALID_INDEX)
                {
                  normal_[ee][k]=-n;
                  outward_[ee][k]=-1;
                }
                else
                {
                    normal_[ee][k]=n;
                    outward_[ee][k]=1;
                }
            }



        }

        n_elements_=mesh.n_elements();
    }
    
    SignedNormal():
    n_elements_(0)
    {}
    
    template< typename MeshT>
    SignedNormal(MeshT& mesh):
    n_elements_(0)
    {init(mesh);}
    
    template< typename MeshT>
    void print(const MeshT& mesh)
    {
        using Elem=typename MeshT::Elem;
        Vector<Vector<Real,Dim>, ManifoldDim+1> elempoints_sort;
        std::array<Integer,ManifoldDim+1> nodes;
        std::vector<Vector<Real,Dim>> points;
        Integer facenodes_carray[ManifoldDim];
        std::array< Integer, ManifoldDim> facenodes_local;
        
        for(Integer ee=0;ee<normal_.size();ee++)
        {
            // std::cout<<std::endl<<"ELEMENT ID == "<<ee<<" WITH NODES"<<std::endl;
            Elem elem=mesh.elem(ee);
            nodes=elem.nodes;
            // for(Integer mm=0;mm<ManifoldDim+1;mm++)
            //     std::cout<< nodes[mm] <<" ";
            // for(Integer mm=0;mm<ManifoldDim+1;mm++)
            //     std::cout<< mesh.points()[nodes[mm]] <<" ";
            // std::cout<<std::endl;
            
            for(Integer mm=0;mm<ManifoldDim+1;mm++)
            {
                Combinations<ManifoldDim + 1,ManifoldDim>::generate(mm,facenodes_carray);
                std::copy(std::begin(facenodes_carray), std::end(facenodes_carray), std::begin(facenodes_local));
                // std::cout<<"facenodes_local== "<<std::endl;
                // for(Integer ii=0;ii<ManifoldDim;ii++)
                //     std::cout<<facenodes_local[ii]<<" ";
                // std::cout<<"normal== ";
                // normal_[ee][mm].describe(std::cout);
                // std::cout<<"outward/inward(+1/-1)== "<<outward_[ee][mm]<<std::endl;
            }
        }
    }
    
 private:
    
    std::vector< NormalType >  normal_;
    std::vector< SignType > outward_;
    Integer n_elements_;   
};


}
#endif