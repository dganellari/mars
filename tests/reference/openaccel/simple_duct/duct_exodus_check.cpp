// A duct_mesh.py Exodus file, read by the production native reader (read_simple_mesh: Exodus
// arrays, side-set selection, exterior coverage), must equal the C++ lattice of the host runs:
// coordinates within Lattice::coordinate_tolerance() (the same operations in both languages, so
// normally bitwise), connectivity exactly, and every boundary face with its inlet/outlet/wall
// kind exactly. Selecting the side sets under other names must be rejected.
//   duct_exodus_check FILE.exo --cells N [--length 7 --width 2 --height 1 --stretch 2]
#include "duct_mesh.hpp"
#include "mars_segregated_simple_native_mesh.hpp"
#include <iostream>
#include <set>
#include <tuple>
using namespace mars::segregated::runtime;

int main(int argc,char** argv) {
    MPI_Init(&argc,&argv); int rank; MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    try {
        ensure(argc>=4,"usage: duct_exodus_check FILE.exo --cells N [--length L --width W --height H --stretch S]");
        int cells=0; double length=7, width=2, height=1, stretch=2;
        for (int i=2;i+1<argc;i+=2) {
            const std::string key=argv[i]; const double value=std::stod(argv[i+1]);
            if (key=="--cells") cells=int(value); else if (key=="--length") length=value; else if (key=="--width") width=value;
            else if (key=="--height") height=value; else if (key=="--stretch") stretch=value; else throw std::runtime_error("unknown option "+key);
        }
        const auto lattice=duct::Lattice::make(cells,length,width,height,stretch);
        const auto expected=duct::global_input(lattice);
        const auto mesh=read_simple_mesh(MPI_COMM_WORLD,argv[1]);
        // Coordinates within the documented rounding bound (the file holds Python's doubles);
        // topology and boundary kinds exactly.
        double worst=0;
        bool same=mesh.x.size()==expected.x.size() && mesh.y.size()==expected.y.size() && mesh.z.size()==expected.z.size();
        for (std::size_t i=0;same && i<mesh.x.size();++i)
            worst=std::max({worst,std::abs(mesh.x[i]-expected.x[i]),std::abs(mesh.y[i]-expected.y[i]),std::abs(mesh.z[i]-expected.z[i])});
        same=same && worst<=lattice.coordinate_tolerance();
        for (int j=0;j<4;++j) same=same && mesh.nodes[j]==expected.nodes[j];
        std::set<std::tuple<int,int,int>> a,b;
        for (const auto& f:mesh.faces) a.insert({f.element,f.ordinal,f.kind});
        for (const auto& f:expected.faces) b.insert({f.element,f.ordinal,f.kind});
        same=same && a==b && a.size()==mesh.faces.size() && b.size()==expected.faces.size();
        simple_collective(MPI_COMM_WORLD,same,"Exodus mesh differs from the C++ duct lattice");
        bool rejected=false;
        try { read_simple_mesh(MPI_COMM_WORLD,argv[1],mars::segregated::SimpleBoundaryNames{"inlet","outlet",{"wall"}}); }
        catch (const std::exception&) { rejected=true; }
        simple_collective(MPI_COMM_WORLD,rejected,"unknown wall side-set name accepted");
        if (!rank) std::cout<<"PASS: production reader returns the C++ duct lattice ("<<mesh.x.size()<<" nodes, "
                            <<mesh.nodes[0].size()<<" Tet4, "<<mesh.faces.size()<<" boundary faces)\n";
    } catch (const std::exception& e) { std::cerr<<"FAIL: "<<e.what()<<'\n'; MPI_Abort(MPI_COMM_WORLD,1); }
    MPI_Finalize();
}
