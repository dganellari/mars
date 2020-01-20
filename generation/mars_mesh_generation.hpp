#ifndef GENERATION_MARS_MESH_GENERATION_HPP_
#define GENERATION_MARS_MESH_GENERATION_HPP_
//libmesh similar approach: https://libmesh.github.io.

#include "mars_mesh.hpp"
#include "mars_globals.hpp"

#ifdef WITH_KOKKOS
#include "mars_mesh_kokkos.hpp"
#endif

namespace mars {
    
	/*different approach which works fine for the unit_cube and would have been faster
	but unfortunately it does not generalize for other cube dimensions.*/
	inline Mesh3 generate_unit_cube()
	{
		Mesh3 mesh;

		constexpr Integer xDim = 1, yDim = 1, zDim = 1, Dim = 3;

		const int n_elements = xDim * yDim * zDim * 24; //24 tetrahedrons on one hex27
		const int n_nodes = (2 * xDim + 1) * (2 * yDim + 1) * (2 * zDim + 1);
		mesh.reserve(n_elements, n_nodes);

		Vector<std::vector<Integer>, hex_n_sides> sides;

		for (Integer k = 0; k <= 2 * zDim; ++k)
		{
			for (Integer j = 0; j <= 2 * yDim; ++j)
			{
				for (Integer i = 0; i <= 2 * xDim; ++i)
				{
					Vector<Real, Dim> p(
							{ static_cast<Real>(i) / static_cast<Real>(2 * xDim),
									static_cast<Real>(j)
											/ static_cast<Real>(2 * yDim),
									static_cast<Real>(k)
											/ static_cast<Real>(2 * zDim), });

					mesh.add_point(p);

					Integer ind = index(1, 1, i, j, k);

					//build the faces on the fly
					if (k == 0)
						add_side(sides(0), i, j, ind);

					if (k == 2)
						add_side(sides(1), i, j, ind);

					//if (k == 2 * (zDim - 1))
					if (j == 0)
						add_side(sides(2), i, k, ind);

					//if (j == 2 * (yDim - 1))
					if (j == 2)
						add_side(sides(3), i, k, ind);

					if (i == 0)
						add_side(sides(4), j, k, ind);

					//if (i == 2 * (xDim - 1))
					if (i == 2)
						add_side(sides(5), j, k, ind);

				}
			}
		}

		//build tetrahedra elements from the hex27 faces.
		for (unsigned int i = 0; i < hex_n_sides; ++i)
		{
			Integer tmp = sides(i)[2]; //swap places between the last point and the midface point.
			sides(i)[2] = sides(i)[4];
			sides(i)[4] = tmp;

			for (unsigned int k = 0; k < 4; k++)
			{
				std::array<Integer, 4> nodes;

				nodes[0] = sides(i)[k];
				nodes[1] = sides(i)[4]; // midface point always the last element.
				nodes[2] = (k == 3 ? sides(i)[0] : sides(i)[k + 1]); // rotation to catch all combinations.
				nodes[3] = index(yDim, zDim, 1, 1, 1); // the center of the cube.

				mesh.add_elem(nodes);
			}
		}

		return mesh;
	}
    
    template<Integer Dim, Integer ManifoldDim>
    bool generate_cube(Mesh<Dim, ManifoldDim>& mesh, const Integer xDim,
                       const Integer yDim, const Integer zDim) {
        
        using Elem     = mars::Simplex<Dim, ManifoldDim>;
        
        assert(ManifoldDim <= Dim);
        assert(Dim <= 3);

        mesh.clear();
        
        switch (ManifoldDim) {
                
                /*case 0: {
                 
                 assert(xDim == 0);
                 assert(yDim == 0);
                 assert(zDim == 0);
                 
                 mesh.reserve(1, 1);
                 
                 Vector<Real, Dim> p( { 0.0, 0.0, 0.0 });
                 
                 int index = mesh.add_point(p);
                 
                 Simplex<Dim, ManifoldDim> e;
                 e.node(index);
                 
                 mesh.add_elem(e);
                 
                 return true;
                 }*/
                
            case 1: {
                
                assert(xDim != 0);
                assert(yDim == 0);
                assert(zDim == 0);
                
                const int n_elements = xDim;
                const int n_nodes = xDim + 1;
                mesh.reserve(n_elements, n_nodes);

                Vector<Real, Dim> p;
                p.zero();
                
                for (Integer i = 0; i <= xDim; ++i) {
                    p[0] = static_cast<Real>(i) / static_cast<Real>(xDim);
                    mesh.add_point(p);
                }
                
                for (Integer i = 0; i < xDim; ++i) {
                    
                    std::array<Integer, ManifoldDim + 1> nodes;
                    
                    nodes[0] = i;
                    nodes[1] = i + 1;
                    
                    mesh.add_elem(nodes);
                }
                
                return true;
            }
                
            case 2: {
                
                assert(xDim > 0);
                assert(yDim > 0);
                assert(zDim == 0);
                
                const int n_elements = 2 * xDim * yDim;
                const int n_nodes = (xDim + 1) * (yDim + 1);
                mesh.reserve(n_elements, n_nodes);
                
                Vector<Real, Dim> p;
                p.zero();

                for(Integer i = 0; i <= xDim; ++i) {
                    for(Integer j = 0; j <= yDim; ++j) {
                        
                        p[0] = static_cast<Real>(i) / static_cast<Real>(xDim);
                        p[1] = static_cast<Real>(j) / static_cast<Real>(yDim);
                        
                        mesh.add_point(p);
                        assert(mesh.n_nodes() <= n_nodes);
                    }
                }
                
                const int offset = yDim + 1;
                
                for (Integer i = 0; i < xDim; ++i) {
                    for (Integer j = 0; j < yDim; ++j) {
                        
                        std::array<Integer, ManifoldDim + 1> nodes;
                        
                        nodes[0] = i * offset + j;
                        nodes[1] = (i + 1) * offset + j;
                        nodes[2] = (i + 1) * offset + (j + 1);
                        
                        mesh.add_elem(nodes);
                        
                        //nodes[0] = i * offset + j;
                        nodes[1] = (i + 1) * offset + (j + 1);
                        nodes[2] = i * offset + (j + 1);

                        mesh.add_elem(nodes);
                        
                    }
                }
                
                return true;
            }
            case 3: { //building the tetra mesh from the hex27 elem
                
                assert(xDim != 0);
                assert(yDim != 0);
                assert(zDim != 0);
                
                const int n_elements = xDim * yDim * zDim * 24; //24 tetrahedrons on one hex27
                const int n_nodes = (2 * xDim + 1) * (2 * yDim + 1) * (2 * zDim + 1);
                
                const int n_tetra_nodes = 5 * (xDim * yDim * zDim)
                + 2 * (xDim * yDim + xDim * zDim + yDim * zDim)
                + (xDim + yDim + zDim) + 1;
                
                std::vector<bool> active_nodes(n_nodes);
                
                // mesh points need resize because the point vector is getting accessed using the [] operator.
                //First the center node is added somewhere at the end of the vector and then all the other nodes are added to it.
                mesh.reserve_elements(n_elements);
                mesh.resize_points(n_tetra_nodes);
                
                std::array<Integer, hex_side_n_nodes> side;

                Vector<Real, Dim> p;
                p.zero();
                
                int el_id = 1;
                for (Integer k = 0; k < 2 * zDim; k += 2) {
                    for (Integer j = 0; j < 2 * yDim; j += 2) {
                        for (Integer i = 0; i < 2 * xDim; i += 2) {
                            
                            //build the hex27 element which serves to generalize the idea to many hex27.
                            //without generating the hex27 element first there is no easy way to create the sides.
                            //Locally, using local element indexing it is possible.
                            std::array<Integer, hex_n_nodes> hex;
                            build_hex27(hex, xDim, yDim, i, j, k);
                            
                            //add center of the hex to the new points.
                            int centerHex = n_nodes / 2 + el_id;
                            
                            p[0] = static_cast<Real>(i+1) / static_cast<Real>(2 * xDim);
                            p[1] = static_cast<Real>(j+1) / static_cast<Real>(2 * yDim);
                            p[2] = static_cast<Real>(k+1) / static_cast<Real>(2 * zDim);
                            
                            mesh.point(centerHex) = p;
                            
                            //build tetrahedra elements from the hex27 faces.
                            for (unsigned int i = 0; i < hex_n_sides; ++i) {
                                
                                //use the hex27 element local connectivity to build the tetrahedra connectivity
                                for (unsigned int j = 0; j < hex_side_n_nodes; ++j) {
                                    side[j] = hex[hex_side_nodes[i][j]];
                                }
                                
                                //build 4 tetrahedra out of one face.
                                for (unsigned int k = 0; k < 4; k++) {
                                    
                                    auto& e = mesh.add_elem();
                                    
                                    e.nodes[0] = side[k] / 2;
                                    e.nodes[1] = side[8] / 2; // mid-face point always the last element.
                                    e.nodes[2] = (k == 3 ? side[0] : side[k + 1]) / 2; // rotation to catch all combinations.
                                    e.nodes[3] = centerHex; // the center of the cube.
                                    
                                    active_nodes[side[k]] = true;
                                    active_nodes[side[8]] = true;
                                    active_nodes[(k == 3 ? side[0] : side[k + 1])] =
                                    true;
                                    //active_nodes[hex[26]]=true;
                                }
                                
                            }
                            ++el_id;
                        }
                    }
                }
                
                //first the element indices as above
                //and then at this moment add the only the needed points to avoid extra node removal.
                for (Integer k = 0; k <= 2 * zDim; ++k) {
                    for (Integer j = 0; j <= 2 * yDim; ++j) {
                        for (Integer i = 0; i <= 2 * xDim; ++i) {
                            
                            Integer in = index(xDim, yDim, i, j, k);
                            
                            if (active_nodes[in]) {
                                
                                mesh.point(in / 2)[0] = static_cast<Real>(i)
                                / static_cast<Real>(2 * xDim);
                                mesh.point(in / 2)[1] = static_cast<Real>(j)
                                / static_cast<Real>(2 * yDim);
                                mesh.point(in / 2)[2] = static_cast<Real>(k)
                                / static_cast<Real>(2 * zDim);
                            }
                            
                        }
                    }
                }
                
                //std::cout<<"Used Memory: "<<generation::memory::get_physical_memory()<<std::endl;
                
                return true;
            }
            default: {
                
                std::cerr << "Not yet implemented for other dimensions." << std::endl;
                return false;
            }
        }
    }


    inline bool generate_line(Mesh1& mesh, const Integer xDim) {
        return generate_cube(mesh, xDim, 0, 0);
    }

    inline bool generate_square(Mesh2& mesh, const Integer xDim, const Integer yDim) {
        return generate_cube(mesh, xDim, yDim, 0);
    }

	//non simplex cube generation.
	template<Integer Dim, Integer ManifoldDim, Integer Type>
	bool generate_cube(
			Mesh<Dim, ManifoldDim, DefaultImplementation, NonSimplex<Type>>& mesh,
			const Integer xDim, const Integer yDim, const Integer zDim)
	{

		using Elem = mars::NonSimplex<Type>;

		assert(ManifoldDim <= Dim);
		assert(Dim <= 3);

		mesh.clear();

		switch (ManifoldDim)
		{
		case 2:
		{
			switch (Type)
			{
			case ElementType::Quad4:
			{
				assert(xDim > 0);
				assert(yDim > 0);
				assert(zDim == 0);

				const int n_elements = xDim * yDim;
				const int n_nodes = (xDim + 1) * (yDim + 1);
				mesh.reserve(n_elements, n_nodes);

				Vector<Real, Dim> p;
				p.zero();

				for (Integer j = 0; j <= yDim; ++j)
				{
					for (Integer i = 0; i <= xDim; ++i)
					{
						p[0] = static_cast<Real>(i) / static_cast<Real>(xDim);
						p[1] = static_cast<Real>(j) / static_cast<Real>(yDim);

						mesh.add_point(p);
						assert(mesh.n_nodes() <= n_nodes);
					}
				}

				const int offset = xDim + 1;

				for (Integer j = 0; j < yDim; ++j)
				{
					for (Integer i = 0; i < xDim; ++i)
					{
						std::array<Integer, Type> nodes;

						nodes[0] = i + offset * j;
						nodes[1] = (i+1) + offset * j;
						nodes[2] = (i+1) + offset * (j+1);
						nodes[3] = i + offset * (j+1);

						mesh.add_elem(nodes);
					}
				}

				return true;
			}
			default:
			{
				std::cerr << "Not yet implemented for other 2D element types."
						<< std::endl;
				return false;
			}
			}
		}
		case 3:
		{
			switch (Type)
			{
			case ElementType::Hex8:
			{
				assert(xDim != 0);
				assert(yDim != 0);
				assert(zDim != 0);

				const int n_elements = xDim * yDim * zDim;
				const int n_nodes = (xDim + 1) * (yDim + 1)	* (zDim + 1);

				Vector<Real, Dim> p;
				p.zero();

				for (Integer k = 0; k <= zDim; ++k)
				{
					for (Integer j = 0; j <= yDim; ++j)
					{
						for (Integer i = 0; i <= xDim; ++i)
						{
							p[0] = static_cast<Real>(i) / static_cast<Real>(xDim);
							p[1] = static_cast<Real>(j) / static_cast<Real>(yDim);
							p[2] = static_cast<Real>(k) / static_cast<Real>(zDim);

							mesh.add_point(p);
							assert(mesh.n_nodes() <= n_nodes);
						}
					}
				}

				const int xoffset = xDim + 1;
				const int yoffset = yDim + 1;

				for (Integer k = 0; k < zDim; ++k)
				{
					for (Integer j = 0; j < yDim; ++j)
					{
						for (Integer i = 0; i < xDim; ++i)
						{
							std::array<Integer, Type> nodes;

							nodes[0] = elem_index(i, j, k, xDim, yDim);
							nodes[1] = elem_index(i+1, j, k, xDim, yDim);
							nodes[2] = elem_index(i+1, j+1, k, xDim, yDim);
							nodes[3] = elem_index(i, j+1, k, xDim, yDim);
							nodes[4] = elem_index(i, j, k+1, xDim, yDim);
							nodes[5] = elem_index(i+1, j, k+1, xDim, yDim);
							nodes[6] = elem_index(i+1, j+1, k+1, xDim, yDim);
							nodes[7] = elem_index(i, j+1, k+1, xDim, yDim);

							mesh.add_elem(nodes);
						}
					}
				}

				return true;
			}
			default:
			{
				std::cerr << "Not yet implemented for other 3D element types."
						<< std::endl;
				return false;
			}
			}
		}
		default:
		{
			std::cerr << "Not yet implemented for other dimensions." << std::endl;
			return false;
		}
		}
	}


    template<Integer Type>
	inline bool generate_square(NSMesh2<Type>& mesh, const Integer xDim, const Integer yDim) {
		return generate_cube(mesh, xDim, yDim, 0);
	}
    
    template<Integer Dim, Integer ManifoldDim, class KokkosImplementation>
    bool generate_cube(Mesh<Dim, ManifoldDim, KokkosImplementation>& mesh, const Integer xDim,
                       const Integer yDim, const Integer zDim) {
        
        using Elem = mars::Simplex<Dim, ManifoldDim>;
        
        assert(Dim <= 3);
        assert(ManifoldDim <= Dim);
        
        bool gen_pts = mesh.generate_points(xDim, yDim, zDim);
        
        bool gen_elm = mesh.generate_elements(xDim, yDim, zDim);

        if(!gen_pts || !gen_elm)
            std::cerr << "Not implemented for other dimensions yet" << std::endl;

        return (gen_pts && gen_elm);
    }
    
		//non simplex cube generation.
	template<Integer Dim, Integer ManifoldDim, Integer Type>
	bool generate_cube(
			Mesh<Dim, ManifoldDim, KokkosImplementation, NonSimplex<Type, KokkosImplementation>>& mesh,
			const Integer xDim, const Integer yDim, const Integer zDim)
	{

		using Elem = mars::NonSimplex<Type, KokkosImplementation>;

		assert(Dim <= 3);
        assert(ManifoldDim <= Dim);

        bool gen_pts = mesh.generate_points(xDim, yDim, zDim, Type);
        
        bool gen_elm = mesh.generate_elements(xDim, yDim, zDim, Type);
		
		if(!gen_pts || !gen_elm)
			std::cerr << "Not implemented for other dimensions yet" << std::endl;

		return (gen_pts && gen_elm);
	}
    
#ifdef WITH_KOKKOS
    
    inline bool generate_square(ParallelMesh2& mesh, const Integer xDim, const Integer yDim) {
        return generate_cube(mesh, xDim, yDim, 0);
    }
    
    inline bool generate_line(ParallelMesh1& mesh, const Integer xDim) {
        return generate_cube(mesh, xDim, 0, 0);
    }
    
#endif
    
}

#endif /* GENERATION_MARS_MESH_GENERATION_HPP_ */
