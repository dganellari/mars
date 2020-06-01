/* Copyright (c) 2016, Eidgenössische Technische Hochschule Zürich
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
   may be used to endorse or promote products derived from this software without
   specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. */

/* This example tries to do the same and compare to the step3 p4est example using MARS instead */

#include "mars_context.hpp"
#include <bits/c++config.h>
#include <exception>
#include <iostream>
#include <tuple>
#include <type_traits>
#include <utility>

#ifdef WITH_MPI

#include "mars_mpi_guard.hpp"

#ifdef WITH_KOKKOS
#include "mars_distributed_mesh_generation.hpp"
#include "mars_distributed_user_data.hpp"
#endif //WITH_KOKKOS
#endif

namespace mars
{

using Data = UserData<DistributedQuad4Mesh, double, double, double, double>;

template <Integer DIM>
struct ProblemDesc
{
    double center[DIM]; /**< coordinates of the center of
                                               the initial condition Gaussian
                                               bump */
    double bump_width;  /**< width of the initial condition
                                               Gaussian bump */
    double max_err;     /**< maximum allowed global
                                               interpolation error */
    double v[DIM];      /**< the advection velocity */

    Integer refine_period;      /**< the number of time steps
                                               between mesh refinement */
    Integer repartition_period; /**< the number of time steps
                                               between repartitioning */
    Integer write_period;       /**< the number of time steps
                                               between writing vtk files */
};

template <typename... T>
using user_tuple = mars::ViewsTuple<T...>;
/*
template <typename... T>
struct functor
{
    user_tuple<T...> tuple;

    functor(user_tuple<T...> t) : tuple(t) {}

    MARS_INLINE_FUNCTION
    void operator()(int i) const
    {
        [>note the use of std get instead<]
        std::get<1>(tuple)(i) = 1;
    }
};
 */
//from p4est step3 example
/* template <Integer DIM, Integer ...args>
MARS_INLINE_FUNCTION double initial_condition(double *x, double **du, const ProblemDesc<DIM> &pd)
{
    const double *c = pd.center;
    double bump_width = pd.bump_width;

    double r2, d[DIM];
    double arg, retval;

    r2 = 0.;

    std::cout<<"Dim: "<<DIM<<std::endl;
    for (int i = 0; i < DIM; ++i)
    {
        printf("i: %li, x: %d\n", i, x[i]);
        d[i] = x[i] - c[i];
        r2 += d[i] * d[i];
    }

    arg = -(1. / 2.) * r2 / bump_width / bump_width;
    retval = exp(arg);

    if (du)
    {
        for (int i = 0; i < DIM; ++i)
        {
            *du[i] = -(1. / bump_width / bump_width) * d[i] * retval;
        }
    }

    return retval;
}
 */
/* template<Integer DIM>
struct UpdateDu
{
    double bump_width;
    double retval;
    double d[DIM];


    user_tuple tuple;

    UpdateDu(double b, double r, double dd[DIM], user_tuple t) :
        bump_width(b), retval(r), d(dd), tuple(t) {}

    template<Integer I>
    void operator()(Integer i) const
    {
        std::get<I>(tuple) = -(1. / bump_width / bump_width) * d[i] * retval;
    }

};
 */
template <Integer DIM, Integer first, Integer second>
MARS_INLINE_FUNCTION double initial_condition(const Data &data, const int index, const ProblemDesc<DIM> &pd, double *x)
{
    const double *c = pd.center;
    double bump_width = pd.bump_width;

    double r2, d[DIM];
    double arg, retval;

    r2 = 0.;

    for (int i = 0; i < DIM; ++i)
    {
        /* printf("i: %i, x: %lf\n", i, x[i]); */
        d[i] = x[i] - c[i];
        r2 += d[i] * d[i];
    }

    arg = -(1. / 2.) * r2 / bump_width / bump_width;
    retval = exp(arg);

    data.get_elem_data<first>(index) = -(1. / bump_width / bump_width) * d[0] * retval;
    data.get_elem_data<second>(index) = -(1. / bump_width / bump_width) * d[1] * retval;
    /* apply_impl(UpdateDu<DIM>(bump_width, retval, d, user_data), std::forward_as_tuple(args...)); */

    printf("p:x %lf, py: %lf,  retval: %lf, du: %lf-%lf\n", x[0], x[1], retval, data.get_elem_data<first>(index), data.get_elem_data<second>(index));

    return retval;
}

//in case you might prefer the variadic version better you can use the Nthvalue equivalent to NthType to get
//the nth value from a nontype variadic template.
template <Integer DIM, Integer ...args>
MARS_INLINE_FUNCTION double initial_condition_variadic(const Data &data, const int index, const ProblemDesc<DIM> &pd, double *x)
{
    //in this case at least the DIM isqual to #args...
    assert(DIM == sizeof...(args));

    const double *c = pd.center;
    double bump_width = pd.bump_width;

    double r2, d[DIM];
    double arg, retval;

    r2 = 0.;

    for (int i = 0; i < DIM; ++i)
    {
        d[i] = x[i] - c[i];
        r2 += d[i] * d[i];
    }

    arg = -(1. / 2.) * r2 / bump_width / bump_width;
    retval = exp(arg);

    /* careful: one can not use a for loop for this as the values needs to be constexpr.
    The recursive approach on a parameter pack is the alternative. Check example at distributed_utils.hpp*/
    constexpr Integer first = NthValue<0, args...>::value;
    constexpr Integer second = NthValue<1, args...>::value;

    data.get_elem_data<first>(index) = -(1. / bump_width / bump_width) * d[0] * retval;
    data.get_elem_data<second>(index) = -(1. / bump_width / bump_width) * d[1] * retval;

    printf("p:x %lf, py: %lf,  retval: %lf, du: %lf-%lf\n", x[0], x[1], retval, data.get_elem_data<first>(index), data.get_elem_data<second>(index));

    return retval;
}


template <Integer I = 0, Integer N, Integer... Args>
typename std::enable_if<I == N, void>::type
MARS_INLINE_FUNCTION
for_each_du(const Data &data, const double bump_width, const double retval, const int index, const double *d)
{
}

template <Integer I = 0, Integer N, Integer... Args>
typename std::enable_if<I < N, void>::type
MARS_INLINE_FUNCTION
for_each_du(const Data &data, const double bump_width, const double retval, const int index, const double *d)
{
    constexpr Integer dataIndex = NthValue<I, Args...>::value;
    printf("val: %li\n", dataIndex);

    data.get_elem_data<dataIndex>(index) = -(1. / bump_width / bump_width) * d[I] * retval;
    for_each_du<I+1, N, Args...>(data, bump_width, retval, index, d);
}

//The third approach shows the recursive way of doing things for the non-type parameter pack
template <Integer DIM, Integer ...args>
MARS_INLINE_FUNCTION double initial_condition_recursive(const Data &data, const int index, const ProblemDesc<DIM> &pd, double *x)
{
    //in this case at least the DIM isqual to #args...
    assert(DIM == sizeof...(args));

    const double *c = pd.center;
    double bump_width = pd.bump_width;

    double r2, d[DIM];
    double arg, retval;

    r2 = 0.;

    for (int i = 0; i < DIM; ++i)
    {
        d[i] = x[i] - c[i];
        r2 += d[i] * d[i];
    }

    arg = -(1. / 2.) * r2 / bump_width / bump_width;
    retval = exp(arg);

    //compile time loop over data.get_elem_data<I>
    for_each_du<0, sizeof...(args), args...>(data, bump_width, retval, index, d);

    //just for printing purposes same as the previous example.
    constexpr Integer first = NthValue<0, args...>::value;
    constexpr Integer second = NthValue<1, args...>::value;

    printf("p:x %lf, py: %lf,  retval: %lf, du: %lf-%lf\n", x[0], x[1], retval, data.get_elem_data<first>(index), data.get_elem_data<second>(index));

    return retval;
}

template <Integer Type>
MARS_INLINE_FUNCTION
void get_midpoint_coordinates(double *point, const Integer sfc, const Integer xDim, const Integer yDim, const Integer zDim)
{
    assert(xDim != 0);
    assert(yDim != 0);

    get_vertex_coordinates_from_sfc<Type>(sfc, point, xDim, yDim, zDim);

    double hx = 1. / xDim;
    double hy = 1. / yDim;

    /* /2 for the midpoint */
    point[0] += hx/2;
    point[1] += hy/2;

    if (Type == ElementType::Hex8)
    {
        assert(zDim != 0);

        double hz = 1 / zDim;
        point[2] += hz/3;
    }
}

void advection(int &argc, char **&argv, const int level)
{

    using namespace mars;
    try
    {
        mars::proc_allocation resources;
        /*
        // try to detect how many threads can be run on this system
        resources.num_threads = marsenv::thread_concurrency();

        // override thread count if the user set MARS_NUM_THREADS
        if (auto nt = marsenv::get_env_num_threads())
        {
            resources.num_threads = nt;
        } */

#ifdef WITH_MPI
        // initialize MPI
        marsenv::mpi_guard guard(argc, argv, false);

        // assign a unique gpu to this rank if available
        /*  resources.gpu_id = marsenv::find_private_gpu(MPI_COMM_WORLD); */

        // create a distributed context
        auto context = mars::make_context(resources, MPI_COMM_WORLD);
        int proc_num = mars::rank(context);
#else
        // resources.gpu_id = marsenv::default_gpu();

        // // create a local context
        // auto context = mars::make_context(resources);
#endif

#ifdef WITH_KOKKOS

        DistributedQuad4Mesh mesh;
        generate_distributed_cube(context, mesh, level, level, 0);

        const Integer xDim = mesh.get_XDim();
        const Integer yDim = mesh.get_YDim();
        const Integer zDim = mesh.get_ZDim();

        constexpr Integer Dim = DistributedQuad4Mesh::Dim;

        using Elem = typename DistributedQuad4Mesh::Elem;
        constexpr Integer Type = Elem::ElemType;

        std::cout << "Type: " << Type<< std::endl;

        ProblemDesc<Dim> pd;
        pd.bump_width = 0.1;
        pd.max_err = 2.e-2;
        pd.center[0] = 0.5;
        pd.center[1] = 0.5;

        pd.v[0] = -0.445868402501118;
        pd.v[0] = -0.895098523991131;

       Data data(&mesh);

        ViewVectorType<Integer> sfc = mesh.get_view_sfc();
        /* another option to use this one with a functor instead of the lamda. Just be careful what you use within
         * the Lambda body since it is copied by value. Kokkos views are no problem but the userdata object has other
         * host containers which should not be used inside the lambda as it will copy the entire container's content.*/
        data.set_init_cond(MARS_LAMBDA(const int i) {
            /* data.get_elem_data<1>(i) = i; */
            double midpoint[Dim];
            get_midpoint_coordinates<Type>(midpoint, sfc(i), xDim, yDim, zDim);

            /* double* du[Dim];
            du[0] = &data.get_elem_data<1>(i);
            du[1] = &data.get_elem_data<2>(i);
            data.get_elem_data<0>(i) - the solution u
            data.get_elem_data<0>(i) = initial_condition<Dim, 1, 2>(midpoint, du, pd); */

            /* data.get_elem_data<0>(i) = initial_condition<Dim, 1, 2>(data, i, pd, midpoint); */
            /* data.get_elem_data<0>(i) = initial_condition_variadic<Dim, 1, 2>(data, i, pd, midpoint); */
            data.get_elem_data<0>(i) = initial_condition_recursive<Dim, 1, 2>(data, i, pd, midpoint);
        });

        /* second possibility to do it using a more general approach using a functor by coping
         * the tuple directly instead of the UserData object*/
        /* data.parallel_for_data(mesh.get_chunk_size(), functor<double, Integer, double>(data.get_user_data())); */

        /* the same as above but using a lambda instead of the functor. This way can be used for any
         * UserData class member accessing it using the data object. */
        /* data.parallel_for_data(
            mesh.get_chunk_size(), MARS_LAMBDA(const int i) {
                data.get_elem_data<1>(i) = 2;
            }); */

        exchange_ghost_user_data(context, data);

        data.print_nth_tuple<1>(proc_num);


        const Integer size_ch = mesh.get_chunk_size();
        ViewVectorType<double> max("max", 1);

        data.elem_iterate(MARS_LAMBDA(const int i) {
            if (data.get_elem_data<1>(i) > max(0))
            {
                max(0) = data.get_elem_data<1>(i);
            }

            if (i == size_ch - 1)
                printf("max: %lf\n", max(0));
        });

#endif
    }
    catch (std::exception &e)
    {
        std::cerr << "exception caught in ring miniapp: " << e.what() << "\n";
    }
}
} // namespace mars
