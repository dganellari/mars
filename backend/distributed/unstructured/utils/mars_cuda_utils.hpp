#define cudaCheckError()                                                                 \
    {                                                                                    \
        cudaError_t e = cudaGetLastError();                                              \
        if (e != cudaSuccess) {                                                          \
            printf("CUDA error %s:%d: %s\n", __FILE__, __LINE__, cudaGetErrorString(e)); \
            exit(EXIT_FAILURE);                                                          \
        }                                                                                \
    }

// Recursive compile-time tuple copy
template<size_t I = 0, typename Tuple1, typename Tuple2>
void copyTupleElements(Tuple1& dst, const Tuple2& src)
{
    if constexpr (I < std::tuple_size_v<Tuple1>)
    {
        std::get<I>(dst) = std::get<I>(src);
        copyTupleElements<I + 1>(dst, src);
    }
}

// Helper for accessing tuple elements by runtime index
template <typename Tuple, std::size_t... Is>
decltype(auto) get_tuple_element_impl(Tuple& t, std::size_t i, std::index_sequence<Is...>) {
    decltype(auto) result = std::get<0>(t); // Default value (will be overwritten)
    
    bool found = false;
    (void)std::initializer_list<int>{
        (i == Is 
            ? (found = true, result = std::get<Is>(t), 0) 
            : 0)...
    };
    
    if (!found) {
        throw std::out_of_range("Tuple index out of range");
    }
    
    return result;
}

// Public interface that creates the right index sequence
template <typename Tuple>
decltype(auto) getTupleElement(Tuple& t, std::size_t i) {
    constexpr std::size_t size = std::tuple_size_v<std::remove_reference_t<Tuple>>;
    return get_tuple_element_impl(t, i, std::make_index_sequence<size>{});
}