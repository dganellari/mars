namespace mars
{
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

template<typename DstTuple, typename SrcTuple, std::size_t I = 0>
void copyTupleElements(DstTuple& dst, const SrcTuple& src)
{
    if constexpr (I < std::tuple_size_v<DstTuple>)
    {
        using DstType = typename std::tuple_element_t<I, DstTuple>::value_type;
        using SrcType = typename std::tuple_element_t<I, SrcTuple>::value_type;
        
        if constexpr (std::is_same_v<DstType, SrcType>)
        {
            // Same types - direct assignment
            std::get<I>(dst) = std::get<I>(src);
        }
        else
        {
            // Different types - convert element by element
            const auto& srcVec = std::get<I>(src);
            auto& dstVec = std::get<I>(dst);
            
            std::vector<DstType> converted(srcVec.size());
            std::transform(srcVec.begin(), srcVec.end(), converted.begin(),
                          [](const auto& val) { return static_cast<DstType>(val); });
            dstVec = converted;
        }
        
        copyTupleElements<DstTuple, SrcTuple, I + 1>(dst, src);
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
} // namespace mars