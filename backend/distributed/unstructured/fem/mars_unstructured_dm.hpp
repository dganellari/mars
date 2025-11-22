#pragma once

#include "mars.hpp"
#include "mars_unstructured_dof_handler.hpp"
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/copy.h>
#include <mpi.h>

namespace mars {
namespace fem {

/**
 * @brief Distributed Data Manager for unstructured FEM
 * 
 * Provides distributed data management similar to the structured DM,
 * but adapted for unstructured meshes using ElementDomain.
 */
template<typename DofHandler, typename RealType = float, typename AcceleratorTag = cstone::GpuTag>
class UnstructuredDM {
public:
    using DeviceVector = typename VectorSelector<RealType, AcceleratorTag>::type;
    
    UnstructuredDM(const DofHandler& dof_handler)
        : dof_handler_(dof_handler)
    {}
    
    /**
     * @brief Add a data field
     */
    template<typename DataType = RealType>
    void add_data_field(size_t size = 0) {
        data_fields_.emplace_back(size);
    }
    
    /**
     * @brief Get reference to data field
     */
    template<size_t Index>
    DeviceVector& get_data() {
        return data_fields_[Index];
    }
    
    template<size_t Index>
    const DeviceVector& get_data() const {
        return data_fields_[Index];
    }
    
    /**
     * @brief Resize all data fields
     */
    void resize(size_t size) {
        for (auto& field : data_fields_) {
            field.resize(size);
        }
    }
    
    /**
     * @brief Fill data field with value
     */
    template<size_t Index>
    void fill(RealType value) {
        thrust::fill(thrust::device_pointer_cast(data_fields_[Index].data()),
                    thrust::device_pointer_cast(data_fields_[Index].data() + data_fields_[Index].size()),
                    value);
    }
    
    /**
     * @brief Copy data from host vector to device field
     */
    template<size_t Index>
    void copy_from_host(const std::vector<RealType>& host_data) {
        thrust::copy(host_data.begin(), host_data.end(),
                    thrust::device_pointer_cast(data_fields_[Index].data()));
    }
    
    /**
     * @brief Copy data from device field to host vector
     */
    template<size_t Index>
    void copy_to_host(std::vector<RealType>& host_data) const {
        host_data.resize(data_fields_[Index].size());
        thrust::copy(thrust::device_pointer_cast(data_fields_[Index].data()),
                    thrust::device_pointer_cast(data_fields_[Index].data() + data_fields_[Index].size()),
                    host_data.begin());
    }
    
    /**
     * @brief Gather ghost data from neighboring ranks
     * This is a simplified version - in practice, this would need proper MPI communication
     */
    void gather_ghost_data() {
        // For now, this is a placeholder
        // In a full implementation, this would use the dof_handler's communication patterns
        // to gather ghost data from neighboring ranks
    }
    
    /**
     * @brief Scatter ghost data to neighboring ranks
     */
    void scatter_ghost_data() {
        // For now, this is a placeholder
        // In a full implementation, this would use the dof_handler's communication patterns
        // to scatter ghost data to neighboring ranks
    }
    
    /**
     * @brief Get the DOF handler
     */
    const DofHandler& get_dof_handler() const {
        return dof_handler_;
    }
    
private:
    const DofHandler& dof_handler_;
    std::vector<DeviceVector> data_fields_;
};

} // namespace fem
} // namespace mars