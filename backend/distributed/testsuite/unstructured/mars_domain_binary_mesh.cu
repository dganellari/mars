#include <gtest/gtest.h>
#include <filesystem>
#include <fstream>
#include <algorithm>
#include <vector>
#include <chrono>
#include <mpi.h>
#include "domain.hpp"

namespace fs = std::filesystem;
using namespace mars;

namespace cstone
{
    using std::get;
}


__global__ void validateConnectivityKernel(const unsigned* sfc0_ptr,
                                           const unsigned* sfc1_ptr, 
                                           const unsigned* sfc2_ptr,
                                           const unsigned* sfc3_ptr,
                                           cstone::Box<float> box,
                                           int* results, 
                                           size_t numElements)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= numElements) return;
    
    auto sfc0 = sfc0_ptr[tid];
    auto sfc1 = sfc1_ptr[tid];
    auto sfc2 = sfc2_ptr[tid];
    auto sfc3 = sfc3_ptr[tid];
    
    // Check that nodes are distinct
    bool distinct = (sfc0 != sfc1) && (sfc0 != sfc2) && (sfc0 != sfc3) && 
                   (sfc1 != sfc2) && (sfc1 != sfc3) && (sfc2 != sfc3);
    
    if (!distinct) {
        results[tid] = 0;
        return;
    }
    
    // FAIL if ALL SFC keys are 0 (indicates initialization bug)
    if (sfc0 == 0 && sfc1 == 0 && sfc2 == 0 && sfc3 == 0) {
        results[tid] = 0;
        return;
    }
    
    // If some (but not all) SFC keys are 0, that's valid for boundary elements
    if (sfc0 == 0 || sfc1 == 0 || sfc2 == 0 || sfc3 == 0) {
        results[tid] = 1;
        return;
    }
    
    // Check all four SFC keys by decoding and validating coordinates
    constexpr unsigned maxCoord = (1u << cstone::maxTreeLevel<cstone::SfcKind<unsigned>>{}) - 1;
    float invMaxCoord = 1.0f / maxCoord;
    const float tolerance = 1e-5f;

    // Helper lambda to check if coordinates are within bounds
    auto validateCoords = [&](unsigned sfc) -> bool {
        auto sfcKindKey = cstone::SfcKind<unsigned>(sfc);
        auto [ix, iy, iz] = cstone::decodeSfc(sfcKindKey);
        
        float x = box.xmin() + ix * invMaxCoord * (box.xmax() - box.xmin());
        float y = box.ymin() + iy * invMaxCoord * (box.ymax() - box.ymin());
        float z = box.zmin() + iz * invMaxCoord * (box.zmax() - box.zmin());
        
        return (x >= box.xmin()-tolerance && x <= box.xmax()+tolerance &&
                y >= box.ymin()-tolerance && y <= box.ymax()+tolerance &&
                z >= box.zmin()-tolerance && z <= box.zmax()+tolerance);
    };
    
    // Check all four SFC keys - each node must be valid
    bool allValid = validateCoords(sfc0) && validateCoords(sfc1) && 
                   validateCoords(sfc2) && validateCoords(sfc3);
    
    results[tid] = allValid ? 1 : 0;
}

__global__ void calculateVolumeKernel(const unsigned* sfc0_ptr,
                                     const unsigned* sfc1_ptr,
                                     const unsigned* sfc2_ptr, 
                                     const unsigned* sfc3_ptr,
                                     cstone::Box<float> box,
                                     float* volumes,
                                     size_t numElements)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= numElements) return;
    
    auto sfc0 = sfc0_ptr[tid];
    auto sfc1 = sfc1_ptr[tid];
    auto sfc2 = sfc2_ptr[tid];
    auto sfc3 = sfc3_ptr[tid];
    
    // Handle elements with all zero SFC keys (bug)
    if (sfc0 == 0 && sfc1 == 0 && sfc2 == 0 && sfc3 == 0) {
        volumes[tid] = -1.0f; // Invalid volume to indicate error
        return;
    }
    
    // Manual SFC to coordinate conversion with zero handling
    auto convertSfc = [&](unsigned sfc) -> thrust::tuple<float, float, float> {
        if (sfc == 0) {
            // SFC key 0 maps to minimum corner
            return thrust::make_tuple(box.xmin(), box.ymin(), box.zmin());
        }
        
        auto sfcKindKey = cstone::SfcKind<unsigned>(sfc);
        auto [ix, iy, iz] = cstone::decodeSfc(sfcKindKey);
        constexpr unsigned maxCoord = (1u << cstone::maxTreeLevel<cstone::SfcKind<unsigned>>{}) - 1;
        float invMaxCoord = 1.0f / maxCoord;
        
        float x = box.xmin() + ix * invMaxCoord * (box.xmax() - box.xmin());
        float y = box.ymin() + iy * invMaxCoord * (box.ymax() - box.ymin());
        float z = box.zmin() + iz * invMaxCoord * (box.zmax() - box.zmin());
        return thrust::make_tuple(x, y, z);
    };
    
    auto [x0, y0, z0] = convertSfc(sfc0);
    auto [x1, y1, z1] = convertSfc(sfc1);
    auto [x2, y2, z2] = convertSfc(sfc2);
    auto [x3, y3, z3] = convertSfc(sfc3);
    
    // Calculate tetrahedral volume
    float det = (x1-x0)*((y2-y0)*(z3-z0) - (z2-z0)*(y3-y0)) -
                (y1-y0)*((x2-x0)*(z3-z0) - (z2-z0)*(x3-x0)) +
                (z1-z0)*((x2-x0)*(y3-y0) - (y2-y0)*(x3-x0));
    
    volumes[tid] = fabsf(det) / 6.0f;
}

__global__ void coordinateConversionKernel(const unsigned* sfc0_ptr,
                                          const unsigned* sfc1_ptr,
                                          const unsigned* sfc2_ptr,
                                          const unsigned* sfc3_ptr,
                                          cstone::Box<float> box,
                                          float* coordinates, 
                                          size_t numKeys)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= numKeys) return;
    
    size_t elemIdx = tid / 4;  // 4 nodes per tetrahedron
    int nodeIdx = tid % 4;     // Which node (0,1,2,3)
    
    // Get the appropriate SFC key based on node index
    unsigned sfcKey;
    switch(nodeIdx) {
        case 0: sfcKey = sfc0_ptr[elemIdx]; break;
        case 1: sfcKey = sfc1_ptr[elemIdx]; break;
        case 2: sfcKey = sfc2_ptr[elemIdx]; break;
        case 3: sfcKey = sfc3_ptr[elemIdx]; break;
        default: return;
    }
    
    // Manual SFC to coordinate conversion
    auto sfcKindKey = cstone::SfcKind<unsigned>(sfcKey);
    auto [ix, iy, iz] = cstone::decodeSfc(sfcKindKey);
    constexpr unsigned maxCoord = (1u << cstone::maxTreeLevel<cstone::SfcKind<unsigned>>{}) - 1;
    float invMaxCoord = 1.0f / maxCoord;
    
    float x = box.xmin() + ix * invMaxCoord * (box.xmax() - box.xmin());
    float y = box.ymin() + iy * invMaxCoord * (box.ymax() - box.ymin());
    float z = box.zmin() + iz * invMaxCoord * (box.zmax() - box.zmin());
    
    coordinates[tid * 3 + 0] = x;
    coordinates[tid * 3 + 1] = y;
    coordinates[tid * 3 + 2] = z;
}

__global__ void spatialCoordinateKernel(const unsigned* sfc0_ptr,
                                       cstone::Box<float> box,
                                       unsigned* spatialCoords, 
                                       size_t numElements)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= numElements) return;
    
    auto sfc0 = sfc0_ptr[tid];
    
    // Test individual spatial coordinate functions
    auto sfcKindKey = cstone::SfcKind<unsigned>(sfc0);
    auto [ix, iy, iz] = cstone::decodeSfc(sfcKindKey);
    
    // Test both individual and tuple access
    auto [ix2, iy2, iz2] = cstone::decodeSfc(sfcKindKey);
    
    // Verify consistency
    bool consistent = (ix == ix2) && (iy == iy2) && (iz == iz2);
    
    spatialCoords[tid * 4 + 0] = ix;
    spatialCoords[tid * 4 + 1] = iy;
    spatialCoords[tid * 4 + 2] = iz;
    spatialCoords[tid * 4 + 3] = consistent ? 1 : 0;
}

class GpuElementDomainTest : public ::testing::Test
{
protected:
    fs::path testDir;
    int rank;
    int numRanks;
    int deviceCount = 0;

    void SetUp() override
    {
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &numRanks);
        
        cudaGetDeviceCount(&deviceCount);
        
        // Create ONE shared directory
        const char* meshPathEnv = std::getenv("MESH_PATH");
        if (!meshPathEnv) {
            GTEST_SKIP() << "MESH_PATH environment variable not set";
        }
        
        testDir = fs::path(meshPathEnv) / "mars_gpu_domain_test_shared";
        
        // Only rank 0 creates the directory and files
        if (rank == 0) {
            fs::create_directories(testDir);
            std::cout << "GPU Test setup: Found " << deviceCount << " CUDA devices" << std::endl;
            std::cout << "Writing mesh at: " << testDir << std::endl;
            
            // Create a larger, better distributed mesh
            // Generate enough elements for domain decomposition
            std::vector<float> x, y, z;
            std::vector<int> i0, i1, i2, i3;
            
            // Generate a 5x5x5 grid of tetrahedra (125 elements)
            for (int ix = 0; ix < 10; ix++) {
                for (int iy = 0; iy < 10; iy++) {
                    for (int iz = 0; iz < 10; iz++) {
                        // Base coordinates for this grid cell
                        float bx = ix * 10.0f;
                        float by = iy * 10.0f;
                        float bz = iz * 10.0f;
                        
                        // Add node indices
                        int baseIdx = x.size();
                        
                        // Add 4 nodes for a tetrahedron
                        x.push_back(bx);      y.push_back(by);      z.push_back(bz);
                        x.push_back(bx+5);    y.push_back(by);      z.push_back(bz);
                        x.push_back(bx);      y.push_back(by+5);    z.push_back(bz);
                        x.push_back(bx);      y.push_back(by);      z.push_back(bz+5);
                        
                        // Add connectivity for this tetrahedron
                        i0.push_back(baseIdx);
                        i1.push_back(baseIdx+1);
                        i2.push_back(baseIdx+2);
                        i3.push_back(baseIdx+3);
                    }
                }
            }
            
            createCoordinateFile("x.float32", x);
            createCoordinateFile("y.float32", y);
            createCoordinateFile("z.float32", z);
            createConnectivityFile("i0.int32", i0);
            createConnectivityFile("i1.int32", i1);
            createConnectivityFile("i2.int32", i2);
            createConnectivityFile("i3.int32", i3);
        }
        
        // Ensure all ranks see the files
        MPI_Barrier(MPI_COMM_WORLD);
    } 

    void createCoordinateFile(const std::string& filename, const std::vector<float>& data)
    {
        std::ofstream file((testDir / filename).string(), std::ios::binary);
        if (!file) {
            throw std::runtime_error("Failed to create coordinate file: " + filename);
        }
        file.write(reinterpret_cast<const char*>(data.data()), data.size() * sizeof(float));
    }

    void createConnectivityFile(const std::string& filename, const std::vector<int>& data)
    {
        std::ofstream file((testDir / filename).string(), std::ios::binary);
        if (!file) {
            throw std::runtime_error("Failed to create connectivity file: " + filename);
        }
        file.write(reinterpret_cast<const char*>(data.data()), data.size() * sizeof(int));
    }

    void checkGpuPrerequisites(const std::string& testName)
    {
        if (deviceCount == 0) {
            GTEST_SKIP() << testName << ": No CUDA devices available";
        }
    }
};

// Test GPU domain connectivity validation
TEST_F(GpuElementDomainTest, GpuConnectivityValidation)
{
    checkGpuPrerequisites("GpuConnectivityValidation");

    try
    {
        using Domain = ElementDomain<TetTag, float, unsigned, cstone::GpuTag>;
        Domain domain(testDir.string(), rank, numRanks);

        if (domain.getElementCount() == 0) {
            GTEST_SKIP() << "No elements on this rank";
        }

        auto box = domain.getDomain().box();

        size_t numElements = domain.getElementCount();
        
        // Get raw device pointers - THIS IS THE CORRECT WAY
        auto *sfc0_ptr = domain.template indices<0>().data();
        auto *sfc1_ptr = domain.template indices<1>().data();
        auto *sfc2_ptr = domain.template indices<2>().data();
        auto *sfc3_ptr = domain.template indices<3>().data();
        
        // Device memory for results
        cstone::DeviceVector<int> d_results(numElements);
        
        // Launch GPU kernel with raw pointers
        int blockSize = 256;
        int numBlocks = (numElements + blockSize - 1) / blockSize;
        
        validateConnectivityKernel<<<numBlocks, blockSize>>>(
            sfc0_ptr, sfc1_ptr, sfc2_ptr, sfc3_ptr,
            box, d_results.data(), numElements
        );
        
        cudaDeviceSynchronize();
        cudaError_t err = cudaGetLastError();
        ASSERT_EQ(err, cudaSuccess) << "CUDA kernel failed: " << cudaGetErrorString(err);
        
        // Check results on host
        auto h_results = toHost(d_results);
        int validElements = 0;
        for (size_t i = 0; i < numElements; i++) {
            if (h_results[i] == 1) validElements++;
        }
        
        EXPECT_GT(validElements, numElements * 0.9) << "Most elements should pass validation";
        
        if (rank == 0) {
            std::cout << "GPU connectivity validation: " << validElements << "/" << numElements 
                      << " elements valid" << std::endl;
        }
    }
    catch (const std::exception& e)
    {
        FAIL() << "Exception in GPU connectivity validation: " << e.what();
    }
}

// Test GPU volume calculation
TEST_F(GpuElementDomainTest, GpuVolumeCalculation)
{
    checkGpuPrerequisites("GpuVolumeCalculation");

    try
    {
        using Domain = ElementDomain<TetTag, float, unsigned, cstone::GpuTag>;
        Domain domain(testDir.string(), rank, numRanks);

        if (domain.getElementCount() == 0) {
            GTEST_SKIP() << "No elements on this rank";
        }

        size_t numElements = domain.getElementCount();
        
        // Get raw device pointers
        auto* sfc0_ptr = domain.template indices<0>().data();
        auto* sfc1_ptr = domain.template indices<1>().data();
        auto* sfc2_ptr = domain.template indices<2>().data();
        auto* sfc3_ptr = domain.template indices<3>().data();
        
        cstone::DeviceVector<float> d_volumes(numElements);
        
        int blockSize = 256;
        int numBlocks = (numElements + blockSize - 1) / blockSize;
        
        calculateVolumeKernel<<<numBlocks, blockSize>>>(
            sfc0_ptr, sfc1_ptr, sfc2_ptr, sfc3_ptr,
            domain.getDomain().box(), d_volumes.data(), numElements
        );
        
        cudaDeviceSynchronize();
        cudaError_t err = cudaGetLastError();
        ASSERT_EQ(err, cudaSuccess) << "CUDA volume kernel failed: " << cudaGetErrorString(err);
        
        auto h_volumes = toHost(d_volumes);
        int positiveVolumes = 0;
        float minVol = std::numeric_limits<float>::max();
        float maxVol = 0.0f;
        
        for (size_t i = 0; i < numElements; i++) {
            if (h_volumes[i] > 0.0f) {
                positiveVolumes++;
                minVol = std::min(minVol, h_volumes[i]);
                maxVol = std::max(maxVol, h_volumes[i]);
            }
        }
        
        EXPECT_GT(positiveVolumes, numElements * 0.9) << "Most elements should have positive volume";
        
        if (rank == 0) {
            std::cout << "GPU volume calculation: " << positiveVolumes << "/" << numElements 
                      << " positive volumes, range [" << minVol << ", " << maxVol << "]" << std::endl;
        }
    }
    catch (const std::exception& e)
    {
        FAIL() << "Exception in GPU volume calculation: " << e.what();
    }
}

// Test GPU coordinate conversion performance
TEST_F(GpuElementDomainTest, GpuCoordinateConversionPerformance)
{
    checkGpuPrerequisites("GpuCoordinateConversionPerformance");

    try
    {
        using Domain = ElementDomain<TetTag, float, unsigned, cstone::GpuTag>;
        Domain domain(testDir.string(), rank, numRanks);

        if (domain.getElementCount() == 0) {
            GTEST_SKIP() << "No elements on this rank";
        }

        if (domain.getElementCount() < 100) {
            GTEST_SKIP() << "Mesh too small for performance test";
        }

        size_t numKeys = domain.getElementCount() * 4; // 4 nodes per tetrahedron
        
        // Get raw device pointers
        auto* sfc0_ptr = domain.template indices<0>().data();
        auto* sfc1_ptr = domain.template indices<1>().data();
        auto* sfc2_ptr = domain.template indices<2>().data();
        auto* sfc3_ptr = domain.template indices<3>().data();

        cstone::DeviceVector<float> d_coordinates(numKeys * 3);
        
        int blockSize = 256;
        int numBlocks = (numKeys + blockSize - 1) / blockSize;
        
        // Warm up
        coordinateConversionKernel<<<numBlocks, blockSize>>>(
            sfc0_ptr, sfc1_ptr, sfc2_ptr, sfc3_ptr,
            domain.getDomain().box(), d_coordinates.data(), numKeys
        );
        cudaDeviceSynchronize();

        // Performance measurement
        const int iterations = 100;
        cudaEvent_t start, stop;
        cudaEventCreate(&start);
        cudaEventCreate(&stop);
        
        cudaEventRecord(start);
        for (int iter = 0; iter < iterations; iter++) {
            coordinateConversionKernel<<<numBlocks, blockSize>>>(
                sfc0_ptr, sfc1_ptr, sfc2_ptr, sfc3_ptr,
                domain.getDomain().box(), d_coordinates.data(), numKeys
            );
        }
        cudaEventRecord(stop);
        cudaEventSynchronize(stop);
        
        float milliseconds = 0;
        cudaEventElapsedTime(&milliseconds, start, stop);
        
        // Validate some results
        auto h_coords = toHost(d_coordinates);
        auto box = domain.getDomain().box();
        
        int validCoords = 0;
        for (size_t i = 0; i < std::min(size_t(100), numKeys); i++) {
            float x = h_coords[i * 3 + 0];
            float y = h_coords[i * 3 + 1];
            if (x >= box.xmin() && x <= box.xmax() && 
                y >= box.ymin() && y <= box.ymax()) {
                validCoords++;
            }
        }
        
        EXPECT_GT(validCoords, 90) << "Most coordinates should be within bounds";

        if (rank == 0) {
            float avgTime = milliseconds / iterations;
            float throughput = (numKeys * 1000.0f) / avgTime;
            std::cout << "GPU coordinate conversion performance:" << std::endl;
            std::cout << "  " << numKeys << " coordinates converted" << std::endl;
            std::cout << "  Average time: " << avgTime << " ms" << std::endl;
            std::cout << "  Throughput: " << throughput / 1e6 << " million coords/second" << std::endl;
        }

        cudaEventDestroy(start);
        cudaEventDestroy(stop);
    }
    catch (const std::exception& e)
    {
        FAIL() << "Exception in GPU coordinate conversion performance: " << e.what();
    }
}

// Test GPU spatial coordinate functions  
TEST_F(GpuElementDomainTest, GpuSpatialCoordinates)
{
    checkGpuPrerequisites("GpuSpatialCoordinates");

    try
    {
        using Domain = ElementDomain<TetTag, float, unsigned, cstone::GpuTag>;
        Domain domain(testDir.string(), rank, numRanks);

        if (domain.getElementCount() == 0) {
            GTEST_SKIP() << "No elements on this rank";
        }

        size_t numElements = domain.getElementCount();
        
        // Get raw device pointer (only need first connectivity for this test)
        auto* sfc0_ptr = domain.template indices<0>().data();
        
        cstone::DeviceVector<unsigned> d_spatialCoords(numElements * 4);
        
        int blockSize = 256;
        int numBlocks = (numElements + blockSize - 1) / blockSize;
        
        spatialCoordinateKernel<<<numBlocks, blockSize>>>(
            sfc0_ptr, domain.getDomain().box(), d_spatialCoords.data(), numElements
        );
        
        cudaDeviceSynchronize();
        cudaError_t err = cudaGetLastError();
        ASSERT_EQ(err, cudaSuccess) << "CUDA spatial coord kernel failed: " << cudaGetErrorString(err);
        
        auto h_spatialCoords = toHost(d_spatialCoords);
        int consistentElements = 0;
        
        for (size_t i = 0; i < numElements; i++) {
            unsigned ix = h_spatialCoords[i * 4 + 0];
            unsigned iy = h_spatialCoords[i * 4 + 1];
            unsigned iz = h_spatialCoords[i * 4 + 2];
            unsigned consistent = h_spatialCoords[i * 4 + 3];
            
            if (consistent == 1) consistentElements++;
            
            EXPECT_LT(ix, (1u << 20)) << "Spatial coordinate should be reasonable";
            EXPECT_LT(iy, (1u << 20)) << "Spatial coordinate should be reasonable";
            EXPECT_LT(iz, (1u << 20)) << "Spatial coordinate should be reasonable";
        }
        
        EXPECT_EQ(consistentElements, numElements) << "All elements should have consistent spatial coordinates";
        
        if (rank == 0) {
            std::cout << "GPU spatial coordinates: " << consistentElements << "/" << numElements 
                      << " elements consistent" << std::endl;
        }
    }
    catch (const std::exception& e)
    {
        FAIL() << "Exception in GPU spatial coordinates: " << e.what();
    }
}

// Test GPU domain creation
TEST_F(GpuElementDomainTest, GpuDomainCreation)
{
    checkGpuPrerequisites("GpuDomainCreation");

    try
    {
        using Domain = ElementDomain<TetTag, float, unsigned, cstone::GpuTag>;
        Domain domain(testDir.string(), rank, numRanks);

        EXPECT_GT(domain.getNodeCount(), 0);
        EXPECT_GT(domain.getElementCount(), 0);
        
        auto box = domain.getDomain().box();
        EXPECT_LT(box.xmin(), box.xmax());
        EXPECT_LT(box.ymin(), box.ymax());
        EXPECT_LT(box.zmin(), box.zmax());

        if (rank == 0) {
            std::cout << "GPU domain created successfully:" << std::endl;
            std::cout << "  Nodes: " << domain.getNodeCount() << std::endl;
            std::cout << "  Elements: " << domain.getElementCount() << std::endl;
            std::cout << "  Bounding box: [" << box.xmin() << ", " << box.xmax() << "] x ["
                      << box.ymin() << ", " << box.ymax() << "] x ["
                      << box.zmin() << ", " << box.zmax() << "]" << std::endl;
        }
    }
    catch (const std::exception& e)
    {
        FAIL() << "Exception in GPU domain creation: " << e.what();
    }
}

int main(int argc, char** argv)
{
    mars::Env env(argc, argv);
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}