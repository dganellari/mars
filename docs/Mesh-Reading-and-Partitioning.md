# Mesh Reading and Partitioning

This module handles the input and initial processing of unstructured meshes, including file reading, format support, and domain decomposition for parallel processing.

## Purpose

The mesh reading system provides:

- Support for multiple mesh file formats
- Automatic mesh partitioning for parallel processing
- Load balancing through space-filling curves
- Memory-efficient mesh storage

## Supported Formats

- **VTK**: Visualization Toolkit format (.vtk)
- **ExodusII**: Sandia National Labs format (.exo)
- **Custom Binary**: MARS-specific optimized format
- **Partitioned Files**: Multi-file mesh representations

## Key Components

### MeshReader Class
```cpp
class MeshReader {
public:
    static std::unique_ptr<MeshData> read_mesh(const std::string& path);
    static std::vector<std::unique_ptr<MeshData>> read_partitioned_mesh(
        const std::string& base_path, int rank, int num_ranks);
};
```

### Partitioning Strategy
- Uses Hilbert space-filling curves for load balancing
- Minimizes inter-partition communication
- Maintains element locality

## Usage Example

```cpp
// Read single mesh file
auto mesh_data = MeshReader::read_mesh("mesh.vtk");

// Create partitioned domain
ElementDomain<TetTag, float, unsigned> domain("mesh_parts", rank, num_ranks);

// The domain automatically handles partitioning during construction
std::cout << "Local elements: " << domain.get_elements().size() << std::endl;
std::cout << "Total elements: " << domain.get_global_element_count() << std::endl;
```

## Partitioning Algorithm

1. **Load Global Mesh**: Read complete mesh on root rank
2. **Compute SFC**: Assign Hilbert curve indices to elements
3. **Distribute Elements**: Partition based on SFC ranges
4. **Build Local Domains**: Create rank-local element collections
5. **Exchange Halos**: Identify and exchange ghost elements

## Performance Considerations

- **Memory Scaling**: O(N) memory usage where N is local elements
- **I/O Optimization**: Parallel file reading for large meshes
- **Load Balancing**: SFC-based partitioning minimizes communication
- **Caching**: Mesh data cached to avoid re-reading

## Error Handling

- Validates mesh integrity during reading
- Checks for supported element types
- Reports partitioning imbalances
- Handles file I/O errors gracefully

## Integration with ElementDomain

The mesh reading system integrates seamlessly with `ElementDomain`:

```cpp
// Domain construction automatically triggers mesh reading
ElementDomain<TetTag, float, unsigned> domain(mesh_path, rank, num_ranks);

// Mesh data is loaded lazily on first access
auto elements = domain.get_elements();  // Triggers reading if not done
```

## See Also

- [ElementDomain Overview](ElementDomain-Overview.md)
- [SFC Mapping](SFC-Mapping.md)
- [Multi-Rank Support](Multi-Rank-Support.md)