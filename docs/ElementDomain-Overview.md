# ElementDomain Overview

The `ElementDomain` class is the central component of MARS's unstructured mesh system, providing a unified interface for managing mesh elements, their properties, and relationships.

## Purpose

`ElementDomain` serves as the main container and manager for unstructured mesh data, handling:

- Element storage and access
- Coordinate management
- Partitioning information
- Lazy initialization of dependent components
- Type-safe element operations

## Template Parameters

```cpp
template<typename ElemTag, typename Real, typename Index>
class ElementDomain;
```

- **ElemTag**: Element type tag (e.g., `TetTag` for tetrahedrons)
- **Real**: Floating-point type for coordinates (e.g., `float`, `double`)
- **Index**: Integer type for indices (e.g., `unsigned`, `size_t`)

## Key Features

### Element Management
- Stores mesh elements with their connectivity
- Provides random access to elements by index
- Supports element iteration and filtering

### Coordinate Handling
- Manages vertex coordinates in 3D space
- Provides coordinate access and transformation
- Supports coordinate caching for performance

### Partitioning Support
- Handles domain decomposition for parallel processing
- Manages local vs. global element indices
- Supports multi-rank MPI operations

## Usage Example

```cpp
// Create domain for tetrahedral mesh
ElementDomain<TetTag, float, unsigned> domain("mesh_directory", rank, numRanks);

// Access elements
auto elements = domain.get_elements();
std::cout << "Number of elements: " << elements.size() << std::endl;

// Access coordinates
auto coords = domain.get_coordinates();
std::cout << "Number of vertices: " << coords.size() / 3 << std::endl;

// Get element at index
auto elem = domain.get_element(0);
std::cout << "Element nodes: ";
for (auto node : elem.nodes) {
    std::cout << node << " ";
}
```

## Lazy Composition

The `ElementDomain` uses lazy initialization for performance:

- Components are created only when first accessed
- Reduces memory usage and startup time
- Thread-safe initialization with mutex protection

## Dependencies

- **Mesh Reading**: For initial mesh loading
- **SFC Mapping**: For load balancing
- **Adjacency Structures**: For neighbor relationships
- **Halo Management**: For parallel ghost cells

## Implementation Notes

- Uses shared pointers for component ownership
- Implements copy-on-write semantics where appropriate
- Provides const-correct access methods
- Supports both host and device (GPU) operations

## See Also

- [Mesh Reading and Partitioning](Mesh-Reading-and-Partitioning.md)
- [SFC Mapping](SFC-Mapping.md)
- [Adjacency Structures](Adjacency-Structures.md)