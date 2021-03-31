#include "mars_base.hpp"
#include "mars_boundary_conditions.hpp"
#include "mars_copy_operator.hpp"
#include "mars_fe_values.hpp"
#include "mars_globals.hpp"
#include "mars_gradient_recovery.hpp"
#include "mars_identity_operator.hpp"
#include "mars_interpolate.hpp"
#include "mars_invert.hpp"
#include "mars_laplace_ex.hpp"
#include "mars_precon_conjugate_grad.hpp"
#include "mars_serial_mesh_type.hpp"
#include "mars_simplex_laplacian.hpp"
#include "mars_umesh_laplace.hpp"

template <class PMesh>
class ImageWriter {
public:
    using VectorReal = mars::ViewVectorType<Real>;
    ImageWriter(VectorReal &x) {}

    bool write(VectorReal &x) { return 0; }
}