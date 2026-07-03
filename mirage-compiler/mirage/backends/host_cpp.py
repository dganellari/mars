"""Host-C++ backend: emit a plain-C++ element apply from the IR.

This is the LOCAL correctness oracle for the whole compiler. It emits the same
Knaus Alg-2 sweep as the CUDA backend, but as ordinary host code, so it can be
compiled and diffed numerically against the hand-written reference
`mars::fem::applyHoCvfemElement` on a CPU -- no GPU needed. If the generated
host apply reproduces the reference to machine precision, the IR has captured
the operator correctly and the CUDA emission (same schedule) is trustworthy up
to FP reduction order.

For the full Laplacian flux it reproduces applyHoCvfemElement operation-for-
operation (bit-identical). For pruned fluxes (no cross terms) it drops the
tangential contractions -- demonstrating the codegen leverage.
"""

from ..ir import ElementApply, jacobian_action_flux  # noqa: F401 (ElementApply used in annotations)
from .common import NAMES_HOST, NAMES_HOST_JVP, render_flux, AUTOGEN_BANNER


def emit(ea: ElementApply) -> str:
    o = ea.op
    flux = render_flux(ea, NAMES_HOST)
    fn = f"{o.name}_apply_host"

    banner = AUTOGEN_BANNER.format(backend="host", name=o.name, flux=o.flux_src)

    # step-1: normal contraction. interp only needed when the flux has tangential
    # terms; deriv only when the flux uses it.
    step1 = []
    if ea.needs_deriv or ea.needs_tangential:
        step1.append("            for (int s = 0; s < n; ++s)")
        step1.append("                for (int r = 0; r < n; ++r) {")
        step1.append("                    double bi = 0, di = 0;")
        step1.append("                    for (int q = 0; q < n; ++q) {")
        if ea.needs_tangential:
            step1.append("                        bi += op.Btil[l * n + q] * u[idx(dir, q, s, r)];")
        if ea.needs_deriv:
            step1.append("                        di += op.Dtil[l * n + q] * u[idx(dir, q, s, r)];")
        step1.append("                    }")
        if ea.needs_tangential:
            step1.append("                    interp[s * n + r] = bi;")
        if ea.needs_deriv:
            step1.append("                    deriv[s * n + r] = di;")
        step1.append("                }")
    step1_src = "\n".join(step1)

    # step-2: tangential contractions + flux.
    s2 = []
    s2.append("            for (int s = 0; s < n; ++s)")
    s2.append("                for (int r = 0; r < n; ++r) {")
    if ea.needs_metric:
        s2.append("                    const auto& g = G[((size_t)(dir * p + l) * n + s) * n + r];")
    if ea.needs_dt2:
        s2.append("                    double dt2 = 0; for (int q = 0; q < n; ++q) dt2 += op.D[r * n + q] * interp[s * n + q];")
    if ea.needs_dt1:
        s2.append("                    double dt1 = 0; for (int q = 0; q < n; ++q) dt1 += op.D[s * n + q] * interp[q * n + r];")
    s2.append(f"                    flux[s * n + r] = {flux};")
    s2.append("                }")
    step2_src = "\n".join(s2)

    # Declarations: only allocate the buffers a step actually uses.
    decls = ["std::vector<double> flux(n * n, 0.0), tmp(n * n, 0.0), intf(n * n, 0.0);"]
    if ea.needs_tangential:
        decls.insert(0, "std::vector<double> interp(n * n, 0.0);")
    if ea.needs_deriv:
        decls.insert(0, "std::vector<double> deriv(n * n, 0.0);")
    decls_src = "\n            ".join(decls)

    return f"""#pragma once
{banner}
#include "mars_cvfem_ho_apply.hpp"
#include <vector>
#include <array>
#include <cstddef>

namespace mars {{
namespace mirage {{

// y = K_e u for one tensor-product hex, generated from operator '{o.name}'.
// Signature and metric layout G match mars::fem::applyHoCvfemElement.
inline void {fn}(const mars::fem::HoCvfemOperators& op,
                 const std::vector<std::array<double, 3>>& G,
                 const std::vector<double>& u, std::vector<double>& y)
{{
    const int p = op.p, n = p + 1, n3 = n * n * n;
    y.assign(n3, 0.0);
    auto idx = [&](int dir, int nrm, int t1, int t2) -> int {{
        if (dir == 0) return nrm * n * n + t1 * n + t2;
        if (dir == 1) return t1 * n * n + nrm * n + t2;
        return t1 * n * n + t2 * n + nrm;
    }};
    for (int dir = 0; dir < 3; ++dir)
        for (int l = 0; l < p; ++l) {{
            {decls_src}
{step1_src}
{step2_src}
            for (int s = 0; s < n; ++s)
                for (int r = 0; r < n; ++r) {{ double v = 0; for (int q = 0; q < n; ++q) v += op.W[r * n + q] * flux[s * n + q]; tmp[s * n + r] = v; }}
            for (int s = 0; s < n; ++s)
                for (int r = 0; r < n; ++r) {{ double v = 0; for (int q = 0; q < n; ++q) v += op.W[s * n + q] * tmp[q * n + r]; intf[s * n + r] = v; }}
            for (int s = 0; s < n; ++s)
                for (int r = 0; r < n; ++r) {{
                    y[idx(dir, l,     s, r)] -= intf[s * n + r];
                    y[idx(dir, l + 1, s, r)] += intf[s * n + r];
                }}
        }}
}}

}}  // namespace mirage
}}  // namespace mars
"""


def emit_jvp(ea: ElementApply) -> str:
    """Emit the Jacobian-action (JVP) host apply: y = J(u) . du, where J is the
    linearization of the operator, derived by SYMBOLIC form differentiation of
    the same authored flux (see ir.jacobian_action_flux). The kernel gathers the
    primal u (for the partials) and the direction du (for the _dot contractions).
    Same integrate/scatter skeleton as the primal apply."""
    o = ea.op
    jvp, _partials = jacobian_action_flux(o)
    fv = jvp.free_vars()
    fn = f"{o.name}_jvp_host"
    flux = jvp.render(NAMES_HOST_JVP)
    banner = AUTOGEN_BANNER.format(backend="host-jvp", name=o.name,
                                   flux=f"J*du of ({o.flux_src})")

    need_p_deriv = "deriv" in fv
    need_p_dt1, need_p_dt2 = "dt1" in fv, "dt2" in fv
    need_p_interp = need_p_dt1 or need_p_dt2
    need_d_deriv = "deriv_dot" in fv
    need_d_dt1, need_d_dt2 = "dt1_dot" in fv, "dt2_dot" in fv
    need_d_interp = need_d_dt1 or need_d_dt2
    need_metric = bool(fv & {"g0", "g1", "g2"})

    # step-1: primal (from u) and perturbation (from du) normal contractions.
    acc = []
    if need_p_interp:
        acc.append("bi")
    if need_p_deriv:
        acc.append("di")
    if need_d_interp:
        acc.append("bid")
    if need_d_deriv:
        acc.append("did")
    s1 = []
    if acc:
        s1.append("            for (int s = 0; s < n; ++s)")
        s1.append("                for (int r = 0; r < n; ++r) {")
        s1.append("                    double " + ", ".join(f"{a} = 0" for a in acc) + ";")
        s1.append("                    for (int q = 0; q < n; ++q) {")
        if need_p_interp:
            s1.append("                        bi  += op.Btil[l * n + q] * u[idx(dir, q, s, r)];")
        if need_p_deriv:
            s1.append("                        di  += op.Dtil[l * n + q] * u[idx(dir, q, s, r)];")
        if need_d_interp:
            s1.append("                        bid += op.Btil[l * n + q] * du[idx(dir, q, s, r)];")
        if need_d_deriv:
            s1.append("                        did += op.Dtil[l * n + q] * du[idx(dir, q, s, r)];")
        s1.append("                    }")
        if need_p_interp:
            s1.append("                    interp[s * n + r] = bi;")
        if need_p_deriv:
            s1.append("                    deriv[s * n + r] = di;")
        if need_d_interp:
            s1.append("                    interp_dot[s * n + r] = bid;")
        if need_d_deriv:
            s1.append("                    deriv_dot[s * n + r] = did;")
        s1.append("                }")
    step1_src = "\n".join(s1)

    # step-2: primal + perturbation tangential D-contractions, then JVP flux.
    s2 = ["            for (int s = 0; s < n; ++s)",
          "                for (int r = 0; r < n; ++r) {"]
    if need_metric:
        s2.append("                    const auto& g = G[((size_t)(dir * p + l) * n + s) * n + r];")
    if need_p_dt2:
        s2.append("                    double dt2 = 0; for (int q = 0; q < n; ++q) dt2 += op.D[r * n + q] * interp[s * n + q];")
    if need_p_dt1:
        s2.append("                    double dt1 = 0; for (int q = 0; q < n; ++q) dt1 += op.D[s * n + q] * interp[q * n + r];")
    if need_d_dt2:
        s2.append("                    double dt2_dot = 0; for (int q = 0; q < n; ++q) dt2_dot += op.D[r * n + q] * interp_dot[s * n + q];")
    if need_d_dt1:
        s2.append("                    double dt1_dot = 0; for (int q = 0; q < n; ++q) dt1_dot += op.D[s * n + q] * interp_dot[q * n + r];")
    s2.append(f"                    flux[s * n + r] = {flux};")
    s2.append("                }")
    step2_src = "\n".join(s2)

    decls = ["std::vector<double> flux(n * n, 0.0), tmp(n * n, 0.0), intf(n * n, 0.0);"]
    if need_d_deriv:
        decls.insert(0, "std::vector<double> deriv_dot(n * n, 0.0);")
    if need_d_interp:
        decls.insert(0, "std::vector<double> interp_dot(n * n, 0.0);")
    if need_p_deriv:
        decls.insert(0, "std::vector<double> deriv(n * n, 0.0);")
    if need_p_interp:
        decls.insert(0, "std::vector<double> interp(n * n, 0.0);")
    decls_src = "\n            ".join(decls)

    return f"""#pragma once
{banner}
#include "mars_cvfem_ho_apply.hpp"
#include <vector>
#include <array>
#include <cstddef>

namespace mars {{
namespace mirage {{

// y = J(u) . du, Jacobian-action of operator '{o.name}', by symbolic form
// differentiation of its flux. G/idx layout match mars::fem::applyHoCvfemElement.
inline void {fn}(const mars::fem::HoCvfemOperators& op,
                 const std::vector<std::array<double, 3>>& G,
                 const std::vector<double>& u, const std::vector<double>& du,
                 std::vector<double>& y)
{{
    const int p = op.p, n = p + 1, n3 = n * n * n;
    y.assign(n3, 0.0);
    auto idx = [&](int dir, int nrm, int t1, int t2) -> int {{
        if (dir == 0) return nrm * n * n + t1 * n + t2;
        if (dir == 1) return t1 * n * n + nrm * n + t2;
        return t1 * n * n + t2 * n + nrm;
    }};
    for (int dir = 0; dir < 3; ++dir)
        for (int l = 0; l < p; ++l) {{
            {decls_src}
{step1_src}
{step2_src}
            for (int s = 0; s < n; ++s)
                for (int r = 0; r < n; ++r) {{ double v = 0; for (int q = 0; q < n; ++q) v += op.W[r * n + q] * flux[s * n + q]; tmp[s * n + r] = v; }}
            for (int s = 0; s < n; ++s)
                for (int r = 0; r < n; ++r) {{ double v = 0; for (int q = 0; q < n; ++q) v += op.W[s * n + q] * tmp[q * n + r]; intf[s * n + r] = v; }}
            for (int s = 0; s < n; ++s)
                for (int r = 0; r < n; ++r) {{
                    y[idx(dir, l,     s, r)] -= intf[s * n + r];
                    y[idx(dir, l + 1, s, r)] += intf[s * n + r];
                }}
        }}
}}

}}  // namespace mirage
}}  // namespace mars
"""
