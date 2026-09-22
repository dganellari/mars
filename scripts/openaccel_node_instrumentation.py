"""Export values at the pinned steady node and post-assembly sites."""


def once(text, old, new):
    if text.count(old) != 1:
        raise ValueError('node instrumentation anchor is not unique: '+old[:80])
    return text.replace(old, new, 1)


def section(text, start, end, edit):
    a, b = text.index(start), text.index(end, text.index(start))
    return text[:a]+edit(text[a:b])+text[b:]


def scope(text, stage, bulk='bulkData', enabled='SPATIAL_DIM == 3'):
    anchor = '    const auto* graph = A.getGraph();' if stage == 'momentum.relaxation' else (
        '    const zone* zonePtr = domain->zonePtr();' if stage == 'momentum.boundary_relaxation' else
        '    stk::mesh::BulkData& bulkData = mesh.bulkDataRef();')
    return once(text, anchor, anchor+'\n    mars_reference::NodeExport node_export('+bulk+'.parallel(), "'+stage+'", '+enabled+');')


def steady(text):
    text = scope(text, 'momentum.node')
    anchor = '            assembler::applyCoeff_('
    capture = '''            if (node_export.active())
                node_export.node(bulkData.identifier(node), {
                    {"density", std::vector<double>{rho}}, {"volume", std::vector<double>{vol}},
                    {"pseudo_dt", std::vector<double>{dt}}, {"mass_divergence", std::vector<double>{div}},
                    {"velocity", Ub + SPATIAL_DIM*iNode, 3},
                    {"pressure_gradient", gradPb + SPATIAL_DIM*iNode, 3},
                    {"force", Fb + SPATIAL_DIM*iNode, 3}, {"source", p_msrc, 3},
                    {"coriolis", p_mat, 9}
                }, {{"lhs", lhs}, {"rhs", rhs}});

'''
    return once(text, anchor, capture+anchor)


def relaxation(text):
    text = scope(text, 'momentum.relaxation', enabled='N == 3')
    anchor = '            for (label k = 0; k < BLOCKSIZE; k++)'
    text = once(text, anchor, '''            std::vector<double> export_before;
            if (node_export.active()) export_before.assign(diag, diag+9);
'''+anchor)
    anchor = '                diag[BLOCKSIZE * k + k] *= urf_inv;\n            }'
    return once(text, anchor, anchor+'''
            if (node_export.active())
                node_export.node(bulkData.identifier(entity), {
                    {"lhs", export_before}, {"alpha", std::vector<double>{urf}}
                }, {{"lhs", diag, 9}});''')


def influence(text):
    text = scope(text, 'momentum.influence')
    anchor = '''                    duTilde[i] = vol / (di + sumOffDiagi + SMALL);
                }
            }'''
    return once(text, anchor, anchor+'''
            if (node_export.active()) {
                std::vector<double> export_tilde(3, 0.0);
                if (consistent) {
                    const double* ptr = stk::mesh::field_data(*duTildeSTKFieldPtr, node);
                    export_tilde.assign(ptr, ptr+3);
                }
                node_export.node(bulkData.identifier(node), {
                    {"volume", std::vector<double>{vol}},
                    {"row_blocks", &rowVals[0], rowVals.size()},
                    {"diagonal_block", std::vector<double>{double(A.diagOffsetRef()[localID])}},
                    {"consistent", std::vector<double>{consistent ? 1.0 : 0.0}},
                    {"fractional_step", std::vector<double>{fsmflag}},
                    {"transient", std::vector<double>{is_transient ? 1.0 : 0.0}},
                    {"small", std::vector<double>{SMALL}}
                }, {{"d", dub + SPATIAL_DIM*iNode, 3}, {"d_tilde", export_tilde}});
            }''')


def boundary_relaxation(text):
    text = scope(text, 'momentum.boundary_relaxation')
    anchor = '                for (int k = 0; k < BLOCKSIZE; k++)'
    text = once(text, anchor, '''                std::vector<double> export_before;
                if (node_export.active()) export_before.assign(rhs_val, rhs_val+3);
'''+anchor)
    anchor = '                    rhs_val[k] *= urf;\n                }'
    return once(text, anchor, anchor+'''
                if (node_export.active())
                    node_export.node(bulkData.identifier(entity), {
                        {"rhs", export_before}, {"factor", std::vector<double>{urf}}
                    }, {{"rhs", rhs_val, 3}});''')


def add_edits(source, edits, reference_dir):
    prefix = 'src/assemble/flow/segregatedFlow/navierStokes/'
    name = prefix+'navierStokesAssemblerNodeTerms.cpp'
    before = (source/name).read_text()
    after = section(before, 'void navierStokesAssembler::assembleNodeTermsFusedSteady_',
                    'void navierStokesAssembler::assembleNodeTermsFusedFirstOrderUnsteady_', steady)
    after = once(after, '#include "flowModel.h"', '#include "flowModel.h"\n#include "../../../../basic/mars_reference/node_export.hpp"')
    edits[name] = (before, after)
    name = prefix+'navierStokesAssembler.cpp'
    before = (source/name).read_text()
    after = section(before, 'void navierStokesAssembler::computeDUCoefficients',
                    'void navierStokesAssembler::postAssemble', influence)
    a = after.index('void navierStokesAssembler::assembleBoundaryRelaxation_')
    b = after.index('\nvoid navierStokesAssembler::', a+1)
    after = after[:a]+boundary_relaxation(after[a:b])+after[b:]
    after = once(after, '#include "flowModel.h"', '#include "flowModel.h"\n#include "../../../../basic/mars_reference/node_export.hpp"')
    edits[name] = (before, after)
    name = 'src/assemble/phiAssembler/phiAssembler.h'
    before = (source/name).read_text()
    after = section(before, 'void phiAssembler<N>::assembleRelaxation_',
                    'void phiAssembler<N>::assembleBoundaryRelaxation_', relaxation)
    after = once(after, '#include "types.h"', '#include "types.h"\n#include "../../basic/mars_reference/node_export.hpp"')
    edits[name] = (before, after)
    edits['src/basic/mars_reference/node_export.hpp'] = ('', (reference_dir/'node_export.hpp').read_text())
