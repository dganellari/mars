"""Capture six pinned Tet4/Tri3 boundary assembly blocks, before global scatter."""
import re
from openaccel_node_instrumentation import once

NAMES = ('InletSpecifiedVelocity_', 'OutletSpecifiedPressure_', 'WallNoSlip_')
PREFIX = 'src/assemble/flow/segregatedFlow/'


def function_range(text, family, suffix):
    name = 'assembleElemTermsBoundary'+suffix
    match = re.search(r'void '+family+r'Assembler::\s*'+name+r'\(', text)
    if not match: raise ValueError('missing boundary function: '+name)
    start = match.start()
    end = text.find('\nvoid '+family+'Assembler::', match.end())
    return start, end if end >= 0 else text.rfind('\n} // namespace')


def captured_fields(stage):
    # Maps are exported in local-node numbering alongside the global connected-node list.
    wall = stage == 5
    fields = [('face_nodes', 'std::vector<double>{0,1,2}' if wall else 'export_face_nodes'),
              ('nearest', 'export_nearest'), ('opposing', 'export_opposing'),
              ('reversal', 'export_flags'), ('area', 'areaVec, 9'),
              ('shape', 'ws_velocity_face_shape_function')]
    if stage in (0, 2):
        fields += [('density', 'ws_rho'), ('boundary_velocity', 'UbcVec, 9')]
    elif stage == 1:
        fields += [('density', 'ws_rho'), ('velocity', 'export_velocity'), ('pressure', 'ws_p'),
                   ('gradient', 'ws_dndx'), ('pressure_gradient', 'ws_Gpdx_elem'),
                   ('influence_lhs', 'ws_du'), ('influence_rhs', 'ws_duRhs'), ('bc_multiplier', 'ws_bcMultiplier')]
    elif stage in (3, 4):
        fields += [('velocity', 'ws_U'), ('gradient', 'ws_dndx'), ('viscosity', 'ws_muEff'), ('stored_flux', 'mDot, 3')]
        if stage == 3:
            fields += [('boundary_velocity', 'UbcVec, 9'), ('bc_multiplier', 'ws_bcMultiplier')]
    else:
        fields += [('velocity', 'ws_U'), ('boundary_velocity', 'UbcVec, 9'), ('wall_coefficient', 'uWallCoeffsBip, 3')]
    return fields


def instrument_function(text, stage):
    name = ('pressure.' if stage < 3 else 'momentum.')+('inlet', 'outlet', 'wall')[stage % 3]
    anchor = '    stk::mesh::MetaData& metaData = mesh.metaDataRef();'
    text = once(text, anchor, anchor+'\n    mars_reference::BoundaryExport export_boundary(bulkData, "'+name+'", SPATIAL_DIM);')
    wall = stage == 5
    capture = '''            if (export_boundary.active()) export_boundary.capture([&] {
                mars_reference::require(nodesPerSide == 3 && numScsBip == 3, "boundary capture requires Tri3");
'''
    if stage == 4:
        # Momentum outlet assembly does not otherwise need this parent/face map.
        capture += '                const auto* faceNodeOrdinals = meSCS->side_node_ordinals(faceOrdinal);\n'
    if not wall:
        capture += '''                mars_reference::require(nodesPerElement == 4, "boundary capture requires Tet4");
                for (int f = 0; f < 3; ++f)
                    mars_reference::require(bulkData.identifier(connectedNodes[faceNodeOrdinals[f]]) == bulkData.identifier(sideNodeRels[f]),
                                            "parent/face node ordering differs");
'''
    capture += '''                mars_reference::require(!domain->isMaterialCompressible() && !domain->zonePtr()->frameRotating()
                                        && !domain->zonePtr()->meshMoving(), "boundary capture requires fixed incompressible flow");
                std::vector<double> export_face_nodes(3), export_nearest(3), export_opposing(3, 0), export_flags(3, 0);
                for (int sample = 0; sample < 3; ++sample) {
'''
    if wall:
        capture += '                    export_face_nodes[sample] = sample;\n                    export_nearest[sample] = faceIpNodeMap[sample];\n'
    else:
        capture += '                    export_face_nodes[sample] = faceNodeOrdinals[sample];\n                    export_nearest[sample] = ipNodeMap[sample];\n                    export_opposing[sample] = meSCS->opposingNodes(faceOrdinal, sample);\n'
    if stage in (0, 1, 4):
        capture += '                    export_flags[sample] = rfflag[sample];\n'
    capture += '                }\n'
    if stage == 1:
        capture += '''                mars_reference::require(cvpgHarm == 0, "harmonic gradient blending not supported in boundary replay");
                for (int f = 0; f < 3; ++f)
                    for (int j = 0; j < 3; ++j)
                        mars_reference::require(ws_Gpdx[3*f+j] == ws_Gpdx_elem[3*faceNodeOrdinals[f]+j],
                                                "face/parent pressure gradient differs");
                for (double value : ws_F) mars_reference::require(value == 0, "body force not supported in boundary replay");
                for (double value : ws_FOrig_elem) mars_reference::require(value == 0, "body force not supported in boundary replay");
                for (double value : ws_F_elem) mars_reference::require(value == 0, "body force not supported in boundary replay");
                std::vector<double> export_velocity(12, 0);
                for (int f = 0; f < 3; ++f)
                    for (int j = 0; j < 3; ++j) export_velocity[3*faceNodeOrdinals[f]+j] = ws_U[3*f+j];
'''
    capture += '                export_boundary.block(side, connectedNodes, '+('3' if stage >= 3 else '1')+', {\n'
    capture += ',\n'.join('                    {"'+key+'", '+value+'}' for key,value in captured_fields(stage))
    capture += '\n                }, lhs, rhs);\n            });\n\n'
    return once(text, '            this->applyCoeff_(', capture+'            this->applyCoeff_(')


def add_boundary_edits(source, edits, reference_dir):
    for family, base in (('pressureCorrection', 0), ('navierStokes', 3)):
        name = PREFIX+family+'/'+family+'AssemblerElemBoundaryConditions.cpp'
        before = (source/name).read_text()
        after = before
        for offset, suffix in enumerate(NAMES):
            a,b = function_range(after, family, suffix)
            after = after[:a]+instrument_function(after[a:b], base+offset)+after[b:]
        after = once(after, '#include "flowModel.h"', '#include "flowModel.h"\n#include "../../../../basic/mars_reference/boundary_export.hpp"')
        edits[name] = (before, after)
    edits['src/basic/mars_reference/boundary_export.hpp'] = ('', (reference_dir/'boundary_export.hpp').read_text())
