"""Capture the pinned SIMPLE update sites without changing their arithmetic."""
from openaccel_node_instrumentation import once

FLOW = 'src/model/flow/flowModel.cpp'
SEQUENCE = 'src/equation/flow/segregatedFlowEquations.cpp'
BASE = 'src/equation/equation.h'
PRESSURE = 'src/equation/flow/pressureCorrection/pressureCorrectionEquation.h'


def replace_function(text, signature, edit, occurrence=0):
    start = -1
    for _ in range(occurrence+1):
        start = text.index(signature, start+1)
    brace = text.index('{', start)
    # The selected functions have balanced braces in their comments/strings too.
    depth, end = 1, brace+1
    while depth:
        depth += (text[end] == '{') - (text[end] == '}')
        end += 1
    return text[:start]+edit(text[start:end])+text[end:]


def emit(stage, entity, sample, inputs, outputs, guard='true'):
    return ('\n                if (mars_reference::update_active() && ('+guard+')) {\n'
            '                    mars_reference::update_record('+str(stage)+', '+entity+', '+sample+',\n'
            '                        {'+inputs+'}, {'+outputs+'});\n                }\n')


def sequence(text):
    anchor='    // predictor step: solve momentum'
    text=once(text,anchor,'    mars_reference::UpdateSession update_session(meshRef().bulkDataRef().parallel(),\n'
              '        SPATIAL_DIM == 3 && !controlsRef().isTransient() && pCorr_eq_->subIters() == 1);\n'
              '    mars_reference::update_phase(0);\n'+anchor)
    text=once(text,'            pCorr_eq_->preSolve();','            mars_reference::update_phase(1);\n            pCorr_eq_->preSolve();\n            mars_reference::update_phase(2);')
    text=once(text,'            FOREACH_DOMAIN(updateMassFlowRate);','            mars_reference::update_phase(3);\n            FOREACH_DOMAIN(updateMassFlowRate);')
    text=once(text,'            FOREACH_DOMAIN(updateFlowReversalFlag);','            mars_reference::update_phase(4);\n            FOREACH_DOMAIN(updateFlowReversalFlag);')
    text=once(text,'            FOREACH_DOMAIN_RAW({','            mars_reference::update_phase(5);\n            FOREACH_DOMAIN_RAW({')
    text=once(text,'                        for (label i = 0; i < SPATIAL_DIM; i++)','''                        double export_velocity[3]{};
                        if (mars_reference::update_active())
                            std::copy(Ub+SPATIAL_DIM*iNode, Ub+SPATIAL_DIM*iNode+3, export_velocity);
                        for (label i = 0; i < SPATIAL_DIM; i++)''')
    anchor='''                                dpCorrdxb[SPATIAL_DIM * iNode + i];
                        }'''
    text=once(text,anchor,anchor+emit(2,'meshRef().bulkDataRef().identifier(nodeBucket[iNode])','0',
        '{export_velocity,3}, {dub+SPATIAL_DIM*iNode,3}, {dpCorrdxb+SPATIAL_DIM*iNode,3}, double(consistent)',
        '{Ub+SPATIAL_DIM*iNode,3}','nodeBucket.owned()'))
    text=once(text,'            FOREACH_DOMAIN(updatePressureGradientField);','            mars_reference::update_phase(6);\n            FOREACH_DOMAIN(updatePressureGradientField);')
    text=once(text,'            // if converged .. break sub-iter loop','            mars_reference::update_phase(7);\n\n            // if converged .. break sub-iter loop')
    anchor='    // predictor step: solve momentum'
    text=once(text,anchor,'''    FOREACH_DOMAIN_RAW({
        if (mars_reference::update_active())
            mars_reference::require(!domain->isMaterialCompressible()
                && !domain->zonePtr()->meshMoving() && !domain->zonePtr()->frameRotating(),
                "update capture requires incompressible fixed-frame flow");
    });
'''+anchor)
    return text


def pressure(text):
    anchor='                    scalar newVal ='
    text=once(text,anchor,'''                    const scalar export_old = fieldVal[i * FIELD_DIM + k];
'''+anchor)
    anchor='                    fieldVal[i * FIELD_DIM + k] = newVal;'
    text=once(text,anchor,anchor+emit(0,'bulkData.identifier(entity)','0',
        'export_old, effectiveCorrection[id * BLOCKSIZE + STRIDE + k], effectiveRelaxValue', 'newVal',
        'mars_reference::update_active(2) && BLOCKSIZE == 1 && FIELD_DIM == 1 && STRIDE == 0 && CLIP == 0 && OFFSET == 0'))
    return text


def increment(text):
    anchor='                pCorrVal[i] = correction[row * BLOCKSIZE + STRIDE];'
    return once(text,anchor,anchor+emit(1,'bulkData.identifier(bucket[i])','0',
        'correction[row * BLOCKSIZE + STRIDE]','pCorrVal[i]'))


def flux(text, stage):
    old='                mDot[ip] = mDotURF * tmDot + (1.0 - mDotURF) * mDot[ip];'
    if stage==3:
        inputs='rhoHR, {p_uIp,3}, {p_duIp,3}, {p_dpdxIp,3}, {p_GpdxIp,3}, {p_FOrigIp,3}, {p_FIp,3}, {p_scs_areav+3*ip,3}, export_old, mDotURF'
        entity='elementBucket[iElement]'
    elif stage==4:
        inputs='rhoBip, {UbcVec+3*ip,3}, {areaVec+3*ip,3}, export_old, mDotURF'
        entity='side'
    else:
        inputs='rhoBip, {p_uBip,3}, {p_duBip,3}, {p_dpdxBip,3}, {p_GpdxBip,3}, {p_FOrigBip,3}, {p_FBip,3}, {areaVec+3*ip,3}, export_old, mDotURF, double(rfflag[ip])'
        entity='side'
        anchor='''                    mDot[ip] = 0.0;
                    continue;'''
        # Reversed samples do not construct interpolation workspaces; canonical zeros are explicit.
        text=once(text,anchor,'''                    const double export_old = mDot[ip];
                    mDot[ip] = 0.0;
                    double export_zero[21]{};'''+emit(5,'bulkData.identifier(side)','ip',
                    '1.0, {export_zero,21}, export_old, mDotURF, 1.0','mDot[ip]',
                    'bulkData.bucket(side).owned()')+'                    continue;')
    return once(text,old,'                const double export_old = mDot[ip];\n'+old+
                emit(stage,'bulkData.identifier('+entity+')','ip',inputs,'mDot[ip]',
                     'bulkData.bucket('+entity+').owned()'))


def reversal(text):
    start=text.index('            case boundaryPhysicalType::outlet:')
    end=text.index('                        case boundaryConditionType::zeroGradient:',start)
    part=text[start:end]
    anchor='                                        // calculate net mass flow for the side'
    part=once(part,anchor,'''                                        double export_old_flux[3]{}, export_old_flags[3]{};
                                        if (mars_reference::update_active()) {
                                            mars_reference::require(numScsBip == 3 && nodesPerSide == 3, "Tri3 updates only");
                                            for (int s=0;s<3;++s) { export_old_flux[s]=mDot[s]; export_old_flags[s]=revf_val[s]; }
                                        }
'''+anchor)
    anchor='''                                                    mDot[ip] = 0.0;
                                            }
                                        }
                                    }'''
    # Match the actual final flux-zero branch, retaining all native arithmetic.
    tail='''                                                mDot[ip] = 0.0;
                                            }
                                        }
                                    }'''
    part=once(part,tail,tail[:-len('                                    }')]+'''
                                        double export_flags[3]{};
                                        if (mars_reference::update_active())
                                            for (int s=0;s<3;++s) export_flags[s]=revf_val[s];
'''+emit(6,'bulkData.identifier(side)','0',
        '{export_old_flux,3}, {export_old_flags,3}, {p_U,9}, {p_p,3}, {pbc,3}, {areaVec,9}, double(ignoreFlagUpdate)',
        '{mDot,3}, {export_flags,3}','bulkData.bucket(side).owned()')+'                                    }')
    return text[:start]+part+text[end:]


def trace(text):
    # Public fixture uses the constant/time-table case; expression/UDS paths are outside this gate.
    end=text.index('        case inputDataType::expression:')
    part=text[:end]
    anchor='                                    p_estimate += pBip * aMag;'
    part=once(part,anchor,emit(7,'bulkData.identifier(side)','ip',
        '{p_p,3}, {p_face_shape_function+offSetSF_face,3}, {areaVec+offSetAreaVec,3}, 0.0',
        'pBip*aMag, aMag')+anchor)
    # Emit zero moments for excluded samples, preserving sample coverage.
    anchor='''                                    // skip
                                }
                                else'''
    loc=part.index(anchor)
    part=part[:loc]+part[loc:].replace(anchor,'''                                    double export_zero[9]{};'''+emit(7,'bulkData.identifier(side)','ip',
        '{export_zero,9}, 1.0','0.0, 0.0')+'                                }\n                                else',1)
    anchor='                    p_estimate /= area;'
    part=once(part,anchor,'                    const double export_moment = p_estimate;\n'+anchor+
              emit(9,'1','0','export_moment, area','p_estimate'))
    anchor='''                                    pbc[ip] =
                                        pAvg + (1.0 - beta) * (p - p_estimate);'''
    part=once(part,anchor,'                                    const double export_old = pbc[ip];\n'+anchor+
        emit(8,'bulkData.identifier(side)','ip','p, pAvg, p_estimate, beta, export_old, 0.0','pbc[ip]',
             'bulkData.bucket(side).owned()'))
    anchor='''                                    // skip
                                }
                                else'''
    part=once(part,anchor,emit(8,'bulkData.identifier(side)','ip',
        '0.0, pAvg, p_estimate, beta, pbc[ip], 1.0','pbc[ip]','bulkData.bucket(side).owned()')+
        '                                }\n                                else')
    return part+text[end:]


def add_update_edits(source, edits, reference_dir):
    for name, edit, include in (
        (SEQUENCE,sequence,'../../basic/mars_reference/update_export.hpp'),
        (BASE,pressure,'../basic/mars_reference/update_export.hpp'),
        (PRESSURE,increment,'../../../basic/mars_reference/update_export.hpp'),
    ):
        before=(source/name).read_text()
        after=edit(before)
        pos=after.index('#include')
        after=after[:pos]+'#include "'+include+'"\n'+after[pos:]
        edits[name]=(before,after)
    before=(source/FLOW).read_text(); after=before
    for name,stage in (('Interior_',3),('BoundaryFieldInletSpecifiedVelocity_',4),('BoundaryFieldOutletSpecifiedPressure_',5)):
        signature='void flowModel::updateMassFlowRate'+name+'('
        after=replace_function(after,signature,lambda text:flux(text,stage),1)
    after=replace_function(after,'void flowModel::updateFlowReversalFlag_(',reversal)
    after=replace_function(after,'void flowModel::updatePressureBoundarySideFieldAverageStaticPressure_(',trace)
    pos=after.index('#include')
    after=after[:pos]+'#include "../../basic/mars_reference/update_export.hpp"\n'+after[pos:]
    edits[FLOW]=(before,after)
    edits['src/basic/mars_reference/update_export.hpp']=('',(reference_dir/'update_export.hpp').read_text())
