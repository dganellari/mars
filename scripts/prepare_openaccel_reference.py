#!/usr/bin/env python3
"""prepare a pinned, opt-in reference patch without editing OpenAccel."""
import argparse
import difflib
import hashlib
import json
from pathlib import Path
import subprocess
import shutil

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "tests/data/public_openaccel_reference/contract_v1.json"
ASSEM = "src/assemble/flow/segregatedFlow/"


def git(source, *args):
    return subprocess.check_output(["git", "-C", str(source), *args], universal_newlines=True).strip()


def replace_once(text, old, new):
    if text.count(old) != 1:
        raise ValueError(f"patch anchor must occur exactly once: {old[:100]!r}")
    return text.replace(old, new, 1)


def frozen_inputs(stage):
    fields = [
        ('coordinates', 'ws_coordinates'), ('velocity', 'ws_U'), ('density', 'ws_rho'),
        ('velocity_shape', 'ws_velocity_shape_function'),
        ('coordinate_shape', 'ws_coordinate_shape_function'), ('shape_gradient', 'ws_dndx'),
        ('compressible', 'std::vector<double>{comp}'),
        ('velocity_shifted', 'std::vector<double>{isUShifted ? 1.0 : 0.0}'),
    ]
    if stage == 'momentum':
        fields += [
            ('viscosity', 'ws_muEff'), ('velocity_blend', 'ws_beta'),
            ('velocity_gradient', 'ws_gradU'), ('stored_flux', 'mDot, std::size_t(numScsIp)'),
            ('gradient_shifted', 'std::vector<double>{isUGradientShifted ? 1.0 : 0.0}'),
            ('nso', 'std::vector<double>{nsoFac}'),
            ('nso_fourth_factor', 'std::vector<double>{fourthFac}'),
        ]
    else:
        fields += [
            ('pressure', 'ws_p'), ('pressure_gradient', 'ws_Gpdx'),
            ('influence_lhs', 'ws_du'), ('influence_rhs', 'ws_duRhs'),
            ('density_blend', 'ws_betaRho'), ('density_gradient', 'ws_gradRho'),
            ('force', 'ws_F'), ('original_force', 'ws_FOrig'), ('mesh_velocity', 'ws_Um'),
            ('gradient_shifted', 'std::vector<double>{isPGradientShifted ? 1.0 : 0.0}'),
            ('consistent', 'std::vector<double>{consistent ? 1.0 : 0.0}'),
            ('harmonic_gradient', 'std::vector<double>{cvpgHarm}'),
            ('mesh_moving', 'std::vector<double>{meshMoving ? 1.0 : 0.0}'),
            ('frame_rotating', 'std::vector<double>{domain->zonePtr()->frameRotating() ? 1.0 : 0.0}'),
        ]
    lines = ',\n'.join('                    {"' + name + '", ' + value + '}' for name, value in fields)
    return ('            if (export_scope.active())\n'
            '                export_scope.inputs(elementBucket[iElement], connectedNodes, lrscv, numScsIp, {\n'
            + lines + '\n                });\n')


def instrument(original, stage):
    text = replace_once(original, '#include "flowModel.h"',
                        '#include "flowModel.h"\n'
                        '#include "../../../../basic/mars_reference/stk_export_adapter.hpp"')
    anchor = '    stk::mesh::MetaData& metaData = mesh.metaDataRef();'
    text = replace_once(text, anchor, anchor + '\n'
                        f'    mars_reference::ExportScope export_scope(bulkData, "{stage}.interior", SPATIAL_DIM);')
    anchor = '            this->applyCoeff_(\n'
    text = replace_once(text, anchor,
                        frozen_inputs(stage) +
                        '            export_scope.block(elementBucket[iElement], connectedNodes, '
                        + ('SPATIAL_DIM' if stage == 'momentum' else '1')
                        + ', lhs, rhs);\n\n' + anchor)
    if stage == "momentum":
        anchor = '                const scalar tmDot = mDot[ip];'
        insert = ('\n                export_scope.sample(elementBucket[iElement], nodeRels[il], nodeRels[ir],\n'
                  '                                    tmDot, &p_scs_areav[ip * SPATIAL_DIM]);')
    else:
        anchor = '                // residual; left and right\n'
        insert = ('                export_scope.sample(elementBucket[iElement], nodeRels[il], nodeRels[ir],\n'
                  '                                    mDot, &p_scs_areav[ip * SPATIAL_DIM]);\n\n')
        return replace_once(text, anchor, insert + anchor)
    return replace_once(text, anchor, anchor + insert)


def dependency_check(original):
    anchor = 'find_package(Trilinos NO_MODULE REQUIRED)'
    addition = '''

# Finding the solver packages alone does not establish an STK-enabled install.
foreach(_mars_header stk_io/FillMesh.hpp stk_io/StkMeshIoBroker.hpp stk_mesh/base/BulkData.hpp)
  unset(_mars_stk_header_dir CACHE)
  find_path(_mars_stk_header_dir NAMES ${_mars_header}
    PATHS ${Trilinos_INCLUDE_DIRS} NO_DEFAULT_PATH)
  if(NOT _mars_stk_header_dir)
    message(FATAL_ERROR
      "OpenAccel reference: ${_mars_header} is absent from Trilinos_INCLUDE_DIRS. "
      "Use a compatible Trilinos with STK IO and Sierra migration enabled, or fix "
      "the include paths for that same installation. Trilinos_DIR=${Trilinos_DIR}; "
      "packages=${Trilinos_PACKAGE_LIST}; includes=${Trilinos_INCLUDE_DIRS}")
  endif()
endforeach()
unset(_mars_stk_header_dir CACHE)
'''
    return replace_once(original, anchor, anchor + addition)


def prepare(source, output, include_nodes=False):
    contract = json.loads(MANIFEST.read_text())
    ref = contract["reference"]
    if git(source, "rev-parse", "HEAD") != ref["revision"]:
        raise ValueError("OpenAccel HEAD is not the contract pin")
    if git(source, "status", "--porcelain", "--untracked-files=no", "--ignore-submodules=untracked"):
        raise ValueError("tracked reference source or a submodule is dirty; use a clean pinned checkout")
    if git(source / "src/solver", "rev-parse", "HEAD") != ref["solver_gitlink"]["revision"]:
        raise ValueError("LibLinSolve HEAD is not the contract pin")
    hashes = {}
    for prefix, entries in [("", ref["inspected_source_sha256"]),
                            ("src/solver/", ref["solver_gitlink"]["inspected_source_sha256"])]:
        for name, expected in entries.items():
            digest = hashlib.sha256((source / (prefix+name)).read_bytes()).hexdigest()
            if digest != expected:
                raise ValueError(f"source differs from inspected pin: {prefix+name}")
            hashes[prefix+name] = digest
    # Do not absorb a concurrent CMake edit into the generated patch.
    cmake = (source / "CMakeLists.txt").read_text()
    if cmake.rstrip() != git(source, "show", ref["revision"]+":CMakeLists.txt"):
        raise ValueError("CMakeLists.txt differs from pinned source")
    hashes["CMakeLists.txt"] = hashlib.sha256(cmake.encode()).hexdigest()
    edits = {"CMakeLists.txt": (cmake, dependency_check(cmake))}
    for stage, directory in [("momentum", "navierStokes"), ("pressure", "pressureCorrection")]:
        name = ASSEM + directory + "/" + directory + "AssemblerElemTerms.cpp"
        original = (source / name).read_text()
        edits[name] = (original, instrument(original, stage))
    for header in ["export_writer.hpp", "stk_export_adapter.hpp"]:
        edits["src/basic/mars_reference/"+header] = (
            "", (ROOT / "tests/reference/openaccel" / header).read_text())
    if include_nodes:
        from openaccel_node_instrumentation import add_edits
        add_edits(source, edits, ROOT / "tests/reference/openaccel")
    patch = "".join("".join(difflib.unified_diff(
        before.splitlines(keepends=True), after.splitlines(keepends=True),
        fromfile="a/"+name if before else "/dev/null", tofile="b/"+name))
        for name, (before, after) in edits.items())
    # A dry run against the actual tree checks paths and context; no files change.
    subprocess.run(["git", "-C", str(source), "apply", "--check", "-"],
                   input=patch, universal_newlines=True, check=True)
    output.mkdir(parents=True, exist_ok=False)
    (output / "instrumentation.patch").write_text(patch)
    provenance = {
        "schema": 1, "reference_revision": ref["revision"],
        "solver_revision": ref["solver_gitlink"]["revision"],
        "source_sha256": hashes, "patch_sha256": hashlib.sha256(patch.encode()).hexdigest(),
        "coverage": "local-interior-and-steady-nodes" if include_nodes else "local-interior-only",
        "frozen_inputs_schema": 2, "executed_reference": False,
        "node_capture": include_nodes,
        "manifest_sha256": hashlib.sha256(MANIFEST.read_bytes()).hexdigest(),
        "pending": ["STK compilation", "effective runtime controls",
                    "node/boundary blocks and coefficient parity",
                    "boundary samples", "full iteration and MPI execution"],
    }
    (output / "provenance.json").write_text(json.dumps(provenance, indent=2)+"\n")
    shutil.copyfile(ROOT / "tests/reference/openaccel/build_reference.py", output / "build_reference.py")
    return output / "instrumentation.patch"


def package_case(case, output, include_nodes=False):
    manifest = json.loads((case / "manifest.json").read_text())
    if manifest["fixture"] != "public_channel" or manifest["iterations"] != 2:
        raise ValueError("capture requires the public two-iteration channel bundle")
    for name in ("input.i", "channel.exo"):
        if hashlib.sha256((case / name).read_bytes()).hexdigest() != manifest["sha256"][name]:
            raise ValueError("public case checksum mismatch: " + name)
        shutil.copyfile(case / name, output / name)
    for name in ("run_openaccel_public.py", "openaccel_reference_check.py"):
        shutil.copyfile(ROOT / "scripts" / name, output / name)
    shutil.copyfile(MANIFEST, output / "contract_v1.json")
    manifest["purpose"] = "actual frozen interior" + (" and steady node" if include_nodes else "") + " capture; full contract and parity pending"
    manifest["require_frozen_inputs"] = True
    manifest.pop("author", None)
    if include_nodes:
        manifest["require_node_capture"] = True
        shutil.copyfile(ROOT / "scripts/openaccel_node_check.py", output / "openaccel_node_check.py")
    manifest["sha256"] = {name: hashlib.sha256((output / name).read_bytes()).hexdigest()
                          for name in ("input.i", "channel.exo", "run_openaccel_public.py",
                                       "openaccel_reference_check.py", "contract_v1.json", "provenance.json")}
    if include_nodes:
        manifest["sha256"]["openaccel_node_check.py"] = hashlib.sha256((output / "openaccel_node_check.py").read_bytes()).hexdigest()
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2)+"\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True, help="pinned public OpenAccel checkout")
    parser.add_argument("--output", type=Path, required=True, help="new bundle directory")
    parser.add_argument("--public-case", type=Path, help="verified v4 public smoke bundle to copy")
    parser.add_argument("--include-nodes", action="store_true", help="also capture steady nodes, relaxation and influence coefficients")
    args = parser.parse_args()
    try:
        print(prepare(args.source.resolve(), args.output.resolve(), args.include_nodes))
        if args.public_case:
            package_case(args.public_case.resolve(), args.output.resolve(), args.include_nodes)
    except (ValueError, OSError, subprocess.CalledProcessError) as error:
        parser.exit(1, f"ERROR: {error}\n")


if __name__ == "__main__":
    main()
