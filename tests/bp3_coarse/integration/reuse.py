"""Resume the reviewed native build without recompiling or modifying its evidence."""
import json
from pathlib import Path
import shutil
from build import mpi_args
from common import CONFIG_SHA, HYPRE, HYPRE_SHA, command, digest, save, verify_bundle, verify_files

OLD_BUNDLE = Path('/capstor/scratch/cscs/gandanie/bp3-coarse-integration-kit')
OLD_MANIFEST_SHA = 'a5f63ebde71b8b34b13aa8d63b864c86bfd7e0ff14e463ab42b01debd9049bfc'
BUILD_SHA = 'b9c150d8739155c7b143e75b2402d4c0e5e85b40b990672b2276143ea3802b3d'
MODULE = 'source/cmake/MarsBP3Coarse.cmake'
REGISTRATION = 'source/tests/bp3_coarse/register.cmake'
OLD_TEST = '''    add_test(NAME marsBP3CoarseSingle COMMAND mars_bp3_coarse_test balanced 511)
    set_tests_properties(marsBP3CoarseSingle PROPERTIES LABELS "bp3-coarse"
        ENVIRONMENT "BP3_CONSOLIDATE_ROWS=100;BP3_FIXTURE_PROFILE=mars;MPICH_GPU_SUPPORT_ENABLED=1")'''


def check_sources(old, bundle):
    before, after = verify_bundle(old), verify_bundle(bundle)
    if digest(old / 'manifest.json') != OLD_MANIFEST_SHA:
        raise ValueError('The original source snapshot changed')
    old_files = {k: v for k, v in before['files'].items() if not k.startswith('kit/')}
    new_files = {k: v for k, v in after['files'].items() if not k.startswith('kit/')}
    changed = {k for k in set(old_files) | set(new_files) if old_files.get(k) != new_files.get(k)}
    if changed != {MODULE, REGISTRATION}:
        raise ValueError('Reusing binaries requires identical numerical and dependency sources')
    text = (old / MODULE).read_text()
    if text.count(OLD_TEST) != 1 or (bundle / MODULE).read_text() != text.replace(
            OLD_TEST, '    include(${CMAKE_SOURCE_DIR}/tests/bp3_coarse/register.cmake)'):
        raise ValueError('Only CTest registration may change in the build module')


def configure_test(bundle, root, fixture, flags):
    project = root / 'ctest-project'
    project.mkdir()
    (project / 'CMakeLists.txt').write_text(
        'cmake_minimum_required(VERSION 3.24)\nproject(BP3CoarseReuse NONE)\nenable_testing()\n'
        'add_executable(mars_bp3_coarse_test IMPORTED)\n'
        'set_target_properties(mars_bp3_coarse_test PROPERTIES IMPORTED_LOCATION "' + str(fixture) + '")\n'
        'include("' + str(bundle / REGISTRATION) + '")\n')
    directory = root / 'ctest-build'
    command(['cmake', '-S', project, '-B', directory] + flags, root / 'ctest-configure.log', timeout=60)
    return directory


def resume(bundle, previous, root):
    check_sources(OLD_BUNDLE, bundle)
    receipt = previous / 'build/build.json'
    if digest(receipt) != BUILD_SHA:
        raise ValueError('This retry requires the independently checked 5266939 build')
    build = json.loads(receipt.read_text())
    if build['bundle_sha256'] != OLD_MANIFEST_SHA:
        raise ValueError('Wrong original build snapshot')
    if digest(HYPRE) != HYPRE_SHA or digest(HYPRE.parent.parent / 'include/HYPRE_config.h') != CONFIG_SHA:
        raise ValueError('The runtime Hypre library or headers changed')
    verify_files(previous / 'build', build['binaries'])
    out = root / 'build'
    out.mkdir()
    for name in build['binaries']:
        target = out / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(str(previous / 'build' / name), str(target))
    # Keep the original build receipt and logs with the copied, identical executables.
    evidence = {}
    for source in sorted((previous / 'build').iterdir()):
        if source.is_file() and source.suffix in ('.json', '.log', '.txt'):
            shutil.copy2(str(source), str(out / source.name))
            evidence[source.name] = digest(source)
    verify_files(out, build['binaries'])
    save(root / 'reused-build.json', dict(previous=str(previous), build_sha256=BUILD_SHA,
         original_manifest_sha256=OLD_MANIFEST_SHA, manifest_sha256=digest(bundle / 'manifest.json'),
         copied_evidence=evidence, binaries=build['binaries'], recompiled=False,
         changed_build_scope='CTest registration only; native binaries copied byte for byte'))
    fixture = out / 'on/examples/distributed/unstructured/mars_bp3_coarse_test'
    return configure_test(bundle, root, fixture, mpi_args())
