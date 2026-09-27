#!/usr/bin/env python3
"""Host-check production timestep factors; does not validate CUDA or boundary operators."""
from pathlib import Path
import re
import subprocess
import tempfile


def coefficients(source, function, end, output):
    section = source[source.index('void '+function+'('):]
    section = section[:section.index(end)]
    declarations = re.findall(r'const (?:bool bdf2Active\w*|RealType (?:invDt|dtEff))\s*=.*?;', section, re.S)
    if not declarations:
        raise ValueError('Missing production timestep coefficients: '+function)
    return 'double '+function+'(State s, double dt) { using RealType=double;\n' + '\n'.join(declarations) + '\nreturn '+output+'; }\n'


def main():
    source = (Path(__file__).resolve().parents[1]/'backend/distributed/unstructured/fem/mars_ns_channel_solver.hpp').read_text()
    cpp = '''#include <vector>
#include <cmath>
#include <cstdio>
struct State { bool useBdf2; int bdfStep; std::vector<double> d_valuesVel_bdf2; };
'''
    cpp += coefficients(source,'runPressureSolveStep','cstone::DeviceVector<RealType> b(', 'invDt')
    cpp += coefficients(source,'runImplicitDiffusionStep','cstone::DeviceVector<RealType> b(', 'invDt')
    cpp += coefficients(source,'runCorrectorStep','// Chorin closes', 'dtEff')
    cpp += '''int main() {
 int checks=0,failures=0;
 for (bool enabled : {false,true}) for (int step : {0,1,2})
 for (int matrix : {0,1}) for (double dt : {.01,.125}) {
  State s{enabled,step,std::vector<double>(matrix)};
  double pressure=runPressureSolveStep(s,dt), diffusion=runImplicitDiffusionStep(s,dt);
  double correction=runCorrectorStep(s,dt);
  double expected=(enabled && step>=1 && matrix)?1.5/dt:1./dt;
  for(bool pass : {std::abs(pressure/expected-1)<1e-14,
                  std::abs(pressure/diffusion-1)<1e-14,
                  std::abs(pressure*correction-1)<1e-14}) {
   ++checks; if(!pass) ++failures;
  }
 }
 std::printf("%s: %d timestep-factor checks; %d failures; boundary/CUDA/MPI not tested\\n",
             failures?"FAIL":"PASS",checks,failures);
 return failures?1:0;
}
'''
    scratch = Path(__file__).resolve().parents[1] / '.local-worktrees' / 'poiseuille-host'
    scratch.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='bdf-', dir=str(scratch)) as tmp:
        cpp_path, exe = Path(tmp)/'check.cpp', Path(tmp)/'check'
        cpp_path.write_text(cpp)
        subprocess.run(['c++','-std=c++17','-O2','-Wall','-Wextra','-Werror',str(cpp_path),'-o',str(exe)],check=True)
        return subprocess.run([str(exe)]).returncode


if __name__ == '__main__':
    raise SystemExit(main())
