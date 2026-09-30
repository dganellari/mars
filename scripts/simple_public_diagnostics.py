#!/usr/bin/env python3
"""Run beside a private SIMPLE log; export fixed labels and booleans only."""

import argparse
import json
import math
import re
import sys


TAGS = ('[simple-linear]', '[HypreGMRES] rejected:', '[hypre-spmv]', '[simple-pressure-audit]')
AUDIT_FLAGS = ('finite', 'zero_row', 'nonpositive_diagonal', 'positive_offdiagonal',
               'constant_mode_detected', 'residual_within_roundoff_bound', 'roundoff_bound_exceeds_limit',
               'compensated_residual_finite', 'compensated_residual_passed')
ERROR_MESSAGES = {
    'outlet_anchor_or_moment_error_seen': (
        'all outlet faces closed: no open pressure anchor or nonfinite outlet moments',
        'all outlet faces closed: no open pressure anchor; cannot solve this prescribed-inflow case'),
    'nonfinite_diagnostics_error_seen': ('nonfinite nonlinear diagnostics',),
    'continuity_consistency_error_seen': ('assembled continuity does not match boundary mass flux',),
    'momentum_assembly_error_seen': ('momentum assembly failed',),
    'flux_update_error_seen': ('flux update failed',),
    'output_error_seen': ('cannot write metrics', 'metric output failed', 'field output failed',
                          'field manifest output failed', 'output exists; choose a fresh prefix',
                          'invalid output node identity or nonfinite field',
                          'field output does not cover every source node',
                          'field gather returned duplicate or missing source nodes'),
    'cuda_error_seen': ('cannot select a CUDA device', 'an illegal memory access was encountered',
                        'out of memory', 'unspecified launch failure', 'device-side assert triggered',
                        'invalid device ordinal', 'invalid configuration argument',
                        'no kernel image is available for execution on the device',
                        'node-field halo: CUDA pack or unpack failed',
                        'prepared Hypre: graph packing CUDA launch failed'),
}
ERROR_PREFIXES = {'hypre_wrapper_error_seen': 'prepared Hypre: ',
                  'halo_error_seen': 'node-field halo: '}
FAILURE_FLAGS = tuple(ERROR_MESSAGES) + tuple(ERROR_PREFIXES) + (
    'application_error_seen', 'unclassified_application_error_seen', 'linear_candidate_missing_seen',
    'scheduler_time_limit_seen', 'scheduler_out_of_memory_seen', 'scheduler_signal_seen', 'mpi_abort_seen')

# These are public source messages, never text copied from an unrecognized error.
SOFTWARE_ERRORS = tuple(message for messages in ERROR_MESSAGES.values() for message in messages) + (
    'std::bad_alloc', 'std::bad_array_new_length', 'native geometry failed',
    'pressure tolerance reduction failed', 'linear verdict reduction failed',
    'pressure audit selection failed', 'run timing reduction failed',
    'field count reduction failed', 'field counts failed', 'device field gather failed',
    'root native mesh read or validation failed',
) + tuple('prepared Hypre: ' + message for message in (
    'invalid or late stopping tolerances', 'stopping tolerances differ between ranks',
    'GMRES/AMG refresh failed', 'fixed graph updates require the device-map overload',
    'nonfinite matrix, inconsistent partition, or empty partition/separate K in fixed graph mode',
    'matrix preparation failed', 'vector creation failed', 'AMG creation failed',
    'AMG setup failed', 'GMRES creation failed', 'GMRES/AMG setup failed',
    'undersized RHS or solution', 'nonfinite RHS or initial guess',
    'undersized or nonfinite RHS/initial guess', 'vector update failed', 'Hypre apply failed',
    'nonfinite result, extraction error, or CUDA failure', 'Krylov residual lookup failed',
    'Krylov residual norm failed', 'Hypre SpMV selection failed', 'residual vector creation failed',
    'true residual evaluation failed', 'residual audit vector extraction failed',
    'residual audit reduction failed', 'residual audit synchronization failed',
    'residual audit RHS copy failed', 'residual audit matvec failed', 'residual audit cleanup failed',
    'fixed graph contains an empty row or a missing diagonal', 'graph packing CUDA launch failed',
    'matrix creation failed', 'nonfinite matrix or numeric packing failure',
    'matrix refresh failed or changed the ParCSR object',
)) + tuple('node-field halo: ' + message for message in (
    'cannot allocate halo buffers', 'cannot create halo events', 'MPI exchange failed',
    'CUDA pack or unpack failed',
))
SOFTWARE_ERROR_LOOKUP = {message: message for message in SOFTWARE_ERRORS}
SCHEDULER_MESSAGES = {
    'Segmentation fault': 'segmentation_fault', 'Bus error': 'bus_error',
    'Killed': 'killed', 'Terminated': 'terminated',
    'DUE TO TIME LIMIT': 'time_limit', 'Out Of Memory': 'out_of_memory',
    'oom_kill': 'out_of_memory', 'oom-kill': 'out_of_memory',
}


def scheduler_error(line):
    return re.match(r'^(?:(?:srun|slurmstepd)(?:\[[0-9]+\])?: error: '
                    r'|\[[0-9T:.\-]+\] error: \*\*\* STEP )', line)


def application_error(line):
    error = re.fullmatch(r'(?:Rank [0-9]+ )?ERROR: (.+)', line)
    return re.sub(r' \(on (?:this|another) rank\)$', '', error.group(1)) if error else None


def failure_markers(line):
    seen = set()
    message = application_error(line)
    if message is not None:
        seen.add('application_error_seen')
        for key, messages in ERROR_MESSAGES.items():
            if message in messages:
                seen.add(key)
        for key, prefix in ERROR_PREFIXES.items():
            if message.startswith(prefix):
                seen.add(key)
        if re.fullmatch(r'(?:momentum|pressure correction) at SIMPLE iteration [0-9]+: '
                        r'linear solver returned no usable candidate', message):
            seen.add('linear_candidate_missing_seen')
        if message.startswith('halo exchange lists rejected on all ranks (global: '):
            seen.add('halo_error_seen')
        if len(seen) == 1 and message not in SOFTWARE_ERROR_LOOKUP:
            seen.add('unclassified_application_error_seen')
    # Only scheduler error lines count, not paths or echoed commands mentioning a limit.
    if scheduler_error(line):
        if 'DUE TO TIME LIMIT' in line:
            seen.add('scheduler_time_limit_seen')
        if re.search(r'\b(?:oom[_-]kill|out of memory)\b', line, re.IGNORECASE):
            seen.add('scheduler_out_of_memory_seen')
        if re.search(r'\b(?:Segmentation fault|Bus error|Killed|Terminated)\b', line):
            seen.add('scheduler_signal_seen')
    if re.match(r'^MPICH ERROR \[Rank [0-9]+\]', line) and 'MPI_Abort' in line:
        seen.add('mpi_abort_seen')
    if re.match(r'^application called MPI_Abort\(', line):
        seen.add('mpi_abort_seen')
    return seen


def fields(text):
    result = {}
    for key, value in re.findall(r'\b([a-z_]+)=([^\s]+)', text):
        # A merged MPI line must not silently replace an earlier value.
        result[key] = None if key in result else value
    return result


def number(record, key):
    try:
        return float(record[key])
    except (KeyError, TypeError, ValueError, OverflowError):
        return None


def finite(record, key):
    value = number(record, key)
    return None if value is None else math.isfinite(value)


def passed(record, residual_key, limit_key):
    residual = number(record, residual_key)
    limit = number(record, limit_key)
    if residual is None or limit is None:
        return None
    return (math.isfinite(residual) and math.isfinite(limit)
            and 0 <= residual <= limit)


def enum_value(record, key, allowed):
    value = record.get(key)
    return value if value in allowed else 'unknown'


def flag(record, key):
    return {'0': False, '1': True}.get(record.get(key))


def cap_reached(record):
    match = re.fullmatch(r'([0-9]+)/([0-9]+)', record.get('iterations') or '')
    if not match:
        return None
    done, cap = map(int, match.groups())
    return done >= cap if cap > 0 else None


def error_nonzero(record):
    value = record.get('solve_error') or ''
    return int(value) != 0 if re.fullmatch(r'[0-9]+', value) else None


def consensus(records, extract, unknown=None):
    values = [extract(record) for record in records]
    return values[0] if values and all(v == values[0] for v in values) else unknown


def pressure_reference(doc, simple, hypre):
    """Compare logged original-row residuals only where the source norm is known."""
    result = dict(requested=doc is not None, resolved=False, family='unknown', solver='unknown',
                  matrix_normalized=None, diagonal_scaled=None, original_norm_comparable=False,
                  target_looser=None, mars_residual_passed=None, hypre_residual_passed=None)
    try:
        library = doc['simulation']['solver']
        settings = library['solver_control']['advanced_options']['linear_solver_settings']
        config = next(settings[key] for key in ('pressure_correction', 'segregated_flow', 'default') if key in settings)
        if 'lookup' in config:
            config = library[config['lookup']]
        family = config['family'].lower()
        result['family'] = {'hypre': 'Hypre', 'petsc': 'PETSc', 'trilinos': 'Trilinos',
                            'amgsolver': 'AMGsolver', 'gmres': 'GMRES'}.get(family, 'unknown')
        for source, target in (('normalize_matrix', 'matrix_normalized'), ('diagonal_scaling', 'diagonal_scaled')):
            value = config.get(source, False)
            result[target] = value if type(value) is bool else None
        result['resolved'] = result['family'] != 'unknown'
        if family != 'hypre':
            return result
        solver = 'gmres' if 'options' not in config else config['options']['solver']['type'].lower()
        result['solver'] = {'gmres': 'GMRES', 'flexgmres': 'FlexGMRES',
                            'boomeramg': 'BoomerAMG', 'mgr': 'MGR'}.get(solver, 'unknown')
        rtol, atol = config.get('rtol', 1e-6), config.get('atol', 1e-16)
        if isinstance(rtol, bool) or isinstance(atol, bool):
            return result
        rtol, atol = float(rtol), float(atol)
        if not (math.isfinite(rtol) and math.isfinite(atol) and rtol > 0 and atol >= 0):
            return result
        if solver not in ('gmres', 'flexgmres', 'boomeramg'):
            return result
        if result['matrix_normalized'] is not False or result['diagonal_scaled'] is not False:
            return result
        # OpenAccel forwards atol to the Krylov solvers, but not to BoomerAMG.
        if solver == 'boomeramg':
            atol = 0.0
        result['original_norm_comparable'] = True
        if not simple or any(r.get('stage') != 'pressure' for r in simple):
            return result

        def compare(record, key, looser=False):
            rhs, value = number(record, 'rhs_norm'), number(record, key)
            if rhs is None or value is None or not (math.isfinite(rhs) and rhs > 0 and math.isfinite(value) and value >= 0):
                return None
            limit = max(atol, rtol * rhs)
            if not math.isfinite(limit):
                return None
            # Logs round the norms; do not decide at a near-equal threshold.
            if abs(value-limit) <= 1e-5 * max(value, limit):
                return None
            return limit > value if looser else value < limit

        result['target_looser'] = consensus(simple, lambda r: compare(r, 'acceptance_limit', True))
        result['mars_residual_passed'] = consensus(simple, lambda r: compare(r, 'mars_absolute_residual'))
        result['hypre_residual_passed'] = consensus(hypre, lambda r: compare(r, 'absolute_residual'))
    except (KeyError, TypeError, ValueError, AttributeError, StopIteration, OverflowError):
        # Arbitrary lookup names and malformed private values must not escape.
        pass
    return result


class DiagnosticState:
    def __init__(self):
        self.records = {tag: [] for tag in TAGS}
        self.completions = []
        self.false_convergence = False
        self.failures = set()
        self.software_errors = []
        self.scheduler_messages = []
        self.solver_started = False
        self.iteration_report_seen = False

    def feed(self, line):
        line = line.strip()
        self.failures.update(failure_markers(line))
        error = SOFTWARE_ERROR_LOOKUP.get(application_error(line))
        if error is not None and error not in self.software_errors:
            self.software_errors.append(error)
        if scheduler_error(line):
            for message, label in SCHEDULER_MESSAGES.items():
                if re.search(r'\b' + re.escape(message) + r'\b', line, re.IGNORECASE):
                    if label not in self.scheduler_messages:
                        self.scheduler_messages.append(label)
        for tag in TAGS:
            if line.startswith(tag):
                record = fields(line[len(tag):])
                if record not in self.records[tag]:
                    self.records[tag].append(record)
        if re.fullmatch(r'CONVERGED iterations=[0-9]+ ranks=[0-9]+ exchange_rounds=[0-9]+', line):
            self.completions.append('converged')
        if re.fullmatch(r'NOT CONVERGED: iteration limit iterations=[0-9]+ ranks=[0-9]+ exchange_rounds=[0-9]+', line):
            self.completions.append('iteration_limit')
        if re.match(r'^false convergence [12](?:,|$)', line):
            self.false_convergence = True
        if re.fullmatch(r'SIMPLE Tet4, [0-9]+ ranks \(ElementDomain/cstone\), '
                        r'(?:upwind|high-resolution), laminar', line):
            self.solver_started = True
        if re.match(r'^\[simple\] iteration=[0-9]+ momentum=', line):
            self.iteration_report_seen = True

    def result(self, exit_text, reference_deck=None):
        return summarize_state(self, exit_text, reference_deck)


def summarize(lines, exit_text, reference_deck=None):
    state = DiagnosticState()
    for line in lines:
        state.feed(line)
    return state.result(exit_text, reference_deck)


def summarize_state(state, exit_text, reference_deck):
    records, completions, failures = state.records, state.completions, state.failures

    exit_text = exit_text.strip()
    exit_code = int(exit_text) if re.fullmatch(r'[0-9]{1,3}', exit_text) else None
    if exit_code is not None and exit_code > 255:
        exit_code = None
    status = 'incomplete'
    if exit_code is not None and exit_code not in (0, 2):
        status = 'failed'
    elif len(completions) == 1:
        if exit_code == 0 and completions[0] == 'converged':
            status = 'converged'
        elif exit_code == 2 and completions[0] == 'iteration_limit':
            status = 'iteration_limit'

    simple = records['[simple-linear]']
    hypre = records['[HypreGMRES] rejected:']
    spmv = records['[hypre-spmv]']
    if simple or hypre or failures:
        # A concatenated successful run cannot hide a rejection.
        status = 'failed'
    result = {
        'schema': 'mars-simple-public-diagnostics-v1',
        'run_status': status,
        'process_exit_code': exit_code,
        'multiple_completions': len(completions) > 1,
        'simple_diagnostic_present': bool(simple),
        'hypre_rejection_present': bool(hypre),
        'linear_stage': consensus(simple, lambda r: enum_value(r, 'stage', ('momentum', 'pressure')), 'unknown'),
        'backend': consensus(hypre, lambda r: enum_value(r, 'backend', ('GMRES', 'FlexGMRES')), 'unknown'),
        'solver_accepted': consensus(simple, lambda r: flag(r, 'solver_accepted')),
        'mars_passed': consensus(simple, lambda r: flag(r, 'mars_passed')),
        'iteration_cap_reached': consensus(hypre, cap_reached),
        'hypre_error_nonzero': consensus(hypre, error_nonzero),
        'reported_residual_finite': consensus(hypre, lambda r: finite(r, 'reported_relative')),
        'reported_residual_passed': consensus(hypre, lambda r: passed(r, 'reported_relative', 'tolerance')),
        'hypre_explicit_residual_finite': consensus(hypre, lambda r: finite(r, 'absolute_residual')),
        'hypre_explicit_residual_passed': consensus(hypre, lambda r: passed(r, 'absolute_residual', 'acceptance_limit')),
        'mars_explicit_residual_finite': consensus(simple, lambda r: finite(r, 'mars_absolute_residual')),
        'mars_explicit_residual_passed': consensus(simple, lambda r: passed(r, 'mars_absolute_residual', 'acceptance_limit')),
        'spmv_backend': consensus(spmv, lambda r: {'0': 'native', '1': 'vendor'}.get(r.get('vendor_requested'), 'unknown'), 'unknown'),
        'false_convergence_message_seen': state.false_convergence,
        'solver_started': state.solver_started,
        'iteration_report_seen': state.iteration_report_seen,
        'software_errors': list(state.software_errors),
        'first_software_error': state.software_errors[0] if state.software_errors else 'unknown',
        'scheduler_messages': list(state.scheduler_messages),
    }
    result.update((key, key in failures) for key in FAILURE_FLAGS)
    audit = records['[simple-pressure-audit]']
    result['pressure_audit_present'] = bool(audit)
    for key in AUDIT_FLAGS:
        result['pressure_audit_' + key] = consensus(audit, lambda r: flag(r, key))
    if reference_deck is not None:
        for key, value in pressure_reference(reference_deck, simple, hypre).items():
            result['reference_pressure_' + key] = value
    return result


class SafeParser(argparse.ArgumentParser):
    def error(self, message):
        # argparse normally echoes unrecognized arguments and their private values.
        self.exit(2, 'ERROR: invalid diagnostic arguments; use --help.\n')


def main():
    parser = SafeParser(description=__doc__)
    parser.add_argument('--log', required=True)
    parser.add_argument('--exit-file', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--reference-deck', help='Optional local OpenAccel deck; only fixed control comparisons are exported')
    args = parser.parse_args()
    try:
        reference = None
        if args.reference_deck:
            from prepare_simple_deck import load_deck
            with open(args.reference_deck, 'rb') as deck:
                reference = load_deck(deck.read())
            if not isinstance(reference, dict):
                raise ValueError('Invalid reference deck')
        with open(args.log, encoding='utf-8', errors='replace') as log:
            with open(args.exit_file, encoding='ascii') as exit_file:
                result = summarize(log, exit_file.read(64), reference)
        # An exclusive create prevents accidental replacement of an input or old summary.
        with open(args.output, 'x', encoding='ascii') as output:
            json.dump(result, output, indent=2, sort_keys=True, allow_nan=False)
            output.write('\n')
    except Exception:
        print('ERROR: diagnostic export failed; check files privately.', file=sys.stderr)
        return 1
    print('Public diagnostic summary written. Share only that JSON file.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
