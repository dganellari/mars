"""Check saved OpenAccel pressure reports without exporting their numerical values."""
from fractions import Fraction
import math
import re

from simple_snapshot_compare import require


ERROR = 'reference_pressure_log'


def printed_interval(text):
    # linearSystem prints two decimals in scientific notation. Keep that rounding
    # interval; its rounded Drop column and the outer Conv. hint are not verdicts.
    match = re.fullmatch(r'(\d\.\d{2})e([+-]\d{2,3})', text)
    require(match is not None, ERROR)
    exponent = int(match.group(2))
    require(-324 <= exponent <= 308 and math.isfinite(float(text)), ERROR)
    value = Fraction(match.group(1)) * Fraction(10)**exponent
    if value == 0:
        require(text == '0.00e+00', ERROR)
        return value, value
    require(value >= Fraction(10)**exponent, ERROR)
    half_unit = Fraction(10)**(exponent - 2) / 2
    return value - half_unit, value + half_unit


def residual_assessment(start, end, rtol, atol):
    low_start, high_start = printed_interval(start)
    low_end, high_end = printed_interval(end)
    relative, absolute = Fraction.from_float(float(rtol)), Fraction.from_float(float(atol))
    if high_end <= relative * low_start:
        return 'printed_relative_target_met'
    # For scalar pressure, RMS = L2/sqrt(N). N >= 1 makes atol an upper
    # bound on the RMS absolute allowance. No mesh or unknown count is needed.
    if low_end > max(relative * high_start, absolute):
        return 'printed_target_exceeded'
    return 'inconclusive_absolute_tolerance_or_print_rounding'


def inspect_pressure_log(text, steps, rtol, atol, maxiter):
    require(type(steps) is int and steps > 0 and type(maxiter) is int and maxiter > 0
            and all(math.isfinite(x) for x in (rtol, atol)) and 0 < rtol < 1 and atol >= 0, ERROR)
    iterations, blocks, rows = [], [], []
    step, pressure, context, count = None, False, None, None
    for raw in text.splitlines():
        line = raw.strip()
        if line.startswith('Iter ='):
            match = re.fullmatch(r'Iter = ([1-9][0-9]*)', line)
            require(match is not None, ERROR)
            step = int(match.group(1))
            iterations.append(step)
            pressure, context, count = False, None, None
        elif line.startswith('Equation name:'):
            pressure = line.split(':', 1)[1].strip() == 'Pressure Correction'
            context, count = None, None
            if pressure:
                require(step is not None, ERROR)
                blocks.append(step)
        elif pressure and line.startswith('Solver context:'):
            require(context is None, ERROR)
            context = line.split(':', 1)[1].strip().lower()
            require(context == 'hypre', ERROR)
        elif pressure and line.startswith('Iterations:'):
            require(count is None, ERROR)
            match = re.fullmatch(r'Iterations:\s+([0-9]+)/([1-9][0-9]*)', line)
            require(match is not None and int(match.group(2)) == maxiter, ERROR)
            count = int(match.group(1))
        elif re.match(r'^\|\s*P-Mass\s*\|', line):
            fields = [s.strip() for s in line.split('|')]
            require(pressure and context == 'hypre' and count is not None and step is not None
                    and len(fields) == 8 and fields[0] == fields[-1] == '' and fields[1] == 'P-Mass', ERROR)
            residuals = fields[4].split('/')
            require(len(residuals) == 2, ERROR)
            rows.append(dict(iteration=step, iteration_cap_reached=count >= maxiter,
                             zero_iterations=count == 0,
                             residual_assessment=residual_assessment(*residuals, rtol, atol)))
    expected = list(range(1, steps + 1))
    require(iterations == blocks == expected and [row['iteration'] for row in rows] == expected, ERROR)
    return dict(scope='rounded_printed_rms_not_independent_residual_or_backend_status',
                all_steps_present=True,
                all_printed_relative_targets_met=all(row['residual_assessment'] == 'printed_relative_target_met'
                                                    for row in rows),
                any_printed_target_exceeded=any(row['residual_assessment'] == 'printed_target_exceeded'
                                               for row in rows),
                any_inconclusive=any(row['residual_assessment'].startswith('inconclusive') for row in rows),
                any_iteration_cap_reached=any(row['iteration_cap_reached'] for row in rows),
                backend_convergence_verified=False, residuals_recomputed=False,
                absolute_tolerance_scope='conservative_bound_using_at_least_one_row', steps=rows)
