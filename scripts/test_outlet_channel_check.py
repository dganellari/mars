#!/usr/bin/env python3
"""Failure-injection checks for the public channel log acceptance tool."""
import unittest
from outlet_channel_check import compare


def log(ranks):
    lines = [f'[outlet-channel-config] ranks={ranks} steps=2 dt=0.01 nu=0.1 rho=1 '
             f'inlet=0.1 ramp_steps=4 beta=0.05 p_ref=0 empty_opening_ranks={int(ranks == 4)}']
    for step, bdf, dt_eff, flux in ((1, 1, .01, .025), (2, 2, .02/3, .05)):
        lines.append(f'[outlet-channel] step={step} ranks={ranks} bdf={bdf} dt_eff={dt_eff:.17g} '
                     f'rms=1e-9 max=1e-8 q_in={-flux} q_out={flux} residual_sum=0 '
                     'u_mean=0.01 v_mean=0 w_mean=0 p_mean=2 u_rms=0.02 p_rms=3 '
                     'trace_mean=0 history_error=0 halo_error=0')
    return '\n'.join(lines)+'\nPASS: public outlet channel integration steps=2\n'


class Acceptance(unittest.TestCase):
    def test_valid(self):
        self.assertEqual(compare([log(4), log(1), log(2)], require_empty=True), (2, True))

    def test_failures(self):
        mutations = (
            ('PASS: public outlet channel integration steps=2\n', ''),
            ('[outlet-channel] step=2', '[missing] step=2'),
            ('step=2 ranks=4 bdf=2', 'step=2 ranks=4 bdf=1'),
            ('dt_eff=0.01', 'dt_eff=0.02'),
            ('u_mean=0.01', 'u_mean=nan'),
            ('u_mean=0.01', 'u_mean=0.1'),
            ('history_error=0', 'history_error=0.01'),
            ('halo_error=0', 'halo_error=0.01'),
            ('max=1e-8', 'max=1e-4'),
            ('q_out=0.025', 'q_out=0'),
            ('trace_mean=0', 'trace_mean=1'),
            ('residual_sum=0', 'residual_sum=1e-5'),
            ('nu=0.1', 'nu=0.2'),
        )
        for before, after in mutations:
            with self.subTest(before=before), self.assertRaises(ValueError):
                compare([log(1), log(2), log(4).replace(before, after)])
        with self.assertRaises(ValueError):
            compare([log(1), log(2), log(4)+'ERROR: injected failure\n'])
        with self.assertRaises(ValueError):
            compare([log(1), log(2), log(2)])

    def test_empty_coverage_is_explicit(self):
        logs = [log(n).replace('empty_opening_ranks=1', 'empty_opening_ranks=0') for n in (1, 2, 4)]
        self.assertEqual(compare(logs), (2, False))
        with self.assertRaises(ValueError):
            compare(logs, require_empty=True)


if __name__ == '__main__':
    unittest.main()
