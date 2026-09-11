#!/usr/bin/env python3
"""Failure-injection checks for the public channel log acceptance tool."""
import unittest
from outlet_channel_check import compare, read_run


def log(ranks):
    lines = [f'[outlet-channel-config] ranks={ranks} steps=2 dt=0.01 nu=0.1 rho=1 '
             f'inlet=0.1 ramp_steps=4 beta=0.05 p_ref=0 empty_opening_ranks={int(ranks == 4)}']
    for step, bdf, dt_eff, flux in ((1, 1, .01, .025), (2, 2, .02/3, .05)):
        lines.append(f'[outlet-channel] step={step} ranks={ranks} bdf={bdf} dt_eff={dt_eff:.17g} '
                     f'rms=1e-9 max=1e-8 q_in={-flux} q_out={flux} residual_sum=0 '
                     'u_mean=0.01 v_mean=0 w_mean=0 p_mean=2 u_rms=0.02 p_rms=3 '
                     'trace_mean=0 history_error=0 halo_error=0')
    return '\n'.join(lines)+'\nPASS: public outlet channel integration steps=2\n'


def corner_log():
    text = log(4).replace('p_ref=0 empty_opening_ranks=1',
                          'p_ref=0 opening_area=0.0625 cut_check=1 empty_opening_ranks=2')
    for flux in (.025, .05):
        text = text.replace(f'q_in={-flux} q_out={flux}', f'q_in={-flux/16} q_out={flux/16}')
    for step, flux in ((1, .025/16), (2, .05/16)):
        text = text.replace(f'[outlet-channel] step={step}',
            f'[outlet-cut] step={step} q25={flux} q50={flux} q75={flux} identity_error=1e-15\n'
            f'[outlet-channel] step={step}')
    return text


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

    def test_additional_coverage(self):
        logs = [log(n).replace('empty_opening_ranks=1', 'empty_opening_ranks=0') for n in (1, 2, 4)]
        self.assertEqual(compare(logs, True, corner_log()), (2, True))
        for old, new in (('empty_opening_ranks=2', 'empty_opening_ranks=0'),
                         ('[outlet-cut] step=1', '[missing] step=1'),
                         ('q25=0.0015625', 'q25=0.0013020833333333333'),
                         ('identity_error=1e-15', 'identity_error=1e-4'),
                         ('opening_area=0.0625', 'opening_area=1'),
                         ('q50=0.003125', 'q50=nan')):
            with self.subTest(old=old), self.assertRaises(ValueError):
                compare(logs, True, corner_log().replace(old, new))
        self.assertEqual(read_run(corner_log())[0]['opening_area'], .0625)
        with self.assertRaises(ValueError):
            compare([log(1), log(2), corner_log()])

    def test_empty_coverage_is_explicit(self):
        logs = [log(n).replace('empty_opening_ranks=1', 'empty_opening_ranks=0') for n in (1, 2, 4)]
        self.assertEqual(compare(logs), (2, False))
        with self.assertRaises(ValueError):
            compare(logs, require_empty=True)


if __name__ == '__main__':
    unittest.main()
