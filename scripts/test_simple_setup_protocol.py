"""Synthetic missing-message checks; no confidential cases or cluster execution."""
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import simple_setup_protocol as protocol


class ProtocolTests(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.root = Path(tmp.name)

    def fixture(self, reply=True, reply_return=True, comm=True, round_mismatch=False):
        traces = []
        for rank in range(2):
            traces.append([dict(event='setup_entry', exchange=0),
                           dict(event='exchange_entry', exchange=12345, round=2 if rank and round_mismatch else 1,
                                rank=rank, ranks=2, comm_world=comm)])

        def message(rank, kind, ident, peer, tag, returned=True):
            traces[rank].append(dict(event=kind, exchange=12345, id=ident, peer=peer, tag=tag,
                                     zero=tag == 1004, same_comm=True))
            if returned:
                traces[rank].append(dict(event='call_return', exchange=12345, id=ident, ok=True))

        message(0, 'irecv', 1, 1, 1002)
        message(0, 'isend', 2, 1, 1000)
        message(1, 'recv', 1, 0, 1000)
        if reply:
            message(1, 'send', 2, 0, 1002, reply_return)
        message(1, 'isend', 3, 0, 1004)
        message(1, 'irecv', 4, 0, 1004)
        traces[0].append(dict(event='testall', exchange=12345, requests=[1]))
        traces[1].append(dict(event='test', exchange=12345, requests=[4]))
        self.traces = traces
        self.save()

    def save(self):
        for rank, events in enumerate(self.traces):
            (self.root / ('protocol-{}.jsonl'.format(rank))).write_text(
                ''.join(json.dumps(value) + '\n' for value in events))

    def result(self):
        result = protocol.report(self.root, 2)
        text = json.dumps(result)
        for private in ('12345', '1002', '1004', str(self.root), 'SECRET'):
            self.assertNotIn(private, text)
        self.assertFalse(result['message_payloads_inspected'])
        self.assertFalse(result['protocol_cause_verified'])
        return result

    def test_missing_reply_is_separate_from_unreceived_contact(self):
        self.fixture(reply=False)
        result = self.result()
        self.assertEqual(result['assessment'], 'protocol_waits_located')
        first = result['rank_waits']['0']
        self.assertEqual(first['phase'], 'contact_replies')
        self.assertTrue(first['matching_contacts_received'])
        self.assertFalse(first['matching_send_calls_seen'])
        self.assertEqual(result['rank_waits']['1']['phase'], 'parent_release')

    def test_sent_reply_does_not_prove_delivery(self):
        self.fixture()
        result = self.result()
        self.assertTrue(result['rank_waits']['0']['matching_send_calls_returned'])
        self.assertEqual(result['assessment'], 'protocol_waits_located')

    def test_unreturned_send_is_visible(self):
        self.fixture(reply_return=False)
        result = self.result()
        self.assertTrue(result['rank_waits']['0']['matching_send_calls_seen'])
        self.assertFalse(result['rank_waits']['0']['matching_send_calls_returned'])
        self.assertTrue(result['rank_waits']['1']['message_call_in_progress'])

    def test_same_ordinal_different_round_is_not_a_matched_context(self):
        self.fixture(round_mismatch=True)
        result = self.result()
        self.assertTrue(result['same_exchange_ordinal'])
        self.assertFalse(result['same_exchange_round'])
        self.assertNotIn('rank_waits', result)

    def test_different_communicator_is_unverified(self):
        self.fixture(comm=False)
        result = self.result()
        self.assertFalse(result['all_communications_use_world'])
        self.assertNotIn('rank_waits', result)

    def test_waiting_for_child_notification(self):
        self.fixture()
        self.traces[0][2]['tag'] = 1004
        self.traces[0][2]['zero'] = True
        self.save()
        result = self.result()
        self.assertEqual(result['rank_waits']['0']['phase'], 'child_completion')
        self.assertTrue(result['rank_waits']['0']['matching_send_calls_returned'])

    def test_unknown_payload_key_is_refused(self):
        self.fixture()
        self.traces[0][2]['payload'] = 'SECRET'
        self.save()
        self.assertFalse(self.result()['metadata_verified'])

    def test_untracked_request_is_refused(self):
        self.fixture()
        self.traces[0][-1]['requests'] = [0]
        self.save()
        self.assertFalse(self.result()['metadata_verified'])

    def test_missing_hook_record_is_refused(self):
        self.fixture()
        self.traces[0].pop(2)
        self.save()
        self.assertFalse(self.result()['metadata_verified'])

    def test_truncated_trace_is_not_a_missing_message(self):
        self.fixture()
        path = self.root / 'protocol-0.jsonl'
        with path.open('a') as stream:
            stream.write('{"event":')
        self.assertFalse(self.result()['metadata_verified'])

    def test_error_and_completion_are_tracked(self):
        self.fixture()
        self.traces[0].append(dict(event='tested', exchange=12345, requests=[1], ok=False, complete=False))
        self.save()
        self.assertTrue(self.result()['mpi_call_error_seen'])

    def test_completed_receive_group_is_not_still_waiting(self):
        self.fixture()
        self.traces[0].append(dict(event='tested', exchange=12345, requests=[1], ok=True, complete=True))
        self.save()
        self.assertEqual(self.result()['rank_waits']['0']['phase'], 'unknown')

    def test_missing_contact_differs_from_missing_reply(self):
        self.fixture(reply=False)
        del self.traces[1][2:4]
        self.save()
        self.assertFalse(self.result()['rank_waits']['0']['matching_contacts_received'])

    def test_extract_installed_declaration_without_gpu_includes(self):
        fragment = ('#ifndef hypre_EXCHANGE_DATA_HEADER\n#define hypre_EXCHANGE_DATA_HEADER\n'
                    'typedef struct { int unused; } hypre_DataExchangeResponse;\n'
                    'int hypre_DataExchangeList(void);\n#endif\n')
        source = self.root / '_hypre_utilities.h'
        source.write_text('#include <cuda_runtime.h>\n' + fragment + '#include <private_cuda.h>\n')
        result = protocol.exchange_header(self.root, self.root)
        self.assertEqual(result.read_text(), fragment)
        for text in (fragment.replace('#endif', ''), fragment + fragment,
                     fragment.replace('typedef', '#if CUDA\ntypedef'), 'unknown header'):
            source.write_text(text)
            with self.assertRaises(ValueError):
                protocol.exchange_header(self.root, self.root)

    def test_existing_preload_refused_before_compilation(self):
        public = {}
        with patch.object(protocol.subprocess, 'run') as compiler, self.assertRaises(ValueError):
            protocol.prepare(self.root, {}, None, {'LD_PRELOAD': 'PRIVATE'}, public)
        compiler.assert_not_called()
        self.assertEqual(public['failed_check'], 'protocol_existing_preload')

    def test_changed_tracer_invalidates_saved_evidence(self):
        from simple_snapshot_compare import digest
        target = self.root / 'shim.so'
        target.write_bytes(b'original')
        (self.root / 'protocol-build.json').write_text(json.dumps(dict(
            inputs={}, libraries={}, binary=str(target), sha256=digest(target))))
        protocol.recheck(self.root)
        target.write_bytes(b'changed')
        with self.assertRaises(ValueError):
            protocol.recheck(self.root)

    def test_between_exchanges_is_valid_evidence_not_a_broken_trace(self):
        self.fixture()
        self.traces[1].append(dict(event='exchange_return', exchange=12345))
        self.save()
        result = self.result()
        self.assertTrue(result['metadata_verified'])
        self.assertEqual(result['assessment'], 'not_all_ranks_inside_exchange')
        self.assertEqual(result['rank_trace_checks']['0']['state'], 'inside_exchange')
        self.assertEqual(result['rank_trace_checks']['1']['state'], 'between_exchanges')
        self.assertNotIn('rank_waits', result)

    def test_setup_before_exchange_and_setup_return_keep_their_states(self):
        self.fixture()
        self.traces[0] = [dict(event='setup_entry', exchange=0)]
        self.traces[1].extend([dict(event='exchange_return', exchange=12345),
                               dict(event='setup_return', exchange=12345)])
        self.save()
        checks = self.result()['rank_trace_checks']
        self.assertEqual(checks['0']['state'], 'setup_before_first_exchange')
        self.assertEqual(checks['1']['state'], 'setup_returned')

    def test_each_rank_keeps_its_fixed_rejection_reason(self):
        self.fixture()
        self.traces[0][-1]['requests'] = [9999]
        self.traces[1][-1]['event'] = ['SECRET']
        self.save()
        result = self.result()
        self.assertFalse(result['metadata_verified'])
        self.assertEqual(result['rank_trace_checks']['0']['rejection'], 'unposted_request')
        self.assertEqual(result['rank_trace_checks']['1']['rejection'], 'unsupported_event')

    def test_missing_and_partial_files_do_not_hide_the_other_rank(self):
        self.fixture()
        (self.root / 'protocol-0.jsonl').unlink()
        with (self.root / 'protocol-1.jsonl').open('a') as stream:
            stream.write('SECRET')
        checks = self.result()['rank_trace_checks']
        self.assertEqual(checks['0']['rejection'], 'missing_file')
        self.assertEqual(checks['1']['rejection'], 'partial_record')

    def test_message_rejections_distinguish_ids_from_peer_range(self):
        for defect, expected in (('zero', 'nonpositive_message_id'),
                                 ('duplicate', 'duplicate_message_id'),
                                 ('peer', 'peer_outside_communicator')):
            with self.subTest(defect=defect):
                self.fixture()
                if defect == 'zero':
                    self.traces[0][2]['id'] = 0
                elif defect == 'duplicate':
                    self.traces[0].insert(3, dict(self.traces[0][2]))
                else:
                    self.traces[0][2]['peer'] = 12345
                self.save()
                result = self.result()
                self.assertFalse(result['metadata_verified'])
                self.assertEqual(result['rank_trace_checks']['0']['rejection'], expected)
                self.assertTrue(result['rank_trace_checks']['1']['records_valid'])

    def test_invalid_peer_context_keeps_return_and_integer_width_checks(self):
        self.fixture()
        self.traces[0][2]['peer'] = 2**40 + 1
        self.save()
        result = self.result()
        check = result['rank_trace_checks']['0']['invalid_peer_checks']
        self.assertEqual(check['call_kinds'], ['irecv'])
        self.assertTrue(check['all_use_exchange_communicator'])
        self.assertTrue(check['all_have_unique_return_records'])
        self.assertTrue(check['all_calls_reported_success'])
        self.assertTrue(check['any_peer_exceeds_signed_32_bit'])
        self.assertTrue(check['all_low_32_bits_are_valid_ranks'])
        self.assertFalse(check['mpi_argument_corruption_verified'])
        self.assertFalse(check['abi_mismatch_verified'])
        self.assertFalse(result['metadata_verified'])
        self.assertNotIn(str(2**40 + 1), json.dumps(result))

    def test_invalid_peer_in_another_communicator_is_not_assumed_a_world_rank(self):
        self.fixture()
        self.traces[0][2].update(peer=12345, same_comm=False)
        self.traces[0][3]['ok'] = False
        self.save()
        result = self.result()
        check = result['rank_trace_checks']['0']['invalid_peer_checks']
        self.assertFalse(check['all_use_exchange_communicator'])
        self.assertFalse(check['all_calls_reported_success'])
        self.assertFalse(check['any_peer_exceeds_signed_32_bit'])
        self.assertFalse(check['all_low_32_bits_are_valid_ranks'])
        self.assertFalse(result['metadata_verified'])

    def test_missing_or_duplicate_return_cannot_certify_invalid_peer_success(self):
        for mode in ('missing', 'duplicate_return', 'duplicate_message'):
            with self.subTest(mode=mode):
                self.fixture()
                self.traces[0][2]['peer'] = 12345
                if mode == 'missing':
                    del self.traces[0][3]
                else:
                    index = 3 if mode == 'duplicate_return' else 2
                    self.traces[0].insert(index, dict(self.traces[0][index]))
                self.save()
                result = self.result()
                check = result['rank_trace_checks']['0']['invalid_peer_checks']
                self.assertFalse(check['all_have_unique_return_records'])
                self.assertIsNone(check['all_calls_reported_success'])
                self.assertFalse(result['metadata_verified'])

    def test_valid_messages_do_not_get_invalid_peer_checks(self):
        self.fixture()
        result = self.result()
        for rank in ('0', '1'):
            self.assertNotIn('invalid_peer_checks', result['rank_trace_checks'][rank])


if __name__ == '__main__':
    unittest.main()
