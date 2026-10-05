import io
import json
import subprocess
import tarfile
import unittest
from unittest.mock import patch

from simple_reference_messages import MessageMatcher, diagnostic_literals, read_catalog, words, REFERENCE_REVISION
from simple_startup_probe import ReferenceLogState


class ReferenceMessageTests(unittest.TestCase):
    def test_adjacent_literals_and_dynamic_values(self):
        source = '''// errorMsg("ignore this diagnostic entirely");
void f() {
    errorMsg("Error in the expression " /* comment */ "provided at boundary `" + name + "`");
    throw std::runtime_error("reference setup could not continue");
}'''
        literals = list(diagnostic_literals(source))
        self.assertEqual(literals, [
            (words('Error in the expression provided at boundary'), 3),
            (words('reference setup could not continue'), 4)])

    def test_quotes_and_parentheses_do_not_end_call(self):
        source = r'''errorMsg("unsupported option (see documentation) \"here\"");
errorMsg("prefix for a private name: " + fn("ordinary string for this argument"));'''
        self.assertIn((words('unsupported option see documentation here'), 1), list(diagnostic_literals(source)))
        self.assertIn((words('prefix for a private name'), 2), list(diagnostic_literals(source)))

    def test_unknown_private_text_does_not_escape(self):
        matcher = MessageMatcher({words('reference setup could not continue'): {'src/public.cpp:8'}})
        matcher.feed('what(): PRIVATE confidential 123.45')
        result = matcher.result()
        self.assertEqual(result['reference_message_candidates'], [])
        self.assertNotIn('PRIVATE', json.dumps(result))
        self.assertNotIn('123', json.dumps(result))

    def test_multiline_error_and_repeated_ranks(self):
        state = ReferenceLogState({words('Error in the expression provided at boundary'): {'src/model/model.cpp:62'}})
        for line in ['Error in the expression provided at boundary PRIVATE',
                     "terminate called after throwing an instance of 'std::runtime_error'",
                     'what(): Error in the expression', 'provided at boundary PRIVATE',
                     'what(): Error in the expression provided at boundary PRIVATE']:
            state.feed(line)
        result = state.result()
        self.assertEqual(result['reference_message_candidates'], ['src/model/model.cpp:62'])
        self.assertNotIn('PRIVATE', json.dumps(result))

    def test_ordinary_output_is_not_scanned(self):
        state = ReferenceLogState({words('normal output refers to boundary'): {'src/public.cpp:8'}})
        state.feed('normal output refers to boundary PRIVATE')
        self.assertEqual(state.result()['reference_message_candidates'], [])

    def test_candidate_limit_and_deterministic_specificity(self):
        catalog = {words('setup failed here'): {'src/public.cpp:' + str(i) for i in range(25)},
                   words('setup failed here after allocation'): {'src/precise.cpp:1'}}
        matcher = MessageMatcher(catalog)
        matcher.feed('what(): setup failed here after allocation PRIVATE')
        result = matcher.result()
        self.assertTrue(result['reference_message_candidates_truncated'])
        self.assertEqual(len(result['reference_message_candidates']), 20)
        self.assertEqual(result['reference_message_candidates'][0], 'src/precise.cpp:1')

    def test_catalog_uses_only_pinned_source_members(self):
        stream = io.BytesIO()
        with tarfile.open(fileobj=stream, mode='w') as archive:
            for name, text in (('src/public.cpp', 'errorMsg("known public failure in setup");'),
                               ('private.cpp', 'errorMsg("private content must stay hidden");'),
                               ('src/../../private.cpp', 'errorMsg("private content must stay hidden");')):
                data = text.encode()
                item = tarfile.TarInfo(name)
                item.size = len(data)
                archive.addfile(item, io.BytesIO(data))
        with patch('simple_reference_messages.subprocess.run', return_value=
                   subprocess.CompletedProcess([], 0, stream.getvalue(), b'')) as run:
            catalog = read_catalog('/synthetic/source')
        self.assertEqual(catalog, {words('known public failure in setup'): {'src/public.cpp:1'}})
        self.assertIn(REFERENCE_REVISION, run.call_args[0][0])

    def test_source_errors_hide_process_output(self):
        with patch('simple_reference_messages.subprocess.run', return_value=
                   subprocess.CompletedProcess([], 1, b'', b'PRIVATE')):
            with self.assertRaisesRegex(ValueError, '^public reference source unavailable$'):
                read_catalog('/synthetic/source')


if __name__ == '__main__':
    unittest.main()
