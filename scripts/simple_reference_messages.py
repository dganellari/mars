"""Match private errors to literal diagnostics in the pinned public reference."""

import ast
import io
import re
import subprocess
import tarfile

REFERENCE_REVISION = '0d69041ba1afda63e9e4328d9e0d9834bba37756'
SOURCE_ROOTS = ('src/', 'external/nalu/', 'external/gplotpp/')
LEXEME = re.compile(r'//[^\n]*|/\*.*?\*/|"(?:\\.|[^"\\])*"|\'(?:\\.|[^\'\\])*\'', re.S)
CALL = re.compile(r'\b(?:errorMsg|STK_Throw\w*|throw\s+std::\w+)\s*\(')


def words(text):
    return tuple(re.findall(r'[A-Za-z0-9_]+', text.lower()))


def diagnostic_literals(text):
    # Mask comments and literals before balancing call parentheses. This is a
    # literal-message index, not a C++ evaluator; dynamically built errors may miss.
    lexemes = list(LEXEME.finditer(text))
    masked = list(text)
    for token in lexemes:
        masked[token.start():token.end()] = ['\n' if c == '\n' else ' ' for c in token.group()]
    masked = ''.join(masked)
    for call in CALL.finditer(masked):
        end, depth = call.end(), 1
        while end < len(masked) and depth:
            depth += (masked[end] == '(') - (masked[end] == ')')
            end += 1
        if depth:
            continue
        line = text.count('\n', 0, call.start()) + 1
        fragments, previous_end = [], call.end()
        for token in lexemes:
            if token.start() < call.end() or token.end() > end:
                continue
            if token.group().startswith(('//', '/*')):
                continue
            if not token.group().startswith('"'):
                previous_end = token.end()
                continue
            try:
                literal = ast.literal_eval(token.group())
            except (ValueError, SyntaxError):
                continue
            if fragments and not masked[previous_end:token.start()].strip():
                fragments[-1] += literal
            else:
                fragments.append(literal)
            previous_end = token.end()
        for fragment in fragments:
            tokens = words(fragment)
            if len(fragment) >= 16 and 3 <= len(tokens) <= 64:
                yield tokens, line


def read_catalog(source):
    # Read immutable, allowlisted Git content, never working-tree/private files.
    result = subprocess.run(['git', '-C', str(source), 'archive', REFERENCE_REVISION,
                             'main.cpp', 'src', 'external/nalu', 'external/gplotpp'],
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if result.returncode:
        raise ValueError('public reference source unavailable')
    catalog = {}
    with tarfile.open(fileobj=io.BytesIO(result.stdout), mode='r:') as archive:
        for item in archive:
            if (not item.isfile() or '..' in item.name.split('/')
                    or not (item.name == 'main.cpp' or item.name.startswith(SOURCE_ROOTS))
                    or not item.name.endswith(('.cpp', '.hpp', '.h', '.C', '.cxx', '.cc'))):
                continue
            text = archive.extractfile(item).read().decode('utf-8', errors='replace')
            for tokens, line in diagnostic_literals(text):
                catalog.setdefault(tokens, set()).add(item.name + ':' + str(line))
    if not catalog:
        raise ValueError('public reference diagnostics unavailable')
    return catalog


class MessageMatcher:
    def __init__(self, catalog):
        self.catalog = catalog
        self.lengths = sorted(set(map(len, catalog)))
        self.tail = ()
        self.matches = {}

    def feed(self, line):
        tokens = self.tail + words(line)
        for size in self.lengths:
            for start in range(max(0, len(self.tail) - size + 1), len(tokens) - size + 1):
                for location in self.catalog.get(tokens[start:start + size], ()):
                    self.matches[location] = max(size, self.matches.get(location, 0))
        self.tail = tokens[-63:]

    def result(self):
        matches = sorted(self.matches, key=lambda location: (-self.matches[location], location))
        return dict(reference_message_candidates=matches[:20],
                    reference_message_candidates_truncated=len(matches) > 20,
                    reference_message_match_scope='literal_fragment_candidates_not_stack_trace')
