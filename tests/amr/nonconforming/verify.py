#!/usr/bin/env python3
"""Check the complete immutable source kit before compilation or submission."""
import hashlib
import json
from pathlib import Path
import sys


def verify(root):
    manifest = json.loads((root / 'manifest.json').read_text())
    expected = manifest['files']
    actual = set()
    for path in root.rglob('*'):
        if path.is_symlink():
            raise ValueError(f'Symlink in kit: {path}')
        if path.is_file() and path != root / 'manifest.json':
            actual.add(path.relative_to(root).as_posix())
    if actual != set(expected):
        raise ValueError(f'Kit inventory mismatch: missing={set(expected)-actual}, extra={actual-set(expected)}')
    for name, digest in expected.items():
        if hashlib.sha256((root / name).read_bytes()).hexdigest() != digest:
            raise ValueError(f'Changed kit file: {name}')
    print(f'PASS: {len(expected)} sealed files')


if __name__ == '__main__':
    verify(Path(sys.argv[1]))
