#!/usr/bin/env python3
"""Compare the exact failure dumps across all ranks of a public-cube case."""
from collections import Counter
import csv
import json
from pathlib import Path
import re
import sys


def compare(directory, ranks):
    inventories = []
    for side in ('before', 'after'):
        files = sorted(directory.glob(side + '-*-rank*.csv'))
        found = [int(re.search(r'-rank(\d+)\.csv$', p.name)[1]) for p in files]
        if sorted(found) != list(range(ranks)):
            raise ValueError(f'{side}: incomplete or duplicated rank inventory: {found}')
        cells = Counter()
        for path in files:
            with path.open() as stream:
                for row in csv.DictReader(stream):
                    cells[tuple(int(row[f'key{i}']) for i in range(4))] += 1
        inventories.append(cells)
    before, after = inventories
    return {'before_count': sum(before.values()), 'after_count': sum(after.values()),
            'removed': [{'keys': list(k), 'count': v} for k, v in sorted((before - after).items())],
            'added': [{'keys': list(k), 'count': v} for k, v in sorted((after - before).items())]}


if __name__ == '__main__':
    directory = Path(sys.argv[1])
    result = compare(directory, int(sys.argv[2]))
    (directory / 'cell-difference.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))
