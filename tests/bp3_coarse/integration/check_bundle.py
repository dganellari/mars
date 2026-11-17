import sys
from pathlib import Path
from common import verify_bundle

record = verify_bundle(Path(sys.argv[1]).resolve())
print('PASS bundle: {} files'.format(len(record['files'])))
