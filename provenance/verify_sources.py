#!/usr/bin/env python3
"""Verify the sealed extraction source manifest using only the standard library.

No original checkout is needed. --source-root additionally checks that the
original workspace files still match the pre-extraction snapshot.
"""
import argparse
import hashlib
import json
from pathlib import Path


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-root', type=Path)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    manifest = json.loads((root / 'provenance/extraction-manifest.json').read_text())
    errors = []
    for entry in manifest['files']:
        path = root / entry['destination']
        if not path.is_file() or digest(path) != entry['destination_sha256']:
            errors.append('Changed or missing extracted file: ' + entry['destination'])
    if args.source_root:
        snapshot = json.loads((root / 'provenance/source-snapshot.json').read_text())
        for repo in snapshot['repositories'].values():
            for entry in repo['files']:
                path = args.source_root / repo['source'] / entry['path']
                if not path.is_file() or digest(path) != entry['sha256']:
                    errors.append('Changed original file: ' + str(path))
    if errors:
        raise SystemExit('\n'.join(errors))
    print('PASS: {} imported files match the sealed extraction manifest.'.format(
        len(manifest['files'])))
    if args.source_root:
        print('PASS: original source snapshots remain unchanged.')


if __name__ == '__main__':
    main()
