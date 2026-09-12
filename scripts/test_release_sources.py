"""Small offline integrity controls for the source-material collector."""
import copy
import hashlib
from pathlib import Path
import tempfile
import unittest

from collect_release_sources import download, sources_in, verify


class IntegrityControls(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.output = Path(self.directory.name)
        (self.output / 'source.tar').write_bytes(b'pinned source')
        (self.output / 'recipe.patch').write_bytes(b'original patch')
        self.index = {
            'artifacts': {'source': {'path': 'source.tar',
                'sha256': hashlib.sha256(b'pinned source').hexdigest()}},
            'material_files': {'recipe.patch': {
                'sha256': hashlib.sha256(b'original patch').hexdigest()}},
        }

    def test_valid_archive_and_patch(self):
        self.assertEqual(verify(self.index, self.output),
                         {'passed': True, 'checks': 2, 'failures': []})

    def test_missing_patch_fails(self):
        (self.output / 'recipe.patch').unlink()
        self.assertEqual(verify(self.index, self.output)['failures'], ['recipe.patch'])

    def test_corrupt_archive_fails(self):
        (self.output / 'source.tar').write_bytes(b'changed source')
        self.assertFalse(verify(self.index, self.output)['passed'])

    def test_outside_path_fails(self):
        self.index['artifacts']['source']['path'] = '../outside-source.tar'
        self.assertFalse(verify(self.index, self.output)['passed'])

    def test_official_sha512_mismatch_fails(self):
        self.index['artifacts']['source']['expected_sha512'] = '0' * 128
        self.assertFalse(verify(self.index, self.output)['passed'])

    def test_sha512_cache_needs_no_download(self):
        entry = copy.deepcopy(self.index['artifacts']['source'])
        del entry['sha256']
        entry['expected_sha512'] = hashlib.sha512(b'pinned source').hexdigest()
        self.assertEqual(download(entry, self.output)['status'], 'verified')

    def test_rattler_nested_sources_preserved(self):
        source = {'url': 'https://example.test/source', 'sha256': '1' * 64}
        self.assertEqual(sources_in({'cache': [{'source': [source]}]}), [source])


if __name__ == '__main__':
    unittest.main()
