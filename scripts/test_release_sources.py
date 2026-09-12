"""Small offline integrity controls for the source-material collector."""
import copy
import hashlib
import io
from pathlib import Path
import tarfile
import tempfile
import unittest

from collect_release_sources import download, pinned_native_git_source, sources_in, verify


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

    def test_git_source_preserves_full_revision_and_rejects_mutable_refs(self):
        source = {'git_url': 'https://github.com/example/library.git', 'git_rev': 'a' * 40}
        self.assertEqual(sources_in({'source': source}), [source])
        mapped, _ = pinned_native_git_source(source, 'library', 'runtime', {}, self.output, self.output)
        self.assertEqual(mapped['url'], 'https://codeload.github.com/example/library/tar.gz/' + 'a' * 40)
        with self.assertRaisesRegex(ValueError, 'full recorded commit'):
            pinned_native_git_source({**source, 'git_rev': 'main'}, 'library', 'runtime', {}, self.output, self.output)

    def test_native_openmm_source_binds_embedded_commit_to_archive_and_recipe(self):
        revision = 'c' * 40
        engines = self.output / 'engines'
        engines.mkdir()
        archive = engines / 'native.tar.gz'
        content = f"git_revision = '{revision}'\n".encode()
        with tarfile.open(archive, 'w:gz') as packed:
            member = tarfile.TarInfo('lib/python3.12/site-packages/openmm/version.py')
            member.size = len(content)
            packed.addfile(member, io.BytesIO(content))
        record = {'archive': archive.name, 'sha256': hashlib.sha256(archive.read_bytes()).hexdigest()}
        manifest = {'engines': {'promod3': record}}
        source = {'git_url': 'https://github.com/openmm/openmm.git', 'git_rev': revision[:7]}
        mapped, _ = pinned_native_git_source(source, 'openmm', 'promod3', manifest, self.output, self.output)
        self.assertEqual(mapped['git_revision'], revision)
        self.assertEqual((self.output / 'shipped-version.py').read_bytes(), content)
        with self.assertRaisesRegex(ValueError, 'conflicts with its recipe'):
            pinned_native_git_source({**source, 'git_rev': 'a' * 7}, 'openmm', 'promod3', manifest, self.output, self.output)
        record['sha256'] = '0' * 64
        with self.assertRaisesRegex(ValueError, 'does not match its manifest'):
            pinned_native_git_source(source, 'openmm', 'promod3', manifest, self.output, self.output)


if __name__ == '__main__':
    unittest.main()
