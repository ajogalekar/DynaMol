"""Job downloads remain private snapshots across concurrent downloads and resume."""
import asyncio
import io
import json
import threading
import zipfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pytest

from backend import config, jobs, main, recovery, storage


@pytest.fixture
def stopped_job(tmp_path, monkeypatch):
    jobs_root = tmp_path / "jobs"
    folder = jobs_root / "snapshot-job"
    folder.mkdir(parents=True)
    monkeypatch.setattr(config, "JOBS_DIR", jobs_root)
    job = {
        "id": folder.name,
        "engine": "gromacs",
        "status": "interrupted",
        "created_at": "2026-09-11T00:00:00+00:00",
        "elapsed_seconds": 12,
        "logs": [],
        "config": {},
    }
    storage.atomic_json(folder / "status.json", job)
    (folder / "trajectory.xtc").write_bytes(b"stopped trajectory")
    (folder / "production.cpt").write_bytes(b"checkpoint fixture; no native engine")
    downloads = tmp_path / "downloads"
    downloads.mkdir()
    original_temporary_file = main.tempfile.NamedTemporaryFile

    def isolated_temporary_file(*args, **kwargs):
        return original_temporary_file(*args, dir=downloads, **kwargs)

    monkeypatch.setattr(main.tempfile, "NamedTemporaryFile", isolated_temporary_file)
    return folder, downloads


def stream_response(response, *, headers=(), expected_status=200, fail_body=False):
    """Exercise FileResponse and its cleanup task, rather than invoking cleanup directly."""
    target = Path(response.path)

    async def consume():
        body = bytearray()

        async def receive():
            raise AssertionError("A normal file response should not read the request body")

        async def send(message):
            assert target.is_file(), "The archive must exist until streaming finishes"
            if message["type"] == "http.response.start":
                assert message["status"] == expected_status
            if message["type"] == "http.response.body":
                if fail_body:
                    raise ConnectionError("Injected client disconnect during streaming")
                body.extend(message["body"])

        await response({"type": "http", "method": "GET", "headers": list(headers)}, receive, send)
        return bytes(body)

    return asyncio.run(consume())


def archive_contents(data):
    with zipfile.ZipFile(io.BytesIO(data)) as archive:
        return {name: archive.read(name) for name in archive.namelist()}


def test_concurrent_downloads_have_independent_snapshots_and_response_cleanup(stopped_job):
    folder, downloads = stopped_job
    start = threading.Barrier(3)

    def request_download():
        start.wait(timeout=5)
        return main.download(folder.name)

    with ThreadPoolExecutor(max_workers=2) as pool:
        pending = [pool.submit(request_download) for _ in range(2)]
        start.wait(timeout=5)
        first, second = [future.result(timeout=5) for future in pending]

    first_path, second_path = Path(first.path), Path(second.path)
    assert first_path != second_path
    assert first_path.parent == second_path.parent == downloads
    assert first_path.is_file() and second_path.is_file()

    # Both responses are pending while the job's original outputs change.
    (folder / "trajectory.xtc").write_bytes(b"new trajectory after resume")
    (folder / "production.cpt").unlink()
    (folder / "later-output.txt").write_text("created after the requests")

    first_files = archive_contents(stream_response(first))
    assert not first_path.exists()
    assert second_path.is_file(), "Finishing one download must not delete another response's file"
    second_files = archive_contents(stream_response(second))
    assert not second_path.exists()
    assert first_files == second_files
    assert first_files["trajectory.xtc"] == b"stopped trajectory"
    assert first_files["production.cpt"] == b"checkpoint fixture; no native engine"
    assert "later-output.txt" not in first_files
    assert not list(downloads.iterdir())


def test_resume_waits_for_archive_creation_but_not_for_response_streaming(stopped_job, monkeypatch):
    folder, downloads = stopped_job
    archive_writing = threading.Event()
    release_archive = threading.Event()
    archive_closed = threading.Event()

    class ObservedLock:
        """Report an actual blocked acquisition without relying on timing or sleeps."""
        def __init__(self):
            self.inner = threading.RLock()
            self.resume_thread = None
            self.resume_blocked = threading.Event()

        def __enter__(self):
            if threading.get_ident() == self.resume_thread:
                if not self.inner.acquire(blocking=False):
                    self.resume_blocked.set()
                    self.inner.acquire()
            else:
                self.inner.acquire()
            return self

        def __exit__(self, *exc):
            self.inner.release()

    observed_lock = ObservedLock()
    monkeypatch.setattr(jobs, "_lock", observed_lock)
    original_write, original_close = zipfile.ZipFile.write, zipfile.ZipFile.close

    def gated_write(archive, filename, *args, **kwargs):
        if Path(filename).name == "trajectory.xtc":
            archive_writing.set()
            assert release_archive.wait(timeout=5), "The test must release archive construction"
        return original_write(archive, filename, *args, **kwargs)

    def observed_close(archive):
        original_close(archive)
        if archive.mode == "w" and Path(archive.filename).parent == downloads:
            archive_closed.set()

    monkeypatch.setattr(zipfile.ZipFile, "write", gated_write)
    monkeypatch.setattr(zipfile.ZipFile, "close", observed_close)
    monkeypatch.setattr(recovery, "validate_manifest", lambda *args: {"engine": "gromacs"})
    launches = []

    def launch_without_engine(job_folder, job):
        assert archive_closed.is_set(), "Resume must not mutate outputs while the ZIP is being built"
        assert len(list(downloads.glob("*.zip"))) == 1
        (job_folder / "trajectory.xtc").write_bytes(b"resumed trajectory")
        storage.atomic_json(job_folder / "status.json", job)
        launches.append(job["id"])
        return job

    monkeypatch.setattr(jobs, "_launch_worker", launch_without_engine)

    def request_resume():
        observed_lock.resume_thread = threading.get_ident()
        return jobs.resume_job(folder.name)

    with ThreadPoolExecutor(max_workers=2) as pool:
        download_future = pool.submit(main.download, folder.name)
        try:
            assert archive_writing.wait(timeout=5)
            resume_future = pool.submit(request_resume)
            assert observed_lock.resume_blocked.wait(timeout=5)
            assert not launches
            assert not archive_closed.is_set()
            assert json.loads((folder / "status.json").read_text())["status"] == "interrupted"
        finally:
            release_archive.set()
        response = download_future.result(timeout=5)
        resumed = resume_future.result(timeout=5)

    assert launches == [folder.name]
    assert resumed["status"] == "queued" and len(resumed["restarts"]) == 1
    assert Path(response.path).is_file(), "Resume must proceed while the private response is still pending"
    assert (folder / "trajectory.xtc").read_bytes() == b"resumed trajectory"
    contents = archive_contents(stream_response(response))
    assert contents["trajectory.xtc"] == b"stopped trajectory"
    assert json.loads(contents["status.json"])["status"] == "interrupted"
    assert not list(downloads.iterdir())


def test_archive_construction_failure_removes_its_partial_temp_file(stopped_job, monkeypatch):
    folder, downloads = stopped_job

    def unreadable_output(*args, **kwargs):
        raise OSError("Injected failure reading a job output")

    monkeypatch.setattr(zipfile.ZipFile, "write", unreadable_output)
    with pytest.raises(OSError, match="Injected failure"):
        main.download(folder.name)
    assert not list(downloads.iterdir())
    assert (folder / "trajectory.xtc").read_bytes() == b"stopped trajectory"


@pytest.mark.parametrize("range_header,status", [(b"invalid-range", 400), (b"bytes=999999999999-", 416)])
def test_rejected_range_response_removes_its_snapshot(stopped_job, range_header, status):
    folder, downloads = stopped_job
    response = main.download(folder.name)
    stream_response(response, headers=[(b"range", range_header)], expected_status=status)
    assert not Path(response.path).exists(), "An unsuccessful range request must not leak a complete job archive"
    assert not list(downloads.iterdir())


def test_disconnected_response_removes_its_snapshot(stopped_job):
    folder, downloads = stopped_job
    response = main.download(folder.name)
    with pytest.raises(ConnectionError, match="Injected client disconnect"):
        stream_response(response, fail_body=True)
    assert not Path(response.path).exists(), "A failed send must still release the private archive"
    assert not list(downloads.iterdir())
