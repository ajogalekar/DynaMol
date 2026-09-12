"""Local browser origins and Host headers are independently checked."""
import pytest
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from fastapi.testclient import TestClient

from backend.main import LocalOriginMiddleware, local_browser_origins


@pytest.mark.parametrize("origin", [
    "https://untrusted.example", "http://127.0.0.1:64059/", "http://127.0.0.1:65536",
    "http://127.0.0.1:0", "http://127.0.0.1:064", "http://127.0.0.1:64059.evil.example",
    "http://localhost:64059", "http://user@127.0.0.1:64059", "http://127.0.0.1:*",
])
def test_packaged_origin_rejects_noncanonical_configuration(origin):
    with pytest.raises(ValueError, match="DYNAMOL_ORIGIN"):
        local_browser_origins(origin)


def test_packaged_origin_serves_browser_and_preflight_without_trusting_other_ports():
    origin = "http://127.0.0.1:64059"
    allowed = local_browser_origins(origin)
    app = FastAPI()
    app.add_middleware(CORSMiddleware, allow_origins=allowed, allow_methods=["GET", "POST"], allow_headers=["Content-Type"])
    app.add_middleware(LocalOriginMiddleware, allowed_origins=allowed)

    @app.api_route("/probe", methods=["GET", "POST"])
    def probe():
        return {"ok": True}

    with TestClient(app, base_url="http://127.0.0.1") as client:
        for method in ("GET", "POST"):
            response = client.request(method, "/probe", headers={"Origin": origin})
            assert response.status_code == 200
            assert response.headers["access-control-allow-origin"] == origin
        response = client.options("/probe", headers={"Origin": origin, "Access-Control-Request-Method": "POST", "Access-Control-Request-Headers": "Content-Type"})
        assert response.status_code == 200
        assert response.headers["access-control-allow-origin"] == origin
        for foreign in ("https://untrusted.example", "null", "http://127.0.0.1:64060"):
            assert client.post("/probe", headers={"Origin": foreign}).status_code == 403
        assert client.get("/probe").status_code == 200


def test_packaged_origin_preserves_source_defaults():
    assert len(local_browser_origins()) == 6
    assert local_browser_origins("http://127.0.0.1:8765") == local_browser_origins()


@pytest.mark.parametrize("host", [
    "127.0.0.1", "localhost", "127.0.0.1:8765", "localhost:5173",
    "127.0.0.1:64059", "localhost:64059",
])
def test_real_app_accepts_local_hosts_without_origin(monkeypatch, host):
    from backend.main import app, jobs

    monkeypatch.setattr(jobs, "health", lambda: {"local_probe": True})
    with TestClient(app, base_url=f"http://{host}") as client:
        response = client.get("/api/health")
    assert response.status_code == 200
    assert response.json() == {"local_probe": True}


@pytest.mark.parametrize("host", [
    "untrusted.example", "untrusted.example:8765", "127.0.0.1.evil.example",
    "localhost.evil.example", "www.localhost", "127.0.0.1@evil.example", "",
])
@pytest.mark.parametrize("method,path", [("GET", "/api/health"), ("POST", "/api/jobs")])
def test_real_app_rejects_foreign_or_missing_host_before_dispatch(monkeypatch, host, method, path):
    from backend.main import app, jobs

    def unexpected_dispatch(*args, **kwargs):
        pytest.fail("Untrusted Host reached the application")

    monkeypatch.setattr(jobs, "health", unexpected_dispatch)
    monkeypatch.setattr(jobs, "submit_job", unexpected_dispatch)
    with TestClient(app, base_url="http://127.0.0.1") as client:
        # A trusted or absent Origin and forwarding headers must not exempt Host.
        for origin in (None, "http://127.0.0.1:8765"):
            headers = {"Host": host, "X-Forwarded-Host": "localhost"}
            if origin:
                headers["Origin"] = origin
            response = client.request(method, path, headers=headers)
            assert response.status_code == 400
            assert "location" not in response.headers


def test_real_app_still_rejects_foreign_origin_on_local_host():
    from backend.main import app

    with TestClient(app, base_url="http://127.0.0.1:8765") as client:
        response = client.get("/api/health", headers={"Origin": "http://untrusted.example"})
    assert response.status_code == 403
