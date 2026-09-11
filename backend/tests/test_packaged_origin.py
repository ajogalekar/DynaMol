"""The packaged app has one dynamically allocated, explicitly trusted local origin."""
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

    with TestClient(app) as client:
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
