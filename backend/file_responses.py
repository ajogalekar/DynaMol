"""Owned temporary downloads must be removed on every response exit path."""
from pathlib import Path
from starlette.responses import FileResponse


class TemporaryFileResponse(FileResponse):
    """Stream an owned snapshot, including ranges, then unlink it unconditionally.

    A background task is insufficient: FileResponse returns before its task on
    rejected Range headers and an interrupted send can also skip that task.
    """

    async def __call__(self, scope, receive, send):
        # Keep lifetime ownership here rather than passing a soon-to-be-deleted
        # pathname to an optional asynchronous server file-offload extension.
        scope = {**scope, "extensions": {key: value for key, value in scope.get("extensions", {}).items() if key != "http.response.pathsend"}}
        try:
            await super().__call__(scope, receive, send)
        finally:
            # One synchronous unlink also runs when the request task is already
            # cancelled; there is no cancellation-sensitive await before cleanup.
            Path(self.path).unlink(missing_ok=True)
