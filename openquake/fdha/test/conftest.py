import os
from pathlib import Path


def _disable_numba_disk_cache_for_editable_openquake():
    """Avoid numba cache locator failures for editable OpenQuake checkouts."""
    try:
        import numba
    except Exception:
        return

    original_njit = numba.njit

    def njit_without_disk_cache(*args, **kwargs):
        kwargs["cache"] = False
        return original_njit(*args, **kwargs)

    numba.njit = njit_without_disk_cache


_disable_numba_disk_cache_for_editable_openquake()
os.environ.setdefault("MPLCONFIGDIR", str(Path("/tmp") / "pfdha-matplotlib"))

def pytest_sessionstart(session):
    """
    Change working directory at session start to this test folder,
    so all relative paths in test modules resolve correctly before collection.
    """
    test_dir = os.path.dirname(__file__)
    os.chdir(test_dir)
