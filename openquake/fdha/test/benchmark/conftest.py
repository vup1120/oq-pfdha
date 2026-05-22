"""
Pytest fixtures for benchmark tests.
"""
import os
import pytest


def pytest_sessionstart(session):
    """
    Change working directory at session start to the test folder,
    so all relative paths in test modules resolve correctly.
    """
    test_dir = os.path.dirname(__file__)
    os.chdir(test_dir)


