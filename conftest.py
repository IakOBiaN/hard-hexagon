"""pytest configuration for the hard hexagon test suite.

Putting the repository root on ``sys.path`` lets the tests do a plain
``import hardhexagon`` regardless of the directory pytest is invoked from.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
