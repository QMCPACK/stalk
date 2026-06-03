#!/usr/bin/env python

from stalk.io.stalk_logger import StalkLogger

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"


def test_StalkLogger(tmp_path):

    # Test default init
    logger = StalkLogger()
    assert logger.log_level == 1
    assert logger.filename is None

    # Test logging to console
    logger.log("Test message 1", level=1)  # Should print
    logger.log("Test message 2", level=2)  # Should not print

    # Test logging to file
    log_file = tmp_path / "test_log.txt"
    logger_file = StalkLogger(log_level=2, filename=str(log_file))
    logger_file.log("Test message 3", level=1)  # Should log
    logger_file.log("Test message 4", level=2)  # Should log
    logger_file.log("Test message 5", level=3)  # Should not log
    with open(log_file, 'r') as f:
        lines = f.read().splitlines()
    assert lines == ["Test message 3", "Test message 4"]

# end def
