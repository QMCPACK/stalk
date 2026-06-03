#!/usr/bin env python3

from pathlib import Path

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"


class StalkLogger:
    log_level: int
    filename: str

    def __init__(self, log_level=1, filename=None):
        self.log_level = log_level
        self.filename = filename

        if filename is not None:
            p = Path(filename)
            p.parent.mkdir(parents=True, exist_ok=True)
            with open(filename, 'w') as f:
                f.write('')  # Clear the file
            # end with
        # end if
    # end def

    def log(self, message, level=1):
        if level <= self.log_level:
            if self.filename is not None:
                with open(self.filename, 'a') as f:
                    f.write(message + '\n')
                # end with
            else:
                print(message)
            # end if
        # end if
    # end def

# end class
