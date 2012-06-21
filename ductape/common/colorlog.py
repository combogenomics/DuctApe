#!/usr/bin/env python
"""
ColorLog

Common library

Colorize the StreamHandler log entries

Thanks to dh82 on Stackoverflow (http://stackoverflow.com/a/2532931/1237531)
"""

import logging
import sys

BLACK, RED, GREEN, YELLOW, BLUE, MAGENTA, CYAN, WHITE = range(8)

COLORS = {
    'WARNING'  : YELLOW,
    'INFO'     : WHITE,
    'DEBUG'    : WHITE,
    'CRITICAL' : MAGENTA,
    'ERROR'    : RED,
    'RED'      : RED,
    'GREEN'    : GREEN,
    'YELLOW'   : YELLOW,
    'BLUE'     : BLUE,
    'MAGENTA'  : MAGENTA,
    'CYAN'     : CYAN,
    'WHITE'    : WHITE,
}

RESET_SEQ = "\033[0m"
COLOR_SEQ = "\033[1;%dm"
BOLD_SEQ  = "\033[1m"

class ColorFormatter(logging.Formatter):

    def __init__(self, *args, **kwargs):
        logging.Formatter.__init__(self, *args, **kwargs)
        if not hasattr(sys.stdout, 'isatty') or not sys.stdout.isatty():
            self.use_it = False
        else:
            self.use_it = True

    def format(self, record):
        if not self.use_it:
            message = logging.Formatter.format(self, record)
            message   = message.replace("$RESET", '')\
                           .replace("$BOLD",  '')\
                           .replace("$COLOR", '')
            return message
        levelname = record.levelname
        color     = COLOR_SEQ % (30 + COLORS[levelname])
        message   = logging.Formatter.format(self, record)
        message   = message.replace("$RESET", RESET_SEQ)\
                           .replace("$BOLD",  BOLD_SEQ)\
                           .replace("$COLOR", color)
        for k,v in COLORS.items():
            message = message.replace("$" + k,    COLOR_SEQ % (v+30))\
                             .replace("$BG" + k,  COLOR_SEQ % (v+40))\
                             .replace("$BG-" + k, COLOR_SEQ % (v+40))
        return message + RESET_SEQ
