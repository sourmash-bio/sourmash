from __future__ import print_function
import sys

_quiet = False
def notify(s, *args, **kwargs):
    "A simple logging function => stderr."
    if not _quiet:
        print(s.format(*args, **kwargs), file=sys.stderr,
              end=kwargs.get('end', '\n'))
        if kwargs.get('flush'):
            sys.stderr.flush()


def error(s, *args, **kwargs):
    "A simple error logging function => stderr."
    print(s.format(*args, **kwargs), file=sys.stderr)
    if kwargs.get('flush'):
        sys.stderr.flush()


