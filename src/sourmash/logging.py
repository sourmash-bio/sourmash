import sys
from io import StringIO

_quiet = False
_debug = False


def set_quiet(val, print_debug=False):
    global _quiet, _debug
    _quiet = bool(val)
    _debug = bool(print_debug)


def print_results(s, *args, **kwargs):
    if _quiet:
        return

    print(s.format(*args, **kwargs), file=sys.stdout)
    sys.stdout.flush()


def notify(s, *args, **kwargs):
    "A simple logging function => stderr."
    if _quiet:
        return

    print("\r\033[K", end="", file=sys.stderr)
    print(s.format(*args, **kwargs), file=sys.stderr, end=kwargs.get("end", "\n"))
    if kwargs.get("flush"):
        sys.stderr.flush()


def debug(s, *args, **kwargs):
    "A debug logging function => stderr."
    if _quiet or not _debug:
        return

    print("\r\033[K", end="", file=sys.stderr)
    print(s.format(*args, **kwargs), file=sys.stderr, end=kwargs.get("end", "\n"))
    if kwargs.get("flush"):
        sys.stderr.flush()


def debug_literal(s, *args, **kwargs):
    "A debug logging function => stderr."
    if _quiet or not _debug:
        return

    print("\r\033[K", end="", file=sys.stderr)
    print(s, file=sys.stderr, end=kwargs.get("end", "\n"))
    if kwargs.get("flush"):
        sys.stderr.flush()


def error(s, *args, **kwargs):
    "A simple error logging function => stderr."
    print("\r\033[K", end="", file=sys.stderr)
    print(s.format(*args, **kwargs), file=sys.stderr)
    if kwargs.get("flush"):
        sys.stderr.flush()


def test_notify():
    global _quiet

    qsave = _quiet
    saveerr, sys.stderr = sys.stderr, StringIO()
    try:
        _quiet = False
        notify("hello, world")
    finally:
        _quiet = qsave
        saveerr, sys.stderr = sys.stderr, saveerr

    print(type(saveerr))
    assert "hello, world\n" in saveerr.getvalue()


def test_notify_flush():
    global _quiet

    qsave = _quiet
    saveerr, sys.stderr = sys.stderr, StringIO()
    try:
        _quiet = False
        notify("hello, world", flush=True)
    finally:
        _quiet = qsave
        saveerr, sys.stderr = sys.stderr, saveerr

    print(type(saveerr))
    assert "hello, world" in saveerr.getvalue()


def test_notify_end():
    global _quiet

    qsave = _quiet
    saveerr, sys.stderr = sys.stderr, StringIO()
    try:
        _quiet = False
        notify("hello, world", end="FOO")
    finally:
        _quiet = qsave
        saveerr, sys.stderr = sys.stderr, saveerr

    print(type(saveerr))
    assert "hello, worldFOO" in saveerr.getvalue()


def test_notify_quiet():
    global _quiet

    qsave = _quiet
    saveerr, sys.stderr = sys.stderr, StringIO()
    try:
        _quiet = True
        notify("hello, world")
    finally:
        _quiet = qsave
        saveerr, sys.stderr = sys.stderr, saveerr

    print(type(saveerr))
    assert "hello, world" not in saveerr.getvalue()


def test_error():
    global _quiet

    qsave = _quiet
    saveerr, sys.stderr = sys.stderr, StringIO()
    try:
        _quiet = False
        error("hello, world")
    finally:
        _quiet = qsave
        saveerr, sys.stderr = sys.stderr, saveerr

    print(type(saveerr))
    assert "hello, world\n" in saveerr.getvalue()


def test_error_flush():
    global _quiet

    qsave = _quiet
    saveerr, sys.stderr = sys.stderr, StringIO()
    try:
        _quiet = False
        error("hello, world", flush=True)
    finally:
        _quiet = qsave
        saveerr, sys.stderr = sys.stderr, saveerr

    print(type(saveerr))
    assert "hello, world" in saveerr.getvalue()


def test_error_quiet():
    # error should still output even if _quiet is True
    global _quiet

    qsave = _quiet
    saveerr, sys.stderr = sys.stderr, StringIO()
    try:
        _quiet = True
        error("hello, world")
    finally:
        _quiet = qsave
        saveerr, sys.stderr = sys.stderr, saveerr

    print(type(saveerr))
    assert "hello, world" in saveerr.getvalue()
