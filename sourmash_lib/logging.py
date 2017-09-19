from __future__ import print_function
import sys
from io import StringIO

def print_results(s, *args, **kwargs):
    print(s.format(*args, **kwargs), file=sys.stdout)
    sys.stdout.flush()


_quiet = False
def set_quiet(val):
    global _quiet
    _quiet = bool(val)


def notify(s, *args, **kwargs):
    "A simple logging function => stderr."
    if not _quiet:
        print(s.format(*args, **kwargs), file=sys.stderr,
              end=kwargs.get('end', u'\n'))
        if kwargs.get('flush'):
            sys.stderr.flush()


def error(s, *args, **kwargs):
    "A simple error logging function => stderr."
    print(s.format(*args, **kwargs), file=sys.stderr)
    if kwargs.get('flush'):
        sys.stderr.flush()


def test_notify():
    global _quiet

    qsave = _quiet
    saveerr, sys.stderr = sys.stderr, StringIO()
    try:
        _quiet = False
        notify(u'hello, world')
    finally:
        _quiet = qsave
        saveerr, sys.stderr = sys.stderr, saveerr

    print(type(saveerr))
    assert 'hello, world\n' in saveerr.getvalue()


def test_notify_flush():
    global _quiet

    qsave = _quiet
    saveerr, sys.stderr = sys.stderr, StringIO()
    try:
        _quiet = False
        notify(u'hello, world', flush=True)
    finally:
        _quiet = qsave
        saveerr, sys.stderr = sys.stderr, saveerr

    print(type(saveerr))
    assert 'hello, world' in saveerr.getvalue()


def test_notify_end():
    global _quiet

    qsave = _quiet
    saveerr, sys.stderr = sys.stderr, StringIO()
    try:
        _quiet = False
        notify(u'hello, world', end=u'FOO')
    finally:
        _quiet = qsave
        saveerr, sys.stderr = sys.stderr, saveerr

    print(type(saveerr))
    assert 'hello, worldFOO' in saveerr.getvalue()


def test_notify_quiet():
    global _quiet

    qsave = _quiet
    saveerr, sys.stderr = sys.stderr, StringIO()
    try:
        _quiet = True
        notify(u'hello, world')
    finally:
        _quiet = qsave
        saveerr, sys.stderr = sys.stderr, saveerr

    print(type(saveerr))
    assert 'hello, world' not in saveerr.getvalue()


def test_error():
    global _quiet

    qsave = _quiet
    saveerr, sys.stderr = sys.stderr, StringIO()
    try:
        _quiet = False
        error(u'hello, world')
    finally:
        _quiet = qsave
        saveerr, sys.stderr = sys.stderr, saveerr

    print(type(saveerr))
    assert 'hello, world\n' in saveerr.getvalue()


def test_error_flush():
    global _quiet

    qsave = _quiet
    saveerr, sys.stderr = sys.stderr, StringIO()
    try:
        _quiet = False
        error(u'hello, world', flush=True)
    finally:
        _quiet = qsave
        saveerr, sys.stderr = sys.stderr, saveerr

    print(type(saveerr))
    assert 'hello, world' in saveerr.getvalue()


def test_error_quiet():
    # error should still output even if _quiet is True
    global _quiet

    qsave = _quiet
    saveerr, sys.stderr = sys.stderr, StringIO()
    try:
        _quiet = True
        error(u'hello, world')
    finally:
        _quiet = qsave
        saveerr, sys.stderr = sys.stderr, saveerr

    print(type(saveerr))
    assert 'hello, world' in saveerr.getvalue()
