import os
import sys

from setuptools import setup


DEBUG_BUILD = os.environ.get("SOURMASH_DEBUG") == "1"
NO_BUILD = os.environ.get("NO_BUILD") == "1"


def find_dylib_no_build(name, paths):
    to_find = None
    if sys.platform == 'darwin':
        to_find = f'lib{name}.dylib'
    elif sys.platform == 'win32':
        to_find = f'{name}.dll'
    else:
        to_find = f'lib{name}.so'

    for path in paths.split(":"):
        for filename in os.listdir(path):
            if filename == to_find:
                return os.path.join(path, filename)

        raise LookupError('dylib %r not found' % name)


def find_dylib(build, target):
    cargo_target = os.environ.get("CARGO_BUILD_TARGET")
    if cargo_target:
        in_path = "target/%s/%s" % (cargo_target, target)
    else:
        in_path = "target/%s" % target
    return build.find_dylib("sourmash", in_path=in_path)


def build_native(spec):
    cmd = ["cargo", "build",
           "--manifest-path", "src/core/Cargo.toml",
            # "--features", "parallel",
           "--lib"]

    target = "debug"
    if not DEBUG_BUILD:
        cmd.append("--release")
        target = "release"

    if NO_BUILD:
        dylib = lambda: find_dylib_no_build("sourmash", os.environ["DYLD_LIBRARY_PATH"])
        header_filename = lambda: "include/sourmash.h"
    else:
        build = spec.add_external_build(cmd=cmd, path=".")
        dylib = lambda: find_dylib(build, target)
        header_filename=lambda: build.find_header("sourmash.h", in_path="include")

    rtld_flags = ["NOW"]
    if sys.platform == "darwin":
        rtld_flags.append("NODELETE")
    spec.add_cffi_module(
        module_path="sourmash._lowlevel",
        dylib=dylib,
        header_filename=header_filename,
        rtld_flags=rtld_flags,
    )

setup(
  milksnake_tasks=[build_native],
  package_dir={"": "src"},
)
