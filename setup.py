from setuptools import setup
from setuptools_rust import RustExtension, Binding

setup(
  package_dir={"": "src"},
  rust_extensions=[
    RustExtension("sourmash._lowlevel__lib",
                  py_limited_api="auto",
                  path="src/core/Cargo.toml",
                  binding=Binding.NoBinding),
  ],
  cffi_modules=["src/sourmash/ffi_build.py:ffibuilder"],
)
