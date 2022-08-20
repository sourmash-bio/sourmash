import tempfile


def to_memmap(array):
    """Write a memory mapped array
    Create a memory-map to an array stored in a binary file on disk.
    Memory-mapped files are used for accessing small segments of
    large files on disk, without reading the entire file into memory.
    :param np.array array to memory map
    :return: np.array large_memmap memory mapped array
    :return: str filename name of the file that memory mapped array is written to
    """
    import numpy as np

    filename = tempfile.NamedTemporaryFile(prefix="array", suffix=".mmap", delete=False).name
    shape = array.shape
    f = np.memmap(filename, mode='w+', shape=shape, dtype=array.dtype)
    f[:] = array[:]
    del f
    large_memmap = np.memmap(filename, dtype=array.dtype, shape=shape)
    return large_memmap, filename
