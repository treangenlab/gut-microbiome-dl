'''
To run this setup file, enter the subdirectory that this package sits in
(i.e. tipp_dev_utils/ext) and run the following:

    python setup2.py build_ext --inplace

That will compile the module into a file called "alpha.<something>" where
the <something> depends on the OS and probably the compiler. On windows,
with Visual Studio, it is "alpha.cp36-win_amd64.pyd". From the main
project folder, you can then use:

    import alpha

to import the module and use the functions. Try help(alpha) for the
docstrings.
'''

# from distutils.core import setup, Extension
from setuptools import setup, Extension, Command
# import numpy.distutils.misc_util
import numpy as np

c_alpha = Extension("devtest", ["_devtest.c", "devtest.c"])
c_kmers = Extension("kmers", ["_kmers.c", "kmers.c"])

setup(
    ext_modules=[c_kmers, c_alpha],
    # include_dirs=numpy.distutils.misc_util.get_numpy_include_dirs(),
    include_dirs=np.get_include(),
)
