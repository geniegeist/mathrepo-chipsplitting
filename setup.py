from setuptools import setup, Extension, find_packages
from Cython.Build import cythonize
import numpy  

# Define the extension module
# The first argument is the full "dotted" path to the module,
# just as you would import it in Python.
ext_modules = [
    Extension(
        "chipsplitting.solver_ext",
        ["chipsplitting/solver_ext.pyx"]
    )
]

setup(
    name="chipsplitting",
    version="0.1.0",
    packages=find_packages(),

    # Compile the Cython files
    ext_modules=cythonize(
        ext_modules,
        compiler_directives={'language_level': "3"}
    ),

    # It tells the C compiler where to find NumPy's header files.
    include_dirs=[numpy.get_include()],

    # It's good practice to keep this for C extensions
    zip_safe=False,
)
