# nmpc

Non-linear Model Predictive Control library for UAV control.


## Structure


## Configuration


## Building

Requires `cmake` version 2.8.7 or higher.

Create a build directory outside the source tree, then use cmake to generate
the makefile.

`mkdir nmpc_build`

`cd nmpc_build`

`cmake /path/to/nmpc`

Now, build the library using the `make` command. An appropriate version of
Eigen will be downloaded automatically.

To build the dynamic library, run `make cnmpc`. A dynamic library appropriate
for the host platform should be built.


## Testing


## Python module installation

Requires `cmake` version 2.8.7 or higher.

Run `python setup.py install` to build the C shared library and install the
Python interface (the `nmpc` module) in your `site-packages` directory.

Alternatively, just run `pip install https://github.com/sfwa/nmpc/archive/master.zip#egg=nmpc-1.0.0`
to download and install.
