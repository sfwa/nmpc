import os
import platform
import ctypes
import subprocess
import distutils.command.build_py
from distutils.core import setup


class build_cnmpc(distutils.command.build_py.build_py):
    description = """Build the CNMPC shared library"""

    def run(self):
        subprocess.call("cmake -DCMAKE_BUILD_TYPE=Debug . && make cnmpc " +
                        "&& cp -r c ./python/nmpc/",
                        shell=True,
                        cwd=os.path.dirname(os.path.abspath(__file__)))
        subprocess.call("cmake -DCMAKE_BUILD_TYPE=Debug . " +
                        " && make c66nmpc && cp -r ccs-c66x ./python/nmpc/",
                        shell=True,
                        cwd=os.path.dirname(os.path.abspath(__file__)))
        self.data_files = self.get_data_files()
        distutils.command.build_py.build_py.run(self)


setup(
    name="nmpc",
    url="https://github.com/sfwa/nmpc",
    author="Daniel Dyer",
    author_email="",
    version="1.0.0",
    description="NMPC library for UAV model predictive control",
    long_description=open("README.md").read(),
    package_dir={"": "python"},
    packages=["nmpc"],
    package_data={"nmpc": ["c/cnmpc.dll", "c/libcnmpc.so", "c/libcnmpc.dylib",
                           "ccs-c66x/c66nmpc.dll", "ccs-c66x/libc66nmpc.so",
                           "ccs-c66x/libc66nmpc.dylib"]},
    license="MIT License",
    classifiers=[
        "Development Status :: 4 - Beta",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
        "Topic :: Software Development :: Libraries"
    ],
    cmdclass={"build_py": build_cnmpc}
)
