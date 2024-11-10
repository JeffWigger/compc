# Contributing to CompIntC

## Set-Up the Project for Development

These instructions are for Fedora Linux.
Make sure you have either clang or gcc installed, as well as support for OpenMP:
```
sudo dnf install clang gcc libomp-devel
export CXX=g++
```
Now, install Google Test:
```
cd .. # Assuming you are in the root directory of this GitHub project
git clone https://github.com/google/googletest.git --branch v1.15.2
cd googletest
cmake -Bbuild -DCMAKE_INSTALL_PREFIX=~/.local # Or another staging area
cmake --build build --config Release
cmake --build build --target install --config Release
```

Lastly, you can test your changes by running:
```
make test
```
