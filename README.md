Pre-requisites
==============

Needs boost headers.  Also the boost-system and boost-python libraries.
Also python and numpy headers.
The nosetests test runner is used for 'make test' if present.

```sh
apt-get install libboost-dev libboost-system-dev \
libboost-python-dev python-dev python-nose python-numpy \
cmake build-essential bison flex git
```

Supports python 2.6, 2.7, and 3.4

Building
========

```sh
git clone https://github.com/mdavidsaver/jmbgsddb.git
cd jmbgsddb
mkdir build
cd build
cmake ..
make
```

To build with a specific python version, call cmake with:

```sh
cmake .. -DPYTHON_EXECUTABLE=/usr/bin/python3.4
```

Running tests
-------------

```sh
make test
```

Running examples
----------------

```sh
./examples/example_linear
```
