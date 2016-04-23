Documentation
=============

Generated [documentation](http://frib-high-level-controls.github.io/FLAME)

Pre-requisites
==============

Needs boost headers.  Also the boost-system and boost-python libraries.
Also python and numpy headers.
The nosetests test runner is used for 'make test' if present.

```sh
apt-get install libboost-dev libboost-system-dev \
libboost-thread-dev python-dev python-nose python-numpy \
cmake build-essential bison flex git libhdf5-dev
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

<p><a href="https://travis-ci.org/mdavidsaver/jmbgsddb">CI status
<img src="https://travis-ci.org/mdavidsaver/jmbgsddb.svg" alt="Build Status Images">
</a></p>

Running tests
-------------

```sh
make test
```
