Pre-requisites
==============

Needs boost headers.  Also the boost-system and boost-python libraries.
Also python and numpy headers.
The nosetests test runner is used for 'make test' if present.

```sh
apt-get install libboost-dev libboost-system-dev libboost-python-dev \
 python-dev cmake python-nose python-numpy git
```

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
