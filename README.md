User Guide
==========
Install via pip: `pip install flame-code`, libpython3.so is also needed, install `libpython3-dev` on Debian (and its derivatives), or `python3-devel` on RPM-based OS. See the following sections for developers' guide.

Documentation
=============

Generated [documentation](http://frib-high-level-controls.github.io/FLAME)
including [getting started](http://frib-high-level-controls.github.io/FLAME/gettingstarted.html) guide.

Report [bugs through github.com](https://github.com/frib-high-level-controls/FLAME/issues)

Pre-requisites
==============

Needs boost headers.  Also the boost-system and boost-python libraries.
Also python and numpy headers.
The nosetests test runner is used for 'make test' if present.

```sh
apt-get install libboost-dev libboost-system-dev \
 libboost-thread-dev libboost-filesystem-dev \
 libboost-regex-dev libboost-program-options-dev \
 libboost-test-dev \
 build-essential cmake bison flex cppcheck git libhdf5-dev \
 python-numpy python-nose python3-numpy python3-nose
```

Supports python 2.7 and 3.4

Building
========

```sh
git clone https://github.com/frib-high-level-controls/FLAME.git flame
mkdir flame/build
cd flame/build
cmake ..
make
```

To build with a specific python version, change the cmake invokation to:

```sh
cmake .. -DPYTHON_EXECUTABLE=/usr/bin/python3.4
```

<p><a href="https://travis-ci.org/frib-high-level-controls/FLAME">CI status
<img src="https://travis-ci.org/frib-high-level-controls/FLAME.svg" alt="Build Status Images">
</a></p>

Running tests
-------------

```sh
make test
```

Please attach ```Testing/Temporary/LastTest.log``` when reporting test failures.
