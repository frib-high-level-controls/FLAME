Installation
=============

Build from source code
----------------------
Git clone from FLAME repository.
::

    $ git clone *repository-address*

Pre-requisites (may need to apt-get with sudo)
::

    $apt-get install libboost-dev libboost-system-dev \
        libboost-thread-dev libboost-filesystem-dev \
        libboost-regex-dev libboost-program-options-dev \
        libboost-test-dev \
        build-essential cmake bison flex cppcheck git libhdf5-dev \
        python-numpy python-nose python3-numpy python3-nose

FLAME supports python 2.7 and 3.4, EPICS interface is optional.

Make ``build`` directory and compile with CMake.
::
    
    $ cd flame
    $ mkdir build
    $ cd build
    $ cmake ..
    $ make

Test FLAME (include the beam dynamics test).
::

    $ make test

Install with proper permissions.
::

    $ make install # may need to install with sudo
