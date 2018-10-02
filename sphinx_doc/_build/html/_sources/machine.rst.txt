Machine class
=============

.. py:class:: Machine(config)

    FLAME Machine class for Python API.

    :parameter: **config**: dict, list of tuples, or byte buffer

                   | Input lattice data.

    .. py:function:: conf(index=None)

        Check configuration of the Machine object.

        :parameter: **index**: int (optional)

                        | Index of the lattice element.

        :returns: dict

                    | Configuration of the lattice element

        .. Note::

            In the case of ``index`` is *None*, :py:func:`conf` returns *initial* configuration of the lattice.


    .. py:function:: allocState(config=None)

        Allocate the beam state object.

        :parameter: **config** : dict

                        | Input lattice data. Empty dict is required as dummy data.

        :returns: :py:class:`State` object

                      | Beam state object (see here)

    .. py:function:: propagate(state, start=0, max=INT_MAX, observe=None)

        Run envelope tracking simulation.

        :parameters: **state**: :py:class:`State` object

                        | Allocated beam state object

                    **start**: int (optional)

                        | Index of the starting lattice element.

                    **max**: int (optional)

                        | Number of elements to advance. Negative value works as backward propagation.
                          (E.g. start = 5 and max = 10 mean propagate from 5th element to 14th element.)

                    **observe**: list of int (optional)

                        | List of indexes for observing the beam state.

        :returns: list

                    | List of the beam states at ``observe`` points. Each tuple has (*index*, *State*).

    .. py:function:: reconfigure(index, config)

            Reconfigure the lattice element configuration.

            :parameters: **index**: int

                            | Index of the lattice element.


                         **config**: dict

                            | New configuration of the lattice element parameter.

    .. py:function:: find(name=None, type=None)

            Find the indexes of the lattice elements by *name* or *type*.

            :parameter: **name**: str or unicode

                            | Name of the lattice element to find.


                        **type**: str or unicode

                            | Type of the lattice element to find.

            :returns: list

                        | List of matched element indexes.

