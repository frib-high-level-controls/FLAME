State Class
===========

.. py:class:: State(object)

    FLAME beam state class for Python API.

    .. py:function:: clone()

        Clone the beam state object.

        :return: :py:class:`State` object

.. _beamstate:

    - **Attributes - reference beam**

        .. list-table::
            :widths: 10 25

            * - :py:attr:`pos`
              - z position [m]
            * - :py:attr:`ref_beta`
              - Lorentz :math:`\beta` [1]
            * - :py:attr:`ref_bg`
              - Lorentz :math:`\beta \gamma` [1]
            * - :py:attr:`ref_gamma`
              - Lorentz :math:`\gamma` [1]
            * - :py:attr:`ref_IonEk`
              - Kinetic energy [eV/u]
            * - :py:attr:`ref_IonEs`
              - Nucleon mass [eV/u]
            * - :py:attr:`ref_IonQ`
              - Macro weight [1]
            * - :py:attr:`ref_IonW`
              - Total energy [eV/u]
            * - :py:attr:`ref_IonZ`
              - Charge to mass ratio [1]
            * - :py:attr:`ref_phis`
              - Absolute phase [rad]
            * - :py:attr:`ref_SampleIonK`
              - Phase speed [rad]
            * - :py:attr:`last_caviphi0`
              - Driven phase of the last rf cavity [deg]
            * - :py:attr:`transmat`
              - Transfer matrix of the last element

    - **Attributes - actual beam**

        .. list-table::
            :widths: 10 25

            * - :py:attr:`beta`
              - Lorentz :math:`\beta` [1]
            * - :py:attr:`bg`
              - Lorentz :math:`\beta \gamma` [1]
            * - :py:attr:`gamma`
              - Lorentz :math:`\gamma` [1]
            * - :py:attr:`IonEk`
              - Kinetic energy [eV/u]
            * - :py:attr:`IonEs`
              - Nucleon mass [eV/u]
            * - :py:attr:`IonQ`
              - Macro weight [1]
            * - :py:attr:`IonW`
              - Total energy [eV/u]
            * - :py:attr:`IonZ`
              - Charge to mass ratio [1]
            * - :py:attr:`phis`
              - Absolute phase [rad]
            * - :py:attr:`SampleIonK`
              - Phase speed [rad]
            * - :py:attr:`moment0`
              - Centroids of the all charge states.
            * - :py:attr:`moment0_env`
              - Weighted average of centroids for the all charge states.
            * - :py:attr:`moment0_rms`
              - Weighted average of rms size for the all charge states.
            * - :py:attr:`moment1`
              - Envelope matrixes of the all charge states.
            * - :py:attr:`moment1_env`
              - Weighted average of envelope matrixes for the all charge states.

    .. py:attribute:: pos

        **float**: z position of the reference beam. [m]

    .. py:attribute:: ref_beta

        **float**: Lorentz :math:`\beta` of the reference beam. [1]

    .. py:attribute:: ref_bg

        **float**: Lorentz :math:`\beta \gamma` of the reference beam. [1]

    .. py:attribute:: ref_gamma

        **float**: Lorentz :math:`\gamma` of the reference beam. [1]

    .. py:attribute:: ref_IonEk

        **float**: Kinetic energy of the reference beam. [eV/u]

    .. py:attribute:: ref_IonEs

        **float**: Nucleon mass of the reference beam. [eV/u]

    .. py:attribute:: ref_IonQ

        **float**: Macro weight of the reference beam. [1]

    .. py:attribute:: ref_IonW

        **float**: Total energy of the reference beam. [eV/u]

    .. py:attribute:: ref_IonZ

        **float**: Charge to mass ratio of the reference beam. [1]

    .. py:attribute:: ref_phis

        **float**: Absolute synchrotron phase of the reference beam. [rad]

    .. py:attribute:: ref_SampleIonK

        **float**: Phase speed of the reference beam. [rad]

    .. py:attribute:: last_caviphi0

        **float**: Driven phase of the last rf cavity. [deg]

    .. py:attribute:: transmat

        **list of matrix[7,7]**: Transfer matrix of the last element. This matrix is applied to moment0 and moment1 directly.


    .. py:attribute:: beta

        **list of float**: Lorentz :math:`\beta` of the all charge states. [1]

    .. py:attribute:: bg

        **list of float**: Lorentz :math:`\beta \gamma` of the all charge states. [1]

    .. py:attribute:: gamma

        **list of float**: Lorentz :math:`\gamma` of the all charge states. [1]

    .. py:attribute:: IonEk

        **list of float**: Kinetic energy of the all charge states. [eV/u]

    .. py:attribute:: IonEs

        **list of float**: Nucleon mass of the all charge states. [eV/u]

    .. py:attribute:: IonQ

        **list of float**: Macro weight of the all charge states. [1]

    .. py:attribute:: IonW

        **list of float**: Total energy of the all charge states. [eV/u]

    .. py:attribute:: IonZ

        **list of float**: Charge to mass ratio of the all charge states. [1]

    .. py:attribute:: phis

        **list of float**: Absolute synchrotron phase of the all charge states. [rad]

    .. py:attribute:: SampleIonK

        **list of float**: Phase speed of the all charge states. [rad]

    .. py:attribute:: moment0

        Centroids of the all charge states.

        **list of vector[7]**: :math:`[x, x', y, y', \phi, E_k, 1]` with [mm, rad, mm, rad, rad, MeV/u, 1].

    .. py:attribute:: moment0_env

        Weighted average of centroids for all charge states.

        **vector[7]**: :math:`[x, x', y, y', \phi, E_k, 1]` with [mm, rad, mm, rad, rad, MeV/u, 1].

    .. py:attribute:: moment0_rms

        Weighted average of rms beam envelopes (2nd order moments) for the all charge states.

        **vector[7]**: rms of :math:`[x, x', y, y', \phi, E_k, 1]` with [mm, rad, mm, rad, rad, MeV/u, 1].

    .. py:attribute:: moment1

        Envelope matrixes of the all charge states.

        **list of matrix[7,7]**:

        Cartisan product of :math:`[x, x', y, y', \phi, E_k, 1]^2` with [mm, rad, mm, rad, rad, MeV/u, 1] :math:`^2`.

    .. py:attribute:: moment1_env

        Weighted average of envelope matrixes for the all charge states.

        **matrix[7,7]**:

        Cartisan product of :math:`[x, x', y, y', \phi, E_k, 1]^2` with [mm, rad, mm, rad, rad, MeV/u, 1] :math:`^2`.

