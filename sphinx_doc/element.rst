Lattice File
============

General parameter
-----------------

Basic format of the general parameters are,

.. code-block:: none

    keyword1 = "value1";
    keyword2 = "value2";
    ...

.. _genpara:

.. list-table::
    :header-rows: 1
    :widths: 2, 2, 6

    * - | keyword
      - | value
      - | description
    * - | **sim_type**
      - | "MomentMatrix"
      - | Simulation mode. FRIB simulation uses the particular
        | mode "MomentMatrix".
    * - | **MpoleLevel**
      - | "0", "1", or "2"
      - | Multipole term controller for the rf cavities.
        | "0" - only include focusing and defocusing effect
        | "1" - include dipole terms
        | "2" - include dipole and quadrupole terms
    * - | **EmitGrowth**
      - | "0" or "1"
      - | Flag for cross-cavity emittance growth effect.
        | "0" - False (no emittance growth)
        | "1" - True (calculate emittance growth)
    * - | **HdipoleFitMode**
      - | "0" or "1"
      - | Flag for auto-adjustment of bending element
        | "0" - use "bg" or "beta" for the bending strength
        | "1" - auto-adjust the bending strength


Beam parameter
--------------

Basic format of the beam parameters are,

.. code-block:: none

    keyword1 = value1;
    keyword2 = [value2, value3]; # list input
    ...

.. _beampara:

.. list-table::
    :header-rows: 1
    :widths: 3, 2 , 6

    * - | keyword
      - | value
      - | description
    * - | **IonEs**
      - | float
      - | Nucleaon mass of the reference beam. [eV/u]
    * - | **IonEk**
      - | float
      - | Initial kinetic energy of the reference beam. [eV/u]
    * - | **IonChargeStates**
      - | list of float
      - | List of charge to mass ratios of the all charge states. [1]
    * - | **NCharge**
      - | list of float
      - | List of macro weights of the all charge states. [1]
    * - | **${vector_variable}${n}**
      - | vector[7]
      - | Initial centroid vector of the **n**-th charge state.
        | *${vector_variable}* is defined in :cpp:type:`source`.
        | :math:`[x, x', y, y', \phi, E_k, 1]` with
        | [mm, rad, mm, rad, rad, MeV/u, 1]
    * - | **${matrix_variable}${n}**
      - | vector[49]
      - | Flattened initial envelope matrix of the **n**-th charge state.
        | *${matrix_variable}* is defined in :cpp:type:`source`.
        | Cartisan product of :math:`[x, x', y, y', \phi, E_k, 1]^2` with
        | [mm, rad, mm, rad, rad, MeV/u, 1] :math:`^2`
    * - | **Eng_Data_Dir**
      - | string
      - | Directory path of the rf cavity data.
        | ``dir(path)`` supports relative path.


Lattice elements
----------------

Basic format of the one lattice element is,

.. code-block:: none

    name_of_element1: element_type, parameter1 = value1, parameter2 = value2, ... ;

After writing down the all lattice elements, user need to specify the lattice cell and the cell to USE.

.. code-block:: none

    # define the cell
    name_of_cell: LINE = ( name_of_element1, name_of_element2, name_of_element3, ... );

    # set the cell to USE
    USE: name_of_cell;

.. list-table::
    :header-rows: 1
    :widths: 3 , 6

    * - **element_type**
      - description
    * - :cpp:type:`source`
      - Starting point of the simulation.
    * - :cpp:type:`marker`
      - Marker element.
    * - :cpp:type:`stripper`
      - Chage stripper element.
    * - :cpp:type:`tmatrix`
      - User input transfer matrix.
    * - :cpp:type:`orbtrim`
      - Orbit trim element.
    * - :cpp:type:`drift`
      - Drift space element.
    * - :cpp:type:`solenoid`
      - Solenoid magnet element.
    * - :cpp:type:`quadrupole`
      - Magnetic quadrupole element.
    * - :cpp:type:`equad`
      - Electrostatic quadrupole element.
    * - :cpp:type:`sbend`
      - Magnetic bend element.
    * - :cpp:type:`edipole`
      - Electrostatic dipole element.
    * - :cpp:type:`rfcavity`
      - RF cavity element.

Special element
^^^^^^^^^^^^^^^

.. cpp:type:: source

    Starting point of the simulation. Initial beam state parameters are set at this element.

    :parameters: **vector_variable**: string

                    | Name key of the initial centroid vector.

                 **matrix_variable**: string

                    | Name key of the initial envelope matrix.

.. cpp:type:: marker

    Marker element. Nothing to do.

.. cpp:type:: stripper

    Stripper element.

    :parameters: **IonChargeStates**: list of float (optional)

                    | List of charge to mass ratios after the charge stripper. [1]

                 **NCharge**: list of float (optional)

                    | List of macro weights after the charge stripper. [1]

                 **Stripper_E1Para**: float (optional)

                    | Constant part of the energy struggling parameter of the charge stripper. [MeV/u]

                 **Stripper_lambda**: float (optional)

                    | Momentum spread factor :math:`\lambda` of the charge stripper. [1]

                 **Stripper_upara**: float (optional)

                    | Momentum spread factor :math:`U` of the charge stripper. [1]

                 **Stripper_E0Para**: vector[3]

                    | Energy loss parameters due to the ionization.
                    | [Constant_part, Energy_dependence, Thickness_depenedence] with [eV/u, 1, 1]

                 **Stripper_Para**: vector[3]

                    | Stripper foil parameters.
                    | [Thickness, Thickness_variation, reference_energy] with [um, %, eV/u]

.. cpp:type:: tmatrix

    User input transfer matrix element.

    :parameter: **matrix**: vector[49]

                    | Flattened :math:`7 \times 7` transfer matrix.


Optical element
^^^^^^^^^^^^^^^

.. cpp:type:: orbtrim

    Orbit trim element. This can be use as steering magnet.

    :parameters: **realpara**: int

                    | Realistic input parameter flag for the beam kick angle.
                    | **0** - use ``theta_x`` and ``theta_y`` for the beam kick.
                    | **1** - use ``tm_xkick`` and ``tm_ykick`` for the beam kick.

                 **theta_x**: float

                    | Horizontal beam kick angle. [rad]

                 **theta_y**: float

                    | Vertical beam kick angle. [rad]

                 **tm_xkick**: float

                    | Magnetic field strength for the horizontal beam kick. [T*m]


                 **tm_xkick**: float

                    | Magnetic field strength for the vertical beam kick. [T*m]

                 **xyrotate**: float

                    | Transverse rotation angle of the beam. [deg]

    .. Note::

        In the case of user puts both "beam kick information" and "transverse rotation angle" to the ONE orbtrim element,
        the process order is, beam kick -> transverse rotation. In other words, the beam kick is effected BEFORE the transverse rotation.

.. cpp:type:: drift

    Drift space element.

    :parameters: **L**: float

                    | Length of the lattice element. [m]

.. cpp:type:: solenoid

    Solenoid magnet element.

    :parameters: **L**: float

                    | Length of the lattice element. [m]

                 **B**: float

                    | Solenoid strength (:math:`B_z`). [T]

                 **dx**: float (default: 0.0)

                    | Misalignment of horizontal shift. [m]

                 **dy**: float (default: 0.0)

                    | Misalignment of vertical shift. [m]

                 **pitch**: float (default: 0.0)

                    | Misaglignment of pitch angle. [rad]

                 **yaw**: float (default: 0.0)

                    | Misaglignment of yaw angle. [rad]

                 **roll**: float (default: 0.0)

                    | Misaglignment of roll angle. [rad]

.. cpp:type:: quadrupole

    Magnetic quadrupole element.

    :parameters: **L**: float

                    | Length of the lattice element. [m]

                 **B2**: float

                    | Quadrupole field gradient. [T/m]
                    | Positive value means horizontal focusing.

                 **dx**, **dy**, **pitch**, **yaw**, **roll**: float

                    | Misalignment parameters. See :cpp:type:`solenoid` case.

.. cpp:type:: equad

    Electrostatic quadrupole element.

    :parameters: **L**: float

                    | Length of the lattice element. [m]

                 **V**: float

                    | Electrostatic quadrupole pole tip voltage. [V]
                    | Positive value means horizontal focusing.

                 **radius**: float

                    | Electrostatic quadrupole pole tip radius. [m]

                 **dx**, **dy**, **pitch**, **yaw**, **roll**: float

                    | Misalignment parameters. See :cpp:type:`solenoid` case.

.. cpp:type:: sbend

    Magnetic bend element.

    :parameters: **L**: float

                    | Length of the lattice element. [m]

                 **phi**: float

                    | Bend angle. [deg]

                 **phi1**: float

                    | Front pole face angle. [deg]

                 **phi2**: float

                    | Back pole face angle. [deg]

                 **bg**: float (optional: Used in the case of :ref:`"HdipoleFitMode" <genpara>` is **0**.)

                    | Lorentz :math:`\beta \gamma` for the reference beam. [1]
                    | This parameter is correspond to the bend field strength.

                 **dx**, **dy**, **pitch**, **yaw**, **roll**: float

                    | Misalignment parameters. See :cpp:type:`solenoid` case.

.. cpp:type:: edipole

    Electrostatic dipole (bend) element.

    :parameters: **L**: float

                    | Length of the lattice element. [m]

                 **phi**: float

                    | Bend angle. [deg]

                 **beta**: float (optional: Used in the case of :ref:`"HdipoleFitMode" <genpara>` is **0**.)

                    | Lorentz :math:`\beta` for the reference beam. [1]
                    | This parameter is correspond to the bend field strength.

                 **fringe_x**: float

                    | Horizontal fringe term. [rad/mm]

                 **fringe_y**: float

                    | Vertical fringe term. [rad/mm]

                 **asymfac**: float

                    | Characteristic parameter of the kinetic energy change
                      due to the middle point potential deviation from ground. [1]

                 **spher**: int

                    | Flag for the electrostatic dipole shape.
                    | **0** - cylindrical electrostatic dipole
                    | **1** - spherical electrostatic dipole

                 **ver**: int

                    | Flag for the bending direction.
                    | **0** - horizontal bend
                    | **1** - vertical bend

                 **dx**, **dy**, **pitch**, **yaw**, **roll**: float

                    | Misalignment parameters. See :cpp:type:`solenoid` case.

.. cpp:type:: rfcavity

    RF cavity element.

    :parameters: **L**: float

                    | Length of the lattice element. [m]

                 **cavtype**: string

                    | Cavity type. Supports "Generic", "0.041QWR", "0.085QWR", "0.29HWR", and "0.53HWR".
                      :ref:`The file format is described here. <cavformat>`

                 **f**: float

                    | RF frequency of the cavity. [Hz]

                 **phi**: float

                    | Input phase of the cavity. [deg]

                 **syncflag**: int

                    | Flag for synchronous phase input (**phi**).
                    |    **0** for driven phase input. 
                    |    **1** for synchronous phase input. (default)

                 **scl_fac**: float

                    | Scaling factor of the field. [1]

                 **datafile**: string (optional: Used in the case of ``cavtype`` = "Generic")

                    | File path of the rf cavity data.

                 **Rm**: float (optional: Used in the case of ``cavtype`` = "Generic")

                    | Characteristic radial length of the multipole expansion. [mm]

                 **dx**, **dy**, **pitch**, **yaw**, **roll**: float

                    | Misalignment parameters. See :cpp:type:`solenoid` case.


.. _cavformat:

Rf cavity data format
---------------------

FLAME using Thin-Lens-Model for rf cavity calculation.
Rf cavity data is composed of "Longitudinal axis data", "Multipole lattice data", "Multipole field data", and "TTF fitting data".

Hard-coded FRIB cavity models
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For typical rf cavity in FRIB, the "TTF fitting data" is hard-coded in FLAME.
Following files are required for each rf cavity type.

.. list-table::
    :header-rows: 1

    * - **cavtype**
      - **Longitudinal axis data**
      - **Multipole lattice data**
      - **Multipole field data**
    * - "0.041QWR"
      - "axisData_41.txt"
      - "Multipole41/thinlenlon_41.txt"
      - "Multipole41/CaviMlp_41.txt"
    * - "0.085QWR"
      - "axisData_85.txt"
      - "Multipole85/thinlenlon_85.txt"
      - "Multipole85/CaviMlp_85.txt"
    * - "0.29HWR"
      - "axisData_29.txt"
      - "Multipole29/thinlenlon_29.txt"
      - "Multipole29/CaviMlp_29.txt"
    * - "0.53HWR"
      - "axisData_53.txt"
      - "Multipole53/thinlenlon_53.txt"
      - "Multipole53/CaviMlp_53.txt"

Generic rf cavity model
^^^^^^^^^^^^^^^^^^^^^^^

FLAME supports *lattice format* input for the generic rf cavity model.

The basic format of the rf cavity data is similar to the main lattice file,

.. code-block:: none

    Rm = value1;

    Ez = [
    z1, Ez1,
    z2, Ez2,
    z3, Ez3,
    ...
    ];

    name_of_element1: element_type, parameter1 = value1, parameter2 = value2, ... ;
    ...

    cell: LINE =(name_of_element1, ...);
    USE: cell;

.. list-table::
    :header-rows: 1
    :widths: 2, 2, 6

    * - | keyword
      - | value
      - | description
    * - | **Rm**
      - | float
      - | Characteristic radial length of the multipole expansion. [mm]
    * - | **Ez**
      - | vector[2*n]
      - | On axis :math:`E_z` data.
        | The odd index (1,3,5,...) is z position. [mm]
        | The even index (2,4,6,...) is Electric field strength. [V/m]

Lattice element for the rf cavity data
""""""""""""""""""""""""""""""""""""""""""

Drift space is the same format as the main lattice but unit of ``L`` is [mm] - :cpp:type:`drift`

.. cpp:type:: EDipole

    Dipole term generated by the electric field.

    :parameters: **L**: float

                    | Length of the lattice element. [mm]
                    | This parameter should be 0.0 in Thin-Lens-Model.

                 **V0**: float

                    | Amplitude of the multipole term. [MV]

                 **attr**: vector[20]

                    | TTF fitting parameter. :ref:`(see here) <ttfnote>`
                    | 1 to 10 - fitting parameter for :math:`T`
                    | 11 to 20 - fitting parameter for :math:`S`

.. cpp:type:: EFocus

    Constant focusing term generated by the electric field.

    Parameters are the same as :cpp:type:`EDipole`.

.. cpp:type:: EQuad

    Quadrupole term generated by the electric field.

    Parameters are the same as :cpp:type:`EDipole`.


.. cpp:type:: HMono

    Dipole term generated by the magnetic field.

    :parameters: **L**: float

                    | Length of the lattice element. [mm]
                    | This parameter should be 0.0 in Thin-Lens-Model.

                 **V0**: float

                    | Amplitude of the multipole term. [MA]

                 **attr**: vector[20]

                    | TTF fitting parameter. :ref:`(see here) <ttfnote>`
                    | 1 to 10 - fitting parameter for :math:`T`
                    | 11 to 20 - fitting parameter for :math:`S`

.. cpp:type:: HFocus

    Constant focusing term generated by the magnetic field.

    Parameters are the same as :cpp:type:`HMono`.

.. cpp:type:: HQuad

    Quadrupole term generated by the magnetic field.

    Parameters are the same as :cpp:type:`HMono`.

.. cpp:type:: AccGap

    Acceleration gap term by the longitudinal electric field.

    :parameters: **L**: float

            | Length of the lattice element. [mm]
            | This parameter should be 0.0 in Thin-Lens-Model.

         **V0**: float

            | Amplitude of the multipole term. [MV]

         **attr**: vector[23]

            | TTF fitting parameter. :ref:`(see here) <ttfnote>`
            | 1 to 10 - fitting parameter for :math:`T`
            | 11 to 20 - fitting parameter for :math:`S`
            | 21 to 23 - fitting parameter for the synchronous phase


.. _ttfnote:

.. Note::

    FLAME is using TTF-calculation acceleration technique to boost cavity modeling speed.
    TTF factor :math:`T` and :math:`S` are pre-calculated and fitted using 9th order polynomial function according
    to different particle phase speed :math:`k`. :math:`n`-th fitting parameter :math:`p_n` is listed as,

    .. math::

        T(k), S(k) = \sum^9_{n=0} p_n k^{(9-n)}.

    The driven-phase calculation is also boosted by similar method.
    The phase transferring factor :math:`\varphi_c` is fitted by :math:`p_{i = 0, 1, 2}`

    .. math::

        \varphi_c = p_0 E^{p_1} + p_2.

    Here, :math:`E` is the kinetic energy. The driven phase :math:`\varphi_d` is calculated by using :math:`\varphi_c`,

    .. math::

        \varphi_d = \varphi_s - \varphi_c - m \varphi_\text{abs}

    where, :math:`\varphi_s` is the synchronous phase in input, :math:`\varphi_\text{abs}`
    is absolute phase in front of the rf cavity, and :math:`m` is the harmonic number.