.. _check-gchp:

###############################
Check GCHP emission diagnostics
###############################

Because `GCHP <https://gchp.readthedocs.io>`_ uses the MAPL library's
History component for diagnostic archiving, all emission entries must
be entered individually in both the `HISTORY.rc
<https://gchp.readthedocs.io/en/stable/user-guide/config-files/HISTORY_rc.html>`_
and `HEMCO_Diagn.rc
<https://gchp.readthedocs.io/en/stable/geos-chem-shared-docs/doc/hemco-diagn.html>`_
template files.  These template files are located in the
:file:`run/GCHP` folder of the GEOS-Chem "science codebase"
repository.

This example demonstrates how you can cross-check the GCHP emission
diagnostic entries.  This is especially important if you are adding
new emissions diagnostics (so that you won't forget to update one file
or the other).

.. _check_gchp_code:

===========
Source code
===========

**Script location:** `gcpy/examples/diagnostics/check_gchp_emission_diags.py
<https://github.com/geoschem/gcpy/blob/main/gcpy/examples/diagnostics/check_gchp_emission_diags.py>`_

.. _check-gchp-usage:

=====
Usage
=====

First, clone the GCHP source code (if you haven't done so already).
Then navigate to the :file:`run/GCHP` folder, which is the top-level
folder for GCHP run directory template files.

.. code-block:: console

   $ git clone --recurse-submodules https://github.com/geoschem/GCHP

   $ cd gchp/run

Activate your GCPy environment with mamba or conda:

.. code-block:: console

   $ mamba activate gcpy_env    # If using mamba, or

   $ conda activate gcpy_env    # If using conda

Use one of these commands to check GCHP emissions diagnostics for a
given simulation type:

.. code-block:: console

   $ python -m gcpy.examples.diagnostics.check_gchp_emission_diags . carbon

   $ python -m gcpy.examples.diagnostics.check_gchp_emission_diags . fullchem

   $ python -m gcpy.examples.diagnostics.check_gchp_emission_diags . TransportTracers

   $ python -m gcpy.examples.diagnostics.check_gchp_emission_diags . tagO3

You will then see output similar to this (fullchem example shown).

.. code-block:: console

   ===============================================================================
   Common to both HISTORY.rc and HEMCO_Diagn.rc
   ===============================================================================
     #EmisCH4_Anthro
     #EmisCH4_BioBurn
     #EmisCH4_Ship
     #EmisCH4_Total
     #EmisSESQ_Biogenic
     #InvAEIC_ACET
     #InvAEIC_ALD2
     #InvAEIC_ALK4
     #InvAEIC_BCPI
     #InvAEIC_C2H6
     #InvAEIC_C3H8
     ... etc ...
     EmisACET_BioBurn
     EmisACET_Biogenic
     EmisACET_Ocean
     EmisACET_Total
     EmisACR_BioBurn
     EmisACR_Total
     EmisACTA_BioBurn
     EmisACTA_Total
     EmisALD2_Anthro
     ... etc ...

   ===============================================================================
   In HISTORY.rc but not in HEMCO_Diagn.rc
   ===============================================================================

   ===============================================================================
   In HEMCO_Diagn.rc but not in HISTORY.rc
   ===============================================================================
     #EmisNO_Fert
     #InvCEDS_ALK6
     #InvGTChlorine_HCl

The output indicates that some diagnostics in :file:`HEMCO_Diagn.rc`
are not present in :literal:`HISTORY.rc`, and should be added there
for consistency.  A comment before the emission entry means that the
diagnostic will be disabled by default.
