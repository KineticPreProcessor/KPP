===============================
KPP-for-GEOS-Chem Documentation
===============================

.. raw:: html

   <p>
   <a href="https://github.com/geoschem/KPP/blob/GC_updates/LICENSE.txt"><img
   src="https://img.shields.io/badge/License-MIT-blue.svg"></a>
   <a href="https://doi.org/10.5281/zenodo.4552707"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.4552707.svg"></a>
   <a href="https://kpp.readthedocs.io/en/latest/"><img src="https://img.shields.io/readthedocs/kpp?label=ReadTheDocs"></a>
   </p>

:program:`KPP-for-GEOS-Chem` is `a clean implementation of
the Kinetic Pre Processor (KPP) <https://github.com/geoschem/kpp/tree/GC_updates>`__ that has
been customized for :program:`GEOS-Chem` v11-01 and later versions. You can use
:program:`KPP-for-GEOS-Chem` to create custom :program:`GEOS-Chem` chemistry mechanisms (or
to edit existing mechanisms).

.. toctree::
   :maxdepth: 20
   :caption: Basic Information

   basic-info/about.rst
   basic-info/requirements.rst
   basic-info/key_references.rst
   
.. toctree::
   :maxdepth: 20
   :caption: Usage Details

   usage-details/installation.rst
   usage-details/generating_f90_code.rst
   usage-details/use_custom_mech.rst

.. toctree::
   :maxdepth: 20
   :caption: Creating & Modifying Mechanisms

   creating-mechanisms/configuration_files.rst 
   creating-mechanisms/species.rst
   creating-mechanisms/reactions.rst
   creating-mechanisms/chemical_families.rst
   creating-mechanisms/solver_parameters.rst
  
.. toctree::
   :maxdepth: 20
   :caption: Help & Reference

   reference/known-bugs.rst 
   reference/SUPPORT.md
   reference/CONTRIBUTING.md
   geos-chem-shared-docs/editing_these_docs.rst
