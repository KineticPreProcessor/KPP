.. |br| raw:: html

   <br />

#############################
The Kinetic PreProcessor: KPP
#############################

| **An Environment for the Simuation of Chemical Kinetic Systems**

| **Adrian Sandu**\ :sup:`1`\ **, Rolf Sander**\ :sup:`2`\ **, Michael S. Long**\ :sup:`3`\ **, Haipeng Lin**\ :sup:`4`\ **,**
| **Robert M. Yantosca**\ :sup:`4`\ **, Lucas Estrada**\ :sup:`4`\ **, Lu Shen**\ :sup:`5`\ **, and Daniel J. Jacob**\ :sup:`4`

| :sup:`1` *Virginia Polytechnic Institute and State University, Blacksburg, VA, USA*
| :sup:`2` *Max-Planck Institute for Chemistry, Mainz, Germany*
| :sup:`3` *Renaissance Fiber, LLC, North Carolina, USA*
| :sup:`4` *Harvard John A. Paulson School of Engineering and Applied Sciences, Cambridge, MA, USA*
| :sup:`5` *School of Physics, Peking University, Beijing, China*

----

.. |latest-stable-release| image:: https://img.shields.io/github/v/release/KineticPreProcessor/KPP?label=Latest%20Stable%20Release
   :target: https://github.com/KineticPreProcessor/KPP/releases/
   :alt: Latest Stable Release

.. |release-date| image:: https://img.shields.io/github/release-date/KineticPreProcessor/KPP
   :target: https://github.com/KineticPreProcessor/KPP/releases/
   :alt: Release Date

.. |license| image:: https://img.shields.io/badge/License-GPLv3-blue.svg
   :target: https://github.com/KineticPreProcessor/KPP/blob/main/LICENSE.txt
   :alt: License: GPLv3

.. |doi| image:: https://zenodo.org/badge/DOI/10.5281/zenodo.6828025.svg
   :target: https://doi.org/10.5281/zenodo.6828025
   :alt: DOI

.. |readthedocs| image:: https://readthedocs.org/projects/kpp/badge/?version=stable
   :target: https://kpp.readthedocs.io/en/latest/
   :alt: ReadTheDocs

.. |ci-tests| image:: https://github.com/KineticPreProcessor/KPP/actions/workflows/run-ci-tests.yml/badge.svg
   :target: https://github.com/KineticPreProcessor/KPP/actions/workflows/run-ci-tests.yml
   :alt: C-I tests

|latest-stable-release| |release-date| |license| |br|
|doi| |readthedocs| |ci-tests|

----

This site provides instructions for :program:`KPP`, the Kinetic PreProcessor.

Contributions (e.g., suggestions, edits, revisions) would be greatly
appreciated. See :ref:`editing_this_user_guide` and our contributing
guidelines.  If you find something hard to understand---let us know!

.. toctree::
   :caption: Getting Started
   :maxdepth: 2

   getting_started/revision_history.rst
   getting_started/installation.rst
   getting_started/running_kpp_sample_mech.rst

.. toctree::
   :caption: Input for KPP
   :maxdepth: 3

   input/input_overview.rst
   input/kpp_sections.rst
   input/kpp_commands.rst
   input/inlined_code.rst
   input/auxiliary_files.rst
   input/filter-passive-spc.rst

.. toctree::
   :caption: Output from KPP
   :maxdepth: 3

   output/code_f90.rst
   output/code_c.rst
   output/code_matlab.rst
   output/makefile.rst
   output/log_file.rst

.. toctree::
   :caption: Numerical Integration Methods
   :maxdepth: 3

   num_methods/rosenbrock-methods.rst
   num_methods/runge_kutta_methods.rst
   num_methods/backward_diff.rst
   num_methods/forward_diff.rst
   num_methods/controlling_the_integrator.rst
   num_methods/output_from_integrators.rst

.. toctree::
   :caption: Technical information
   :maxdepth: 3

   tech_info/kpp_dir_structure.rst
   tech_info/kpp_env_vars.rst
   tech_info/kpp_internal_modules.rst
   tech_info/adding_new_commands.rst
   tech_info/about_ci_tests.rst
   tech_info/bnf_description_of_kpp_lang.rst

.. toctree::
   :caption: KPP Reference
   :maxdepth: 2

   citations/acknowledgments.rst
   citations/kpp_references.rst

.. toctree::
   :caption: Help and Support
   :maxdepth: 2

   reference/known-bugs.rst
   reference/support.rst
   reference/contributing.rst
   reference/editing_these_docs.rst
