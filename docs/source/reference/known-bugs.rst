#####################
Known bugs and issues
#####################

Please see our `Issue tracker on GitHub
<https://github.com/KineticPreProcessor/KPP/issues>`_ for a list of recent
bugs and fixes.

===================
Current bug reports
===================

These `bug reports (on GitHub)
<https://github.com/KineticPreProcessor/KPP/issues?q=is%3Aissue%20state%3Aopen%20label%3Abug>`_
are currently unresolved. We hope to fix these in future releases.

LSODE integrator is not thread-safe
-----------------------------------

We have discovered that the current implementation of the LSODE
integrator is not thread-safe for `OpenMP parallelization
<https://www.openmp.org/>`_.  When LSODE is called from within an
OpenMP parallel loop, the integration will fail because key internal
variables in LSODE will be overwritten by concurrent threads.

============================
Bugs that have been resolved
============================

These `bugs (reported on GitHub) <https://github.com/KineticPreProcessor/KPP/issues?q=is%3Aissue%20state%3Aclosed%20label%3Abug>`_ have been resolved.
