.. _kpp-env-vars:

#########################
KPP environment variables
#########################

In order for KPP to find its components, it has to know the path to the
location where the KPP distribution is installed. This is achieved by
setting the :envvar:`$KPP_HOME` environment variable to the path where
KPP is installed.

The :file:`$KPP_HOME/bin` directory. should be added to the
:envvar:`PATH` variable.

There are also several optional environment variable that control the places
where KPP looks for module files, integrators, and drivers:

========
KPP_HOME
========

Required, stores the absolute path to the KPP distribution.

Default setting: none.

================
KPP_FLEX_LIB_DIR
================

Optional. Use this to specify the path to the :ref:`flex library
file <installation-flex>` (:file:`libfl.so` or :file:`libfl.a`) that are
needed to :ref:`build the KPP executable <installation-build>`. The KPP
build sequence will use the path contained in
:envvar:`KPP_FLEX_LIB_DIR` if the flex library file cannot be found
in  :file:`/usr/lib`, :file:`/usr/lib64`, and similar standard
library paths.

=========
KPP_MODEL
=========

Optional, specifies additional places where KPP will look for model
files before searching the default location.

Default setting: :file:`$KPP_HOME/models`.

=======
KPP_INT
=======

Optional, specifies additional places where KPP will look for
integrator files before searching the default.

Default setting: :file:`$KPP_HOME/int`.

=======
KPP_DRV
=======

Optional specifies additional places where KPP will look for driver
files before searching the default directory.

Default setting: :file:`$KPP_HOME/drv`.
