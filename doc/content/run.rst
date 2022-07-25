.. _run_run:

Run mkCOInr
===========================================

The scripts were written and tested on Linux. MacOS users will probably be able to run them. I suggest `WSL <https://docs.microsoft.com/en-us/windows/wsl/>`_  for Windows users.

**To run any of the scripts, open a terminal and type**

.. code-block:: bash

	perl PATH_TO_SCRIPTS/name_of_the_script.pl -argument_name argument_value

For example if you are in the mkCOInr directory a command will look like this

.. code-block:: bash

	perl scripts/download_taxonomy.pl -outdir COInr_new/taxonomy


The list of arguments and options can be obtained by typing

.. code-block:: bash

	perl PATH_TO_SCRIPTS/name_of_the_script.pl -help


For example

.. code-block:: bash

	perl scripts/download_taxonomy.pl -help


