Developer's Guide
=================

Description of branches
-----------------------

main
____
This is the parent branch which
consolidates the current development in the repository.

release-vX.Y.Z
______________
These represent stable release branches.
Users should always check out the most recent stable release branch
when cloning the repository.

.. _dev-install-instructions:

How to incorporate updates to CATChem
-------------------------------------

In order to contribute code to CATChem, you will need to fork the
repository, make changes on your fork, and submit a pull request with your
changes.

(a) Fork the GitHub repository to your own GitHub account
    using the "Fork" button near the top right:

    https://github.com/ufs-community/CATChem

    .. note::
       You can pull updates from the main repository
       by using the "Fetch Upstream" button on your fork.
       Alternatively: [#clone]_ ::

          $ git remote add upstream git@github.com:ufs-community/CATChem.git
          $ git pull upstream main
          $ git push origin main

(b) Navigate on your working machine
    to where you would like to keep the CATChem code
    (e.g. in your work location) and clone [#clone]_ your fork::

       $ git clone git@github.com:$GitHubUsername/$ForkName.git

(c) Checkout the develop branch --- you need to do this with the remote branch
    as well as create a local tracking branch::

       $ git checkout origin/develop
       $ git checkout develop

(d) Make changes to your fork.

(e) Submit a pull request back to the main CATChem repository with your
    changes.

(f) Select two code reviewers (see list under development team section).

(g) Once those two reviewers approve the code, it can be merged into the develop branch.

.. _clone-notes:
.. [#clone] Note that in order to do an SSH clone,
   e.g. ::

      $ git clone git@github.com:ufs-community/CATChem.git

   you must have already
   `added an SSH key <https://docs.github.com/en/authentication/connecting-to-github-with-ssh/adding-a-new-ssh-key-to-your-github-account>`__
   to your GitHub account for your current machine.
   If you are new to GitHub, check out
   `this GitHub tutorial <https://jlord.us/git-it/>`__.
   We recommend the SSH method, but if you don't add an SSH key
   you can still clone the repositories via HTTPS, e.g. ::

       $ git clone https://github.com/ufs-community/CATChem.git

How to add a new *process*
--------------------------

CATChem is developed to be able to be easily extensible with new processes.
There are just a few steps that are required to be able to include a new process.

- First create a new directory under src/process for your new process::

    $ cd src/process
    $ mkdir src/process/<new process>

- Each process should include a process driver named ``CCPr_<new process>_Mod.F90``.
  You can find a template under ``src/process/Process_driver_template.F90``.
  In it, each process driver contains three phases:

  * Init: Processes the config and initializes process defaults if activated
  * Run: Runs the process and adds to ``DiagState`` and ``ChemState`` for any process
  * Finalize: Deallocate any arrays that were allocated.

- Each process should include a common module for any functions that may be used by schemes
  (sub-parameterizations or common calculations between different schemes in that process family).
  It also houses the process type and data information.

- Each process then can have one or more schemes within them.
  An example can be seen under ``src/process/dust`` and a template.
