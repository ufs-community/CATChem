Developer's Guide
=================

Description of Branches
-----------------------

There are several primary repository branches
for specific development tracks.

main
____
This is the stable release branch.
It is updated from the *develop* branch
and tagged prior to each release.

develop
_______
This is the parent branch in which
to consolidate the various development tracks.
General practice, when working within the central repository,
is to create a working branch, branched from *develop*,
with a naming convention such as *develop_topic*.
The *develop_topic* branches are merged back to *develop*
via pull requests.
Forks of the central repository should follow the same
branching conventions.

.. _dev-install-instructions:

How to incorporate updates to CATChem
--------------------------------------------

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

