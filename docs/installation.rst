Installation
============

Installation packages
---------------------
::

Installation on Ubuntu or Debian
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Installation on Mac OSX
~~~~~~~~~~~~~~~~~~~~~~~

Additional items
-  Xcode

::

  free on mac app store

-  Homebrew

::

  ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"

-  pip

::

  http://pip.readthedocs.org/en/stable/installing/

-  gfortran

::

  https://cran.r-project.org/bin/macosx/tools/

-  bedtools2

::

  https://github.com/arq5x/bedtools2

-  Cython

::

  http://cython.org/


Installation from source
------------------------
::

You can download the lastest version of fluff at:

https://github.com/simonvh/fluff/releases

Start by unpacking the source archive

::

  tar xvzf fluff-<version>.tar.gz
  cd fluff-<version>

You can build fluff with the following command:

::

  python setup.py build


If you encounter no errors, go ahead with installing fluff:

- root privileges required

::

  sudo python setup.py install


- install in user site-package

::

  sudo python setup.py install --user



Using Conda
-----------

::

You can download the lastest version of fluff at:

https://github.com/simonvh/fluff/releases

Start by unpacking the source archive

::

  tar xvzf fluff-<version>.tar.gz
  cd fluff-<version>


Make a copy of the environment:

::

  conda env create -f conda/environment.yml


**Activate the new environment from Linux, OS X:**

::

  source activate fluff

**Activate the new environment from Windows:**

::

  activate fluff