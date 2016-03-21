Installation
============

The straightforward way to install
----------------------------------

The most straightforward way to install fluff is by using `conda
<https://docs.continuum.io/anaconda>`_.  
*Note: we are working on creating a conda package, this option 
is not functional at this moment.*


::

    $ conda install biofluff -c bioconda


Alternative: using pip
----------------------

You can use pip to install fluff, 
either as root user or in a `virtal environment
<http://docs.python-guide.org/en/latest/dev/virtualenvs/>`_.

:: 

    $ pip install biofluff


Prerequisites for installation on Mac OS X
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For installation on Mac OS X you might need some additional items:

- Xcode (free on mac app store)
- Homebrew_
- pip_
- gfortran_
- bedtools2_
- Cython_

.. _Homebrew: http://brew.sh
.. _pip: http://pip.readthedocs.org/en/stable/installing/
.. _gfortran: https://cran.r-project.org/bin/macosx/tools/
.. _bedtools2: https://github.com/arq5x/bedtools2
.. _Cython: http://cython.org/

Installation from source
------------------------

You can check out the development version of fluff using git:

::

    # option 1
    $ git clone https://github.com/simonvh/fluff.git
    $ cd fluff

Alternatively, you can download the lastest version of fluff at:

https://github.com/simonvh/fluff/releases

In this case, start by unpacking the source archive

::

  # option 2
  $ tar xvzf fluff-<version>.tar.gz
  $ cd fluff-<version>

Now you can build fluff with the following command:

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
