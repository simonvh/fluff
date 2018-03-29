Installation
============

The straightforward way to install
----------------------------------

**Please note**: from version 3.0 on, fluff only works on Python 3.6+ and Python 2.x will no longer be supported.
For a Python 2 version you can install fluff 2.1.4.
Keep in mind that most scientific Python software will stop supporting Python 2 in 2020: python3_.

The most straightforward way to install fluff is with conda_ 
using the bioconda_ channel:

::

    $ conda config --add channels defaults
    $ conda config --add channels conda-forge
    $ conda config --add channels bioconda

    $ conda install biofluff

Or, in a seperate environment:

::

    $ conda create -n fluff python=3 biofluff
    # Before using fluff activate the environment:
    $ source activate fluff


.. _conda: https://docs.continuum.io/anaconda
.. _bioconda: https://bioconda.github.io/
.. _python3: https://python3statement.org/

Alternative: using pip
----------------------

You can use pip to install fluff, 
either as root user or in a `virtal environment
<http://docs.python-guide.org/en/latest/dev/virtualenvs/>`_ .

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
