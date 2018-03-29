from setuptools import setup
from setuptools.command.test import test as TestCommand
from fluff.config import FL_VERSION
import sys

class PyTest(TestCommand):
    def finalize_options(self):
        TestCommand.finalize_options(self)
        self.test_args = []
        self.test_suite = True

    def run_tests(self):
        #import here, cause outside the eggs aren't loaded
        import pytest
        errno = pytest.main(self.test_args)
        sys.exit(errno)

DESCRIPTION = "fluff : exploratory analysis and visualization of high-throughput sequencing data"

setup(name='biofluff',
      version=FL_VERSION,
      description=DESCRIPTION,
      author='Georgios Georgiou',
      author_email='g.georgiou@science.ru.nl',
      url = 'https://github.com/simonvh/fluff/',
      license='MIT',
      packages=[
          'fluff',
          'fluff/commands'
      ],
      entry_points={'console_scripts': ['fluff = fluff.parse:main'],},
      data_files=[],
      install_requires=["pysam",
                        "HTSeq",
                        "numpy",
                        "scipy",
                        "scikit-learn",
                        "matplotlib",
                        "palettable",
                        "pybedtools",
                        "pyBigWig",
                        ],
        tests_require=['pytest'],
        cmdclass = {'test': PyTest},
    )
