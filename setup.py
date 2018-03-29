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

def readme():
    with open('README.md') as f:
        return f.read()

DESCRIPTION = "fluff : exploratory analysis and visualization of high-throughput sequencing data"

setup(name='biofluff',
      version=FL_VERSION,
      description=DESCRIPTION,
      long_description=readme(),
      long_description_content_type='text/markdown',
      author='Simon van Heeringen',
      author_email='simon.vanheeringen@gmail.com',
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
      setup_requires=['setuptools>=38.6.0'], 
      tests_require=['pytest'],
      cmdclass = {'test': PyTest},
	  classifiers=[
            'Development Status :: 4 - Beta',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: MIT License',
            'Operating System :: MacOS :: MacOS X',
            'Operating System :: POSIX :: Linux',
            'Programming Language :: Python :: 3.6',
            'Programming Language :: Python :: 3.7',
			'Programming Language :: Python :: 3 :: Only',
            'Topic :: Scientific/Engineering :: Bio-Informatics',
            ],

    )
