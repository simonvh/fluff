from distutils.core import setup

from fluff.config import FL_VERSION

DESCRIPTION = """
fluff - plots and graphs 
"""

setup(name='fluff',
      version=FL_VERSION,
      description=DESCRIPTION,
      author='Georgios Georgiou',
      author_email='g.georgiou@science.ru.nl',
      license='MIT',
      packages=[
          'fluff',
          'fluff/commands'
      ],
      scripts=[
          "scripts/fluff_bandplot.py",
          "scripts/fluff_profile.py",
          "scripts/fluff_heatmap.py",
          "scripts/fluff",
      ],
      data_files=[],
      install_requires=[
        "pysam",
        "pybedtools",
        "HTSeq",
        "matplotlib >= 1.1.1",
        "numpy >= 1.6.0",
        "scipy >= 0.9.0",
        "colorbrewer",
        "Pycluster >= 1.52",
        "pp >= 1.6",
      ],
      )
