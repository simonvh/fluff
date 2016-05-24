from setuptools import setup

from fluff.config import FL_VERSION

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
      scripts=[
          "scripts/fluff_bandplot.py",
          "scripts/fluff_profile.py",
          "scripts/fluff_heatmap.py",
          "scripts/fluff",
      ],
      data_files=[],
      install_requires=["pysam",
                        "HTSeq",
                        "numpy",
                        "scipy",
                        "matplotlib",
                        "colorbrewer",
                        "pybedtools",
                        "Pycluster"
                       ],
    #   dependency_links = [
    #     "http://bonsai.hgc.jp/~mdehoon/software/cluster/Pycluster-1.52.tar.gz",
    # ],
      )
