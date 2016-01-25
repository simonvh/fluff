from pip.req import parse_requirements
from setuptools import setup

from fluff.config import FL_VERSION

DESCRIPTION = """
fluff - plots and graphs 
"""

install_reqs = parse_requirements('requirements.txt', session=False)
reqs = [str(ir.req) for ir in install_reqs]

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
      install_requires=reqs
      )