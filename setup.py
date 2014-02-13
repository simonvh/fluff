from distutils.core import setup
import setuptools
from fluff.config import FL_VERSION

DESCRIPTION = """
fluff - plots and graphs 
"""

setup (name = 'fluff',
        version = FL_VERSION,
        description = DESCRIPTION,
        author='Simon van Heeringen',
        author_email='s.vanheeringen@ncmls.ru.nl',
        license='MIT',
        packages=[
            'fluff'
        ],
        scripts=[
            "scripts/fluff_bandplot.py",
            "scripts/fluff_profile.py",
            "scripts/fluff_heatmap.py",
        ],
        data_files=[]
)
