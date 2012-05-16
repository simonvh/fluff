from distutils.core import setup

VERSION = "1.0"
DESCRIPTION = """
fluff - plots and graphs 
"""

setup (name = 'fluff',
		version = VERSION,
		description = DESCRIPTION,
		author='Simon van Heeringen',
		author_email='s.vanheeringen@ncmls.ru.nl',
		license='MIT',
		packages=[
			'fluff'
		],
		scripts=[
			"scripts/profile_screenshot.py"
			"scripts/heatmap.py"
			"scripts/cluster_graph.py"
		],
		data_files=[]
)
