import setuptools
from pathlib import Path
hih2_home = Path(__file__).parent
pypi_descrip = (hih2_home / "README.md").read_text()

setuptools.setup(
	name = "hih2",
	version = "1.0",
	author = "Ava Polzin",
	author_email = "apolzin@uchicago.edu",
	description = "Atomic-to-molecular hydrogen transition models for astrophysical simulations.",
	packages = ["hih2", "hih2/proj", "hih2/vol"],
	url = "https://github.com/avapolzin/hih2",
	license = 'MIT',
	classifiers = [
		"Development Status :: 5 - Production/Stable",
		"Intended Audience :: Science/Research",
		"License :: OSI Approved :: MIT License",
		"Operating System :: OS Independent",
		"Programming Language :: Python",
		"Topic :: Scientific/Engineering :: Astronomy",
		"Topic :: Scientific/Engineering :: Physics"],
	python_requires = ">=3",
	install_requires = ["astropy", "numpy"],
	long_description=pypi_descrip,
    long_description_content_type='text/markdown'
)