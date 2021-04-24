
from setuptools import setup

with open("README.md", 'r') as f:
    long_description = f.read()

with open('requirements.txt') as f:
    required = f.read().splitlines()

setup(
    name='ksrates',
    version='0.3',
    packages=['ksrates', 'wgd_ksrates'],
    url='https://github.com/VIB-PSB/ksrates',
    license='GNU GPL v3.0',
    author='Cecilia Sensalari, Steven Maere, Rolf Lohaus',
    author_email='steven.maere@ugent.vib.be, rolf.lohaus@ugent.vib.be',
    description='Command line tool for substitution rate-adjustment in mixed paralog-ortholog Ks distributions.',
    long_description=long_description,
    package_data={
        'ksrates': ['ks.mplstyle']
    },
    py_modules=['ksrates_cli'],
    install_requires=[required],
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent"
    ],
    python_requires='>=3.6',
    entry_points={
        "console_scripts": [
            "ksrates = ksrates_cli:cli"
        ]
    },
)
