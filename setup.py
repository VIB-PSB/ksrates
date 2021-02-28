
from setuptools import setup

with open("README.md", 'r') as f:
    long_description = f.read()

with open('requirements.txt') as f:
    required = f.read().splitlines()

setup(
    name='ksrates',
    version='0.1',
    packages=['ksrates', 'wgd'],
    url='xxx',
    license='xxx',
    author='xxx',
    author_email='xxx@psb.vib-ugent.be',
    description='Command line tool for substitution rate-correction of mixed Ks distributions.',
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