from os import walk
from os.path import abspath, normpath, splitext
from os.path import join as pj

from setuptools import setup, find_packages

pkg_data = {}

nb = len(normpath(abspath("src/crop_energy_balance"))) + 1
data_rel_pth = lambda pth: normpath(abspath(pth))[nb:]

data_files = []
for root, dnames, fnames in walk("src/crop_energy_balance"):
    for name in fnames:
        if splitext(name)[-1] in ['.json', '.ini', '.csv']:
            data_files.append(data_rel_pth(pj(root, name)))

pkg_data['crop_energy_balance'] = data_files

setup(
    name='crop_energy_balance',
    version="1.0.0",
    description="crop energy balance",
    long_description="A model for simulating energy balance in the soil-crop-atmosphere continuum",
    author="Rami Albasha",
    author_email="rami.albacha@yahoo.com",
    url='https://github.com/RamiALBASHA/crop_energy_balance',
    license='private',
    zip_safe=False,

    packages=find_packages('src'),
    package_dir={'': 'src'},

    package_data=pkg_data,
    setup_requires=[
        "pytest-runner",
    ],
    install_requires=[
        "numpy",
    ],
    tests_require=[
        "coverage",
        "pytest",
        "pytest-cov",
        "pytest-mock",
    ],
    entry_points={},

    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: BSD License",
        "Natural Language :: English",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3 :: Only",
        "Programming Language :: Python :: 3.7",
    ],
)