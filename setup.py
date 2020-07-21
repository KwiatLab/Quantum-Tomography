from setuptools import setup

"""
Copyright 2020 University of Illinois Board of Trustees.
Licensed under the terms of an MIT license
"""

__author__ = 'Quoleon/Turro'
"""CHECK OUT THE REFERENCE PAGE ON OUR WEBSITE :
http://research.physics.illinois.edu/QI/Photonics/Quantum-Tomography_lib_Ref/"""

def readme():
    with open('README.md') as f:
        README  = f.read()
    return README


setup(
    name = "Quantum-Tomography",
    version = "1.0.3.0",
    description = "A python library to help perform tomography on a quantum state.",
    long_description = readme(),
    long_description_content_type = "text/markdown",
    url = "https://github.com/KwiatQIM/Quantum-Tomography",
    author = "Scott Turro",
    author_email = "turroscott@gmail.com",
    license = "MIT",
    classifiers = [
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
    ],
    packages = ["src/QuantumTomography"],
    include_package_data = True,
    install_requires = ["scipy", "numpy", "matplotlib"],
    entry_points = {
        "console_scripts": [
            "Quantum-Tomography = QuantumTomography.cli:main",
        ]
    },
)
