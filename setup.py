from setuptools import setup

def readme():
    with open('README.md') as f:
        README = f.read()
    return README


setup(
    name="Quantum-Tomography",
    version="1.0.2",
    description="A python library to help perform tomography on a quantum state.",
    long_description=readme(),
    long_description_content_type="text/markdown",
    url="https://github.com/KwiatQIM/Quantum-Tomography",
    author="Scott Turro",
    author_email="turroscott@gmail.com",
    license="MIT",
    classifiers=[
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
    ],
    packages=["QuantumTomography"],
    include_package_data=True,
    install_requires=["scipy","numpy","matplotlib"],
    entry_points={
        "console_scripts": [
            "Quantum-Tomography=QuantumTomography.cli:main",
        ]
    },
)