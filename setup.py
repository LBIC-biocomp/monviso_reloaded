"""Python setup.py for monviso_reloaded package"""
import io
import os
from setuptools import find_packages, setup


def read(*paths, **kwargs):
    """Read the contents of a text file safely.
    >>> read("monviso_reloaded", "VERSION")
    '0.1.0'
    >>> read("README.md")
    ...
    """

    content = ""
    with io.open(
        os.path.join(os.path.dirname(__file__), *paths),
        encoding=kwargs.get("encoding", "utf8"),
    ) as open_file:
        content = open_file.read().strip()
    return content


def read_requirements(path):
    return [
        line.strip()
        for line in read(path).split("\n")
        if not line.startswith(('"', "#", "-", "git+"))
    ]


setup(
    name="monviso",
    python_requires='>=3.9',
    version=read("monviso_reloaded", "VERSION"),
    description="MoNvIso is a comprehensive software tool designed for the analysis and modeling of protein isoforms. It automates the process of identifying canonical and additional isoforms, assessing their modeling propensity, mapping mutations accurately, and building structural models of proteins.",
    url="https://github.com/alisamalb/monviso_reloaded",
    long_description=read("README.md"),
    long_description_content_type="text/markdown",
    author="S. Albani",
    packages=find_packages(exclude=["tests", ".github"]),
    install_requires=read_requirements("requirements.txt"),
    entry_points={
        "console_scripts": ["monviso = monviso_reloaded.__main__:main"]
    },
    extras_require={"test": read_requirements("requirements-test.txt")},
)
