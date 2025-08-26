import setuptools
import os

with open(f"{os.path.dirname(__file__)}/README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="cplt3d",
    version="0.1",
    author="Benjamin Cohen",
    author_email="benmckcohen@uchicago.edu",
    description="A package for nifty 3d plotting in matplotlib",
    packages=["cplt3d"],
    install_requires = ['numpy','scipy','matplotlib','imageio'],
    long_description=long_description,
    long_description_content_type="text/markdown",
    license='MIT'
)