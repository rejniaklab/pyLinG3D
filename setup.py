from setuptools import find_packages
from setuptools import setup

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()
    
with open("requirements.txt", "r", encoding="utf-8") as fh:
    required = fh.read()

setup(
    name="pyLinG3D",
    version="0.3",
    author="A. Hu, A.M.E. Ojwang, K.D. Olumoyin, and K.A. Rejniak",
    author_email="kasia.rejniak@moffitt.org, kayode.olumoyin@moffitt.org",
    description="Visualizing the Spatio-Temporal Dynamics of Clonal Evolution",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=find_packages(),
    install_requires=required,  # Add any dependencies here
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
    ],
    license="MIT",
    python_requires=">=3.8",
)
