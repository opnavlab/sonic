![SONIC](https://github.com/opnavlab/sonic/blob/main/logos/sonic_banner.png)
# **S**oftware for **O**ptical **N**avigation and **I**nstrument **C**alibration
This open-sourced, object oriented MATLAB package is a toolkit built on years of research in
optical navigation (OPNAV), with a foundation in principled projective geometry. Many 
capabilities needed for space-flight OPNAV are made transparent and configurable within the SONIC classes, 
with minimal assumption of the analyst's intentions. This allows each user to construct their own 
unique workflows for a specific research or analysis problems. SONIC was designed to help simplify navigation solutions such as (but not limited to):

- Horizon-Based OPNAV
- Triangulation
- Star Identification
- 3D-Reconstruction
- And more

SONIC was created by the researchers and students of the [Space Exploration Analysis 
Laboratory](https://seal.ae.gatech.edu/) at the Georgia Institute of Technology for student, 
academic researchers, space science professionals, or anyone curious to learn more about OPNAV.

For **API documentation**, please visit [https://opnavlab.github.io/sonic/](https://opnavlab.github.io/sonic/).

# Getting Started
Minimal preparation is needed to begin using SONIC. Simply clone the repository locally,
and add the sonic directory to the search path that you're working in (i.e. using addpath). 

A few MATLAB live tutorials are provided under [+examples](https://github.com/opnavlab/sonic/tree/main/%2Bexamples)
for demonstrative purposes:
- **ScanLinesTutorial**: Demonstrates the generic ability to extract image values along a specified scanline, which is useful for horizon-based OPNAV.
- **SyntheticStarImgTutorial**: Demonstrates creating a synthetic image of the star Arcturus, and its surrounding starfield in SONIC. Being able to construct an expected image is an essential operation for OPNAV.
- **TriangVestaReconTutorial**: Illustrates the use of the triangulation in SONIC to perform 3D structure reconstruction on the asteroid Vesta from features extracted in images. A few different approaches to triangulation are provided in SONIC, which can be used both for estimating the spacecraft state or reconstruction.

To verify your installation of SONIC, please run these examples locally and compare the outputs to the pre-run demos in the SONIC API documentation, linked below.

Note: SONIC requires the following MATLAB toolboxes installed:
- Image Processing Toolbox version 24.1
- Computer Vision System Toolbox version 24.1

# Contributing
If you are interesting in contributing, have a feature request, or have found a bug in SONIC, please checkout [CONTRIBUTING.md](https://github.com/opnavlab/sonic/blob/main/CONTRIBUTING.md) for more information.

# References
Please see [REFERENCES.md](https://github.com/opnavlab/sonic/blob/main/REFERENCES.md) for a list of references used in the development of SONIC.
