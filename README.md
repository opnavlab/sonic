![SONIC](https://github.com/opnavlab/sonic/blob/main/logos/sonic_banner.png)
# **S**oftware for **O**ptical **N**avigation and **I**nstrument **C**alibration
This open-sourced, object oriented MATLAB package is a toolkit built on years of research in
optical navigation (OPNAV), with a foundation in principled projective geometry. Many 
capabilities needed for OPNAV are made transparent and configurable within the SONIC classes, 
with minimal assumption of the analyst's intentions, such that full workflows can be 
created and applied to mission specific or individual problems. 

SONIC was created by the researchers and students of the [Space Exploration and Design 
Laboratory](https://seal.ae.gatech.edu/) at the Georgia Institute of Technology for academics, 
industry projects, or anyone curious to learn more about OPNAV.

# Getting Started
Minimal preparation is needed to begin using SONIC. Simply clone the repository locally,
and point your working directory to the SONIC path. 

A few MATLAB live tutorials are provided under [+examples](https://github.com/opnavlab/sonic/tree/main/%2Bexamples)
for demonstrative purposes:
- **ScanLinesTutorial**: Demonstrates the generic ability to extract image values along a specified scanline, which is useful for horizon-based OPNAV.
- **SyntheticStarImgTutorial**: Demonstrates creating a synthetic image of the star Arcturus, and its surrounding starfield in SONIC. Being able to construct an expected image is an essential operation for OPNAV.
- **TriangVestaReconTutorial**: Illustrates the use of the triangulation in SONIC to perform 3D structure reconstruction on the asteroid Vesta from features extracted in images. A few different approaches to triangulation are provided in SONIC, which can be used both for estimating the spacecraft state or reconstruction.

Note: As of 5/16/2024, SONIC has a few MATLAB toolbox dependencies which will
require installation prior to full operation including:
- Image Processing Toolbox
- Computer Vision System Toolbox

# Feedback
If you have a suggestion, please use the GitHub Issues tab to let us know!

# References
Please see REFERENCES.md for a list of references used in the development of SONIC.
