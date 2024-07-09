
Installing SONIC
=================================

Minimal preparation is needed to begin using SONIC. Simply clone the repository locally, and add the sonic directory to the search path that you're working in (i.e. using addpath).

A few MATLAB live tutorials are provided under +examples for demonstrative purposes:

    ScanLinesTutorial: Demonstrates the generic ability to extract image values along a specified scanline, which is useful for horizon-based OPNAV.
    SyntheticStarImgTutorial: Demonstrates creating a synthetic image of the star Arcturus, and its surrounding starfield in SONIC. Being able to construct an expected image is an essential operation for OPNAV.
    TriangVestaReconTutorial: Illustrates the use of the triangulation in SONIC to perform 3D structure reconstruction on the asteroid Vesta from features extracted in images. A few different approaches to triangulation are provided in SONIC, which can be used both for estimating the spacecraft state or reconstruction.

To verify your installation of SONIC, please run these examples locally and compare the outputs to the pre-run demos in the SONIC API documentation, linked below.

Note: SONIC requires the following MATLAB toolboxes installed:

    Image Processing Toolbox version 24.1
    Computer Vision System Toolbox version 24.1
