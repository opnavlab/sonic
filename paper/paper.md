---
title: 'SONIC: Software for Optical Navigation and Instrument Calibration'
tags:
  - optical navigation
  - camera calibration
  - projective geometry
  - MATLAB
authors:
  - name: Ava C. Thrasher
    affiliation: 1
  - name: Michael Krause
    affiliation: 1
  - name: Sébastien Henry
    affiliation: 1   
  - name: Michela Mancini
    affiliation: 1
  - name: Priyal Soni
    affiliation: 1
  - name: John A. Christian
    affiliation: 1 
affiliations:
 - name: Georgia Institute of Technology, USA
   index: 1
date: TBD 2024
bibliography: paper.bib
---

# Summary

Optical navigation (OPNAV) represents an increasingly common solution for spaceflight navigation in which information on the state of the spacecraft is extracted from images. A wide variety of OPNAV techniques allow the estimation of the position, velocity, orientation, and calibration of the viewing camera, many of which rely on a host of similar image processing and geometry principles. Essential components for successful OPNAV include projective geometry, camera models, and image processing routines. `SONIC` is a MATLAB package that integrates these tools into a common framework for easy access.

# Statement of need

`SONIC` is an object-oriented MATLAB toolkit for OPNAV, designed to simplify and accelerate the production of new research results in this field. The `SONIC` toolkit provides an organized and transparent set of utilities to make this research readily available to projects looking to include OPNAV in their mission design and analysis.

Projective geometry [@semple1952] [@richtergebert2011] and geometric algebra [@lengyel2024] is the foundation for `SONIC`. Geometric objects such as points, lines, planes, conics, and quadrics are conveniently represented in ways that enable quick analysis of their meets, joins, and projections. Image processing routines such as image scans, edge detection, estimating background noise, and centroiding are included and easily configurable depending on the desired use. `SONIC` also intuitively handles the representation of camera models and can account for optical distortions and stellar aberration. In addition, `SONIC` provides a direct interface to the Hipparcos star catalog, allowing the user to easily access and filter its entries. Finally, `SONIC` offers solutions to the commonly encountered problems of star identification [@mortari2004][@lang2010], triangulation [@henry2023], Horizon-based OPNAV [@christian2021], and more. These individual functionalities in `SONIC` can be used to synthesize larger capabilities. `SONIC` was deliberately created as a toolbox and does not include complete workflows for end-to-end capabilities. This design choice allows `SONIC` to focus on the foundational building blocks of OPNAV while making minimal assumptions about the user's application. Therefore, the user is given the freedom and flexibility to assemble the `SONIC` building blocks in their own way to produce a unique OPNAV capability or address a particular research problem. 

As demonstrations for applying SONIC to OPNAV problem, multiple example workflows have been provided in the SONIC repository. These cover only a few examples, which are shown in Figs. 1 through 3. For more information on how these results were obtained using SONIC, the reader is encouraged to run these examples locally.

![](https://raw.githubusercontent.com/opnavlab/sonic/main/paper/scan_lines.png)
**Figure 1.** SONIC can be used to quickly perform a scan and do analysis on an image.

![](https://raw.githubusercontent.com/opnavlab/sonic/main/paper/smear_removed.png)
**Figure 2.** SONIC image processing techniques were used to remove image smear from this image of the asteroid Bennu.

![](https://raw.githubusercontent.com/opnavlab/sonic/main/paper/synth_img.png)
**Figure 3.** Synthetic image of the star Arcturus, and its surrounding starfield sized by magnitude. Access to star catalogs in SONIC and projective geometry tools allows for easy generation of synthetic images.

There are a number of existing software packages that address similar problems as `SONIC`. Of particular note is the Goddard Image Analysis and Navigation Toolkit (`GIANT`), which is a well-known OPNAV API [@giant]. `GIANT` is scene-oriented, allowing the user to manipulate objects' positions and orientations within a scene, and then perform a multitude of OPNAV routines. `SONIC` seeks to complement this by design by approaching the OPNAV problem from a camera-oriented and projective geometry point of view, while also adding a few different solutions for state estimation. Other software packages have solutions for different aspects of OPNAV. For example, `GABLE` [@gable], the `Clifford Multivector Toolbox` [@sangwine2016], and NASA’s `Rigid Geometric Algebra` [@rga] all provide generic computing frameworks for geometric algebra. While `SONIC` does use concepts from geometric algebra, generic geometric algebera manipulations are not required. Furthermore, there are well-known image processing toolkits for general applications such as image processing (i.e. `OpenCV` [@opencv]) or astronomical image registration (i.e. `Tetra3` [@tetra3] and `Astrometry.net` [@lang2010][@astronet]). `SONIC` captures the intersection of principles in projective geometry, image processing, and optics needed for the growing multi-disciplinary field of OPNAV. 

# References

