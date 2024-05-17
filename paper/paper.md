---
title: 'SONIC: Software for Optical Navigation and Instrument Calibration'
tags:
  - optical navigation
  - camera calibration
  - projective geometry
  - MATLAB
authors:
  - name: John A. Christian
    affiliation: 1 # (Multiple affiliations must be quoted)
  - name: Michael Krause
    affiliation: 1
  - name: Ava C. Thrasher
    affiliation: 1
  - name: Michela Mancini
    affiliation: 1
  - name: Sébastien Henry
    affiliation: 1
  - name: Priyal Soni
    affiliation: 1
affiliations:
 - name: Georgia Institute of Technology, USA
   index: 1
date: TBD 2024
bibliography: paper.bib
---

# Summary

Optical navigation (OPNAV) is an increasingly common solution for spaceflight 
navigation in which information about the state of the spacecraft can be 
extracted from images. There are a wide variety of OPNAV techniques which can 
estimate position, velocity, and orientation information of the viewing camera. 
Many of these techniques rely on a host of similar image processing and 
geometry principles. Some of the necessities for successful OPNAV include 
fundamental camera models, projective geometry, and image processing abilities. 
`SONIC` is a MATLAB package which neatly integrates these tools for easy access. 

# Statement of need

`SONIC` is an object-oriented MATLAB toolkit for OPNAV, designed to simplify the 
development of commonly used and novel algorithms alike. It was developed by the 
researchers and students at the Space Exploration 
and Analysis Laboratory and is informed by extensive research in the field. 
`SONIC` provides an organized and transparent set of utilities to make this research 
readily available to projects looking to include OPNAV in their mission design 
and analysis.  There are many features included in `SONIC` to facilitate this. A 
strong basis of geometry and mathematics is the primary foundation for `SONIC` [@Hartley2003]. 
Geometric objects such as points, lines, planes, conics, and quadrics are 
represented in convenient ways that enable quick analysis of their meets, joins, 
and projections. Image processing routines directly useful for OPNAV such as image scans, 
edge detection, estimating background noise, and centroiding are included, and are
easily configurable depending on the desired use. `SONIC` also intuitively 
handles representation of camera models and can account for image distortion and 
aberration. In addition, a direct interface with the Hipparcos star catalog is 
included, allowing the user to easily access and filter its entries. 

These individual functionalities in `SONIC` can be used to synthesize larger 
capabilities. `SONIC` was deliberately created as a toolbox, and is not inclusive of 
complete workflows. This design allows for a multitude of useful tools that do 
not make assumptions about the users intentions.

There are existing software packages which separately address many of the required pieces 
for OPNAV or a particular type of OPNAV. GABLE, the 
Clifford Multivector Toolbox, and NASA’s Rigid Geometric Algebra [@RGA2023] all provide 
solutions for projective geometry (essential mathematics for OPNAV). There are 
also many well-known image processing toolkits (i.e. OpenCV) and those specific 
to certain OPNAV problems, such as Tetra3 or Astrometry.net [@tetra3][@astronet]. `SONIC` aims to 
capture the intersection of these mathematical principles, image processing, and 
optics needed for the growing, multi-disciplinary field of OPNAV.  Of particular note is
the Goddard Image Analysis and Navigation Toolkit (GIANT), which is a well-known optical 
navigation API [@GIANT]. The basis of GIANT  is scene-oriented, allowing the user to manipulate
objects' positions and orientations within a scene. The `SONIC` design seeks to complement
this and contribute a different approach to the OPNAV problem from a camera-oriented and 
projective geometry point of view, while adding a few different solutions to the problems of
star identification, triangulation [@Henry2023], and more.

# References
