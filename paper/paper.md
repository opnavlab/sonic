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

Optical navigation (OPNAV) represents an increasingly common solution for spaceflight navigation in which information on the state of the spacecraft is extracted from images. A wide variety of OPNAV techniques allow the estimation of the position, velocity, orientation, and calibration of the viewing camera, many of which rely on a host of similar image processing and geometry principles. Essential components for a successful OPNAV include fundamental camera models, projective geometry, and image processing abilities. `SONIC` is a MATLAB package that cohesively integrates these tools for easy access.

# Statement of need

`SONIC` is an object-oriented MATLAB toolkit for OPNAV, designed to simplify the development of both commonly used and novel algorithms. It was developed by researchers and students at the Space Exploration and Analysis Laboratory and thus benefits from extensive research experience in the field. `SONIC` provides an organized and transparent set of utilities to make this research readily available to projects looking to include OPNAV in their mission design and analysis.

A strong basis of geometry and mathematics is the primary foundation for `SONIC` [@Hartley2003]. Geometric objects such as points, lines, planes, conics, and quadrics are conveniently represented in ways that enable quick analysis of their meets, joins, and projections. Image processing routines such as image scans, edge detection, estimating background noise, and centroiding are included and easily configurable depending on the desired use. `SONIC` also intuitively handles the representation of camera models and can account for image distortion and aberration. In addition, `SONIC` provides a direct interface to the Hipparcos star catalog, allowing the user to easily access and filter its entries. Finally, `SONIC` offers solutions to the commonly encountered problems of star identification, triangulation [@Henry2023], Horizon-based OPNAV [@Christian2021], and more. These individual functionalities in `SONIC` can be used to synthesize larger capabilities. `SONIC` was deliberately created as a toolbox and does not include complete workflows. This design allows for a multitude of useful tools that do not make assumptions about the user's intentions.

Existing software packages separately address many of the required pieces for OPNAV or particular aspects of OPNAV. Of particular note is the Goddard Image Analysis and Navigation Toolkit (GIANT), which is a well-known optical navigation API [@GIANT]. GIANT is scene-oriented, allowing the user to manipulate objects' positions and orientations within a scene, and then perform a multitude of OPNAV routines. `SONIC` seeks to complement this by design by approaching the OPNAV problem from a camera-oriented and projective geometry point of view, while also adding a few different solutions for state estimation. Other software packages have solutions for different aspects of OPNAV. For example, GABLE [@GABLE], the Clifford Multivector Toolbox [@Sangwine2016], and NASA’s Rigid Geometric Algebra [@RGA] among others all provide solutions for projective geometry (essential mathematics for OPNAV). Furthermore, there are well-known image processing toolkits for general applications such as image processing (i.e. OpenCV [@OpenCV]) or astronomical image recognition (i.e. Tetra3 [@tetra3] and Astrometry.net [@astronet]). `SONIC` captures the intersection of principles in projective geometry, image processing, and optics needed for the growing multi-disciplinary field of OPNAV. 

# References

