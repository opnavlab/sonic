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

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
aas-journal: Astrophysical Journal <- The name of the AAS journal.
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

`SONIC` is a MATLAB toolkit for OPNAV, designed with an object-oriented structure 
to simplify the development of commonly used and novel algorithms alike. This 
toolkit was developed by the researchers and students at the Space Exploration 
and Analysis Laboratory and is informed by extensive research in the field. 
`SONIC` provides an organized and transparent toolkit to make this research 
readily available to projects looking to include OPNAV in their mission design 
and analysis.  There are many features included in `SONIC` to facilitate this. A 
strong basis for geometry and mathematics is the primary foundation for `SONIC`. 
Geometric objects such as points, lines, planes, conics, and quadrics are 
represented in convenient ways that allow quick analysis for their meets, joins, 
and projections. Image processing directly useful for OPNAV such as image scans, 
edge detection, estimating background noise, and centroiding are included, and 
easily configurable depending on the desired use. `SONIC` also intuitively 
handles representation of camera models and can account for image distortion and 
aberration. In addition, a direct interface with the Hipparcos star catalog is 
included, allowing the user to easily access and filter its entries. These 
individual functionalities in `SONIC` can be used to synthesize larger 
capabilities. 

There are existing toolkits which separately address many of the required pieces 
for OPNAV or provide solutions for a particular type of OPNAV. GABLE, the 
Clifford Multivector Toolbox, and NASA’s Rigid Geometric Algebra all provide 
solutions for projective geometry (essential mathematics for OPNAV). There are 
also many well-known image processing toolkits (i.e. OpenCV) and those specific 
to certain OPNAV problems, such as Tetra3 or Astrometry.net. `SONIC` aims to 
capture the intersection of the mathematical principles, image processing, and 
optics needed for the growing, multi-disciplinary field of OPNAV. 


# Mathematics
Include paragraph or two of mathematical philosphy. Projective geometry, geometric algebra, meets and join

# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% }

# Acknowledgements

# References
