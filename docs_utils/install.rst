
Installing SONIC
=================================

SONIC requires little set up to begin using. Please ensure you have MATLAB installed on your computer. As of 5/16/2024, please install the following toolkit dependencies:
- Image Processing Toolbox
- Computer Vision System Toolbox

Once you are ready to begin using SONIC, clone the `SONIC repository <https://github.com/opnavlab/sonic>`_ to the local directory you wish to access it from. After cloning, 
you will see the sonic folder, which contains the package directory +sonic and the accompanying +examples.

The +sonic directory contains the collection of classes that make up SONIC, while te +examples directory contains a few MATLAB 
live demonstrations of SONIC use-cases.

To begin using SONIC, all that is left to do is open a new script, and add SONIC to the working path (i.e. using addpath). SONIC is an object oriented
software package, so from here, you can follow the typical object-oriented programming (OOP) style. For example, you can instantiate objects using:
obj = sonic.ClassNameHere(inputs)
access their properties via:
prop = obj.propertyNameHere
and call class methods as:
output = obj.methodNameHere(inputs)
or static methods as
output = sonic.ClassNameHere.methodNameHere(inputs)

For more information on MATLAB OOP format, please visit their documentation `here <https://www.mathworks.com/products/matlab/object-oriented-programming.html>`_.