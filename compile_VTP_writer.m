% Script to compile the .vtp C++ writer 
%
% You will need a working MEX environment in MATLAB
%
% It needs to link against the VTK libraries provide by VMTK. Picking up 
% the header files from the version of VTK installed by homebrew
%
% The -D option is needed to get around an incompatibility between MATLAB
% R2013b and XCode 5.1.1
%
mex -Dchar16_t=uint16_T -I/usr/local/include/vtk-5.10/ -L/Library/Python/2.7/site-packages/vmtk-1.2-py2.7-macosx-10.9-intel.egg/vmtk/lib/ -lvtkGraphics -lvtkFiltering -lvtkIO -lvtkCommon VTP_writer.cpp
