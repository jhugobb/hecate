# Hecate
Hecate is a program that creates volume data of triangle meshes and writes that data in a lossless compressed format efficiently.

## What was the purpose of this?
As part of my Master's thesis "Compressed Representation of 3D Models for 3D Printing" I created this C++ implementation from the ground-up

## What's the code quality?
It's not the cleanest code ever written, but it's mostly documented and should be simple to follow. The core is written in Grid.cpp, where there are the coloring functions. The compressed encoding formats are in HecateFiles.cpp.

## Are there any configurable parameters? 
Yes! Almost anything that requires parameters is configurable, as well as deciding what if any format to write in. You can also decide if you want full coloring or not.

## What are these python scripts?
I created some scripts to streamline a bash version of the program and make the result generation process faster. There are also scripts for compression and decompression using 7zip and Gzip.

## How do I build the repository?
Firstly, this is a Linux only implementation, and was developed in a Linux Mint OS. Also you are going to need some dependencies:
 - Boost
 - OpenGL
 - GLUT
 - GLM
 - OpenMP
 
 Then you need CMake to build it. Create a build folder and inside it run:
 
 `cmake .. && make`
 
 And the binary `hecate` should be built. The only argument is the PLY of the mesh you want to use. Make sure that the PLY is not in binary format.
 
 ## Anything I should keep an eye on?
 Yes, the program only processes one slice of the model at a time, in order to not have a huge memory consumption. However, it's heavily parallelized with openMP, so having resolutions higher than 4096 may not be reccomended. Also PNG format takes a lot of time to write, keep that in mind.
 
 
