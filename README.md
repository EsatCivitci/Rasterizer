# Rasterizer
Software Rasterizer

Software Rasterizer is a C++ application that implements a forward rendering pipeline for 3D scenes. Developed as part of the CENG 477 Introduction to Computer Graphics course at Middle East Technical University, this project converts 3D scenes described in an XML file into 2D images.

Features

Modeling Transformations: Applies translations, rotations, and scalings to meshes.
Viewing Transformations: Transforms 3D world coordinates to camera space using both orthographic and perspective projections.
Rasterization: Renders triangles in:
Wireframe Mode: Uses the midpoint algorithm for edges with clipping (using an algorithm like Cohen-Sutherland or Liang-Barsky).
Solid Mode: Uses barycentric coordinates for triangle filling and includes depth buffering.
Backface Culling: Optionally discards triangles facing away from the camera.
Color Interpolation: Smoothly interpolates vertex colors along edges and across surfaces.
Multiple Camera Support: Allows several camera configurations per scene.
File Structure

Source Files:
Scene.cpp and Helpers.cpp handle input parsing and basic mathematical operations.
Additional modules implement transformations, rasterization, clipping, and culling.
Makefile:
Contains build instructions to compile the project.
