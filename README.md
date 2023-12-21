# Determine Knot Crossings

## Table of Contents

1. [Overview](#overview)
2. [Requirements](#requirements)
3. [Run the Code](#run-the-code)
4. [Vect File Format](#vect-file-format)
5. [Structure](#structure)
6. [Algorithm](#algorithm-overview)
7. [Results](#results)

## Overview
This code reads in input containing the information on a knot with the x,y,z coords of all the vertices of the shape. This is then used to calculate where the knot crosses on itself and label the resulting knot

## Requirements
A C++ compiler (e.g., g++)
Basic knowledge of the "vect" file format

## Run the Code
To compile the code, use any c compiler, then to run the code run ensure to pipe in the input file in the "vect" format.

## Vect File Format
The first line should be "VECT," followed by the number of components and the total number of lines. The subsequent lines should contain the 3D coordinates of points and their connections to form line segments. Comments are only allowed at the end of the line starting with the # symbol.


## Structure
The code is structured into several sections:

Data Structures: Defines structures for points, edges, crossings, and a custom priority queue.
Sorting Functions: Implements merge sort for points and edges.
Intersection Testing: Provides functions to test for intersections between line segments.
Main Function: Reads input, processes line segments, detects intersections, and prints results.

## Algorithm Overview
Read input points and edges from the "vect" file.
Sort the points from left to right using merge sort.
Process line segments and detect intersections using a custom priority queue.
Identify and label intersections, storing the information for later analysis.
Output the results, including the type of intersection and relevant codes.

## Results
The program prints the crossing code information, labelling the knot and where the different strands of the crossing connect to the other knot strands

