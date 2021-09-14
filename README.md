<p align="center">
  <h3 align="center">🌊 Fluid Simulation 🌊</h3>
  <p align="center">
    C++ fluid simulation solver implemented from scratch !
    <br />
    <a href="#screenshots"><strong>Screenshots »</strong></a>
    <br />
  </p>
</p>



<!-- TABLE OF CONTENTS -->
<details open="open">
  <summary>Table of Contents</summary>
  <ol>
    <li>
      <a href="#about-the-project">About The Project</a>
      <ul>
        <li><a href="#built-with">Built With</a></li>
      </ul>
    </li>
    <li>
      <a href="#getting-started">Getting Started</a>
      <ul>
        <li><a href="#prerequisites">Prerequisites</a></li>
        <li><a href="#installation">Installation</a></li>
      </ul>
    </li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#screenshots">Screenshots</a></li>
    <li><a href="#contact">Contact</a></li>
    <li><a href="#references">References</a></li>
    <li><a href="#license">License</a></li>
  </ol>
</details>



## About The Project


During my 6 months NII International Internship Program 2020 supervised by [Ryoichi Ando](https://ryichando.graphics/), I studied and implemented from scratch in C++ some research papers for fluid simulation. Here is a recap of what I've been working on :

#### Fluid Simulation
- **Advection**
	- Semi-Lagrangian scheme
	- MacCormack
- **Projection**
	- Pressure solver with Neumann boundary condition
- **Linear solver**
	- Gauss-Seidel relaxation
	- Conjugate Gradient Method
		- Modified Incomplete Cholesky Level Zero preconditionning
- **Grid structure**
	- Basic
	- Staggered
- **Level-Set**
	- With reinitialisation

#### Rendering 
- **OpenGL** : Volume Ray Marching
- **Mitsuba** : Volume File export


#### Built With
* [Eigen](https://gitlab.com/libeigen/eigen) C++ template library for linear algebra: matrices, vectors, numerical solvers, and related algorithms.
* [GLM](https://github.com/g-truc/glm) Header only C++ mathematics library for graphics software based on the OpenGL Shading Language (GLSL) specifications.
* [GLFW](https://github.com/glfw/glfw) Open Source, multi-platform library for OpenGL. It provides a simple API for creating windows, contexts and surfaces, receiving input and events.
* [glad](https://github.com/Dav1dde/glad) GL/GLES/EGL/GLX/WGL Loader-Generator based on the official specs.
* [inipp](https://github.com/mcmtroffaes/inipp) Simple header-only C++ .ini parser and generator.
* [stb](https://github.com/nothings/stb) Image writing to disk: PNG, TGA, BMP.
* [tinyply](https://github.com/ddiakopoulos/tinyply) A single-header, zero-dependency (except the C++ STL) public domain implementation of the PLY mesh file format.

## Getting Started
### Prerequisites

### Installation
1. Clone the repo with the submodules
   ```sh
   git clone --recurse-submodules git@github.com:Trietch/fluid-simulation.git
   ```
2. Compile
   ```sh
   mkdir build
   cd build
   cmake -DCMAKE_BUILD_TYPE=Release ..
   make -j
   ```

## Usage

1. Edit the configuration file `config.ini`
2. Execute the `fluid-simulation` binary file you just compiled
   ```sh
   ./fluid-simulation
   ```

## Screenshots

## Contact
 » [contact@tmarrec.dev](mailto:contact@tmarrec.dev)

## References
1.  **Stable fluids**  
    Jos Stam. 1999. In Proceedings of the 26th annual conference on Computer graphics and interactive techniques - SIGGRAPH'99. ACM Press.  [doi.org/10.1145/311535.311548](https://doi.org/10.1145/311535.311548)
2.  **Fluid simulation for computer graphics**  
    Robert Bridson. 2016. CRC Press, Taylor & Francis Group, CRC Pressis an imprint of the Taylor & Francis Group, an informa Business, Boca Raton.  [doi.org/10.1201/9781315266008](https://doi.org/10.1201/9781315266008)
3.  **Fluid simulation**  
    Robert Bridson and Matthias Müller-Fischer. 2007. In ACM SIGGRAPH 2007 courses on SIGGRAPH'07. ACM Press.  [doi.org/10.1145/1281500.1281681](https://doi.org/10.1145/1281500.1281681)
4.  **Practical animation of liquids**  
    Nick Foster and Ronald Fedkiw. 2001. In Proceedings of the 28th annual conference on Computer graphics and interactive techniques - SIGGRAPH'01. ACM Press.  [doi.org/10.1145/383259.383261](https://doi.org/10.1145/383259.383261)
5.  **Adventures in Fluid Simulation**  
    Miles Macklin. 2010.  [blog.mmacklin.com/2010/11/01/adventures-in-fluid-simulation](https://blog.mmacklin.com/2010/11/01/adventures-in-fluid-simulation/)
6.  **An Unconditionally Stable MacCormack Method**  
    Andrew Selle, Ronald Fedkiw, ByungMoon Kim, Yingjie Liu, and Jarek Rossignac. 2007. Journal of Scientific Computing 35, 2-3 (Nov. 2007), 350–37.  [doi.org/10.1007/s10915-007-9166-4](https://doi.org/10.1007/s10915-007-9166-4)
7.  **Real-Time Fluid Dynamics for Games**  
    Jos Stam. 2003.  [http://graphics.cs.cmu.edu/nsp/course/15-464/Spring11/papers/StamFluidforGames.pdf](http://graphics.cs.cmu.edu/nsp/course/15-464/Spring11/papers/StamFluidforGames.pdf)
8.  **Real-Time Simulation and Rendering of 3D Fluids**  
    Keenan Crane, Ignacio Llamas, and Sarah Tariq. 2007. in GPU Gems 3. Chapter 30.  [dl.acm.org/doi/10.5555/1407436](https://dl.acm.org/doi/10.5555/1407436)
9.  **Visual simulation of smoke**  
    Ronald Fedkiw, Jos Stam, and Henrik Wann Jensen. 2001. In Proceedings of the 28th annual conference on Computer graphics and interactive techniques - SIGGRAPH'01. ACM Press.  [doi.org/10.1145/383259.383260](https://doi.org/10.1145/383259.383260)
10.  **Realistic Animation of Liquids**  
    Nick Foster and Dimitri Metaxas. 1996. Graphical Models and Image Processing 58, 5 (Sept. 1996), 471–483.  [doi.org/10.1006/gmip.1996.0039](https://doi.org/10.1006/gmip.1996.0039)
11.  **Mitsuba 2**  
    [mitsuba-renderer.org](https://www.mitsuba-renderer.org/)
12.  **OpenMP**  
    [openmp.org](https://www.openmp.org/)

## License
Distributed under the MIT License. See `LICENSE` for more information.
