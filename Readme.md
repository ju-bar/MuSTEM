#  μSTEM <img src="https://github.com/ju-bar/MuSTEM/blob/master/Manual/Figures/350x350_inelastic_cbed.png" width="64" height="64" />

μSTEM is a transmission electron microscopy (TEM) simulation suite, in particular for scanning transmission electron microscopy (STEM) images, that was developed mainly at the University of Melbourne. The computing suite is based on the multislice
method. More detail can be found in the [manual](https://github.com/ju-bar/MuSTEM/blob/master/Manual/muSTEM_manual.pdf) and the following scientific paper:

[Modelling the inelastic scattering of fast electrons,
L.J. Allen, A.J. D'Alfonso and S.D. Findlay,
Ultramicroscopy, Vol. 151, pp. 11-22, (2015).](http://www.sciencedirect.com/science/article/pii/S0304399114002034)

<img src="https://github.com/ju-bar/MuSTEM/blob/master/Manual/Figures/512x512_out_PACBED.png" width="512" height="512" />

This is a fork of the [original µSTEM repository](https://github.com/HamishGBrown/MuSTEM).

## Prerequisites

### GPU version

* 64 bit windows OS
* A CUDA enabled GPU with compute class 3.0 or greater
* Version 10.1 of the [Nvidia CUDA toolkit](https://developer.nvidia.com/cuda-toolkit-archive), this installs the .dll files necessary to run the GPU versions of μSTEM.

#### To compile source code

* [PGI Visual Fortran and Microsoft Visual Studio](https://www.pgroup.com/products/pvf.htm)
* [Intel oneAPI Fortran Compiler](https://www.intel.com/content/www/us/en/developer/tools/oneapi/fortran-compiler.html)
* [Intel Math Kernel Library](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl.html)
* [Microsoft Visual Studio (Community)](https://visualstudio.microsoft.com/de/vs/community/)

### CPU version

* 64 bit Windows OS

#### To compile source code

* Any Fortran 90 compiler, FFTW libraries, or the Intel MKL

## Precompiled executables

Precompiled executables can be found in the [Executables](https://github.com/ju-bar/MuSTEM/tree/master/Executables) folder of the repository or by clicking the links below:

### Current versions

**Version 6.0**
* [**CPU** single precision](https://github.com/ju-bar/MuSTEM/blob/master/Executables/CPU_muSTEM_x64_v6.0_single_precision.zip)
* [**GPU** single precision](https://github.com/ju-bar/MuSTEM/blob/master/Executables/CUDA_muSTEM_x64_v6.0_single_precision.zip)

**Remark:** The GPU version is distributed without the CUDA dynamic link libraries. The executable links to cudart64_101.dll and to cufft64_10.dll, which are distributed with the [NVIDIA CUDA Toolit Version 10.1](https://developer.nvidia.com/cuda-10.1-download-archive-base).
The CPU version is linked against the Intel OpenMP libraries. It requires libiomp5md.dll, which is provided in the zip archive with the executable.

## Tutorial

Click [here](https://minhaskamal.github.io/DownGit/#/home?url=https://github.com/ju-bar/MuSTEM/tree/master/Tutorial) to download the μSTEM tutorial calculations, powerpoint and activity sheet.

### Compiling source code


The GPU version of MuSTEM is built using the [PGI Fortran compiler and Microsoft Visual Studio 2015](https://www.pgroup.com/products/pvf.htm), please make sure that this software is correctly installed before proceeding. Create a new Visual Studio project and add the source code contained in this repository. The GPU version of the code requires the source files in the GPU_routines folder and the CPU only version of the code requires the source files in the CPU_routines folder. Modify the project properties so that Microsoft Visual Studio passes the following commands to the PGI compiler:

Build commands:

-Mpreprocess -Bstatic -Mbackslash -mp -Mcuda=cuda9.0 -I"C:\Program Files\PGI\win64\17.3\include" -I"c:\program files\pgi\win64\17.3\include" -I"C:\Program Files\PGI\Microsoft Open Tools 14\include" -I"C:\Program Files (x86)\Windows Kits\10\Include\shared" -I"C:\Program Files (x86)\Windows Kits\10\Include\um" -fast -ta=tesla -Minform=warn 

To build the single precision version add the command -Dsingle_precision. For the Double precision version add the command -Ddouble_precision. To build the GPU version (requires the PGI compiler) add the command -Dgpu. The code also requires recursive routines to be enabled (The  /recursive command in the Intel Visual Fortran compiler) for correct calculation of the absorptive form factors.

Linker commands:

-Bstatic -mp -Mcuda=cuda8.0 -ta=tesla -o"$(Outdir)\MuSTEM_Open.exe" -Wl,/libpath:"C:\Program Files\PGI\win64\17.3\lib" -Wl,/libpath:"C:\Program Files\PGI\win64\2017\cuda\9.0\lib64" cufft.lib 

The links to the PGI CUDA libraries and Windows kits may need to be modified depending on the install directories. $(Outdir) represents the output directory of the Microsoft Visual Studio build.

Some example simulations are included on the [MuSTEM website](http://tcmp.ph.unimelb.edu.au/mustem/download.php).

The CPU version of MuSTEM can also be built using the
* [Intel oneAPI Fortran Compiler](https://www.intel.com/content/www/us/en/developer/tools/oneapi/fortran-compiler.html)
* [Intel Math Kernel Library](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl.html)
* [Microsoft Visual Studio (Community)](https://visualstudio.microsoft.com/de/vs/community/)
Solution and project files are contained in this fork repository, ready to compile.

## Contributing

Please contact [Dr. Hamish Brown](https://github.com/HamishGBrown) with project suggestions and contributions.

## Authors
* **Professor Les J Allen**
* [**Dr. Hamish G Brown**](https://github.com/HamishGBrown) - *Maintainer*
* **Dr. Adrian J D'Alfonso**
* **Dr. Scott D Findlay**
* **Dr. Ben D Forbes**
* [**Dr. Juri Barthel**](https://github.com/ju-bar) - *Fork Maintainer*


## License

This project is licensed under the GNU GPLv3.0 - see the [LICENSE.txt](LICENSE.txt) file for details.


## Acknowledgments

Les Allen, Adrian D'Alfonso and Scott Findlay originally took the initiative to make this code publicly available, with Adrian D'Alfonso taking on the considerable task of setting up the GPU code. Hamish Brown and Ben Forbes have subsequently made substantial refinements and extensions, with Ben Forbes responsible for several efficiency improvements. The code was developed mainly at the University of Melbourne. We would like to acknowledge the contributions of Mark Oxley and Chris Rossouw to the theoretical and numerical aspects of collaborative research that underpins μSTEM

In particular, the code IONIZER, whose primary author has been Mark Oxley, was used to set up the parametrized atomic scattering factors for ionization used in μSTEM. Eireann Cosgriff, Torgny Josefsson, Nathan Lugg, Andrew Martin, Gary Ruben and Chris Witte have been members of the group at Melbourne University at various stages and all made contributions to our understanding of inelastic scattering and the simulation thereof. 

