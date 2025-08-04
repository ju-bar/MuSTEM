!--------------------------------------------------------------------------------
!
!  Copyright (C) 2017  L. J. Allen, H. G. Brown, A. J. D’Alfonso, S.D. Findlay, 
!                      B. D. Forbes, J. Barthel
!
!  modified:
!  - J. Barthel, 2019-12-12 - added plasmonmc switch
!  - J. Barthel, 2025-06-26 - note on grid allocations
!
!  This program is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!  
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!   
!  You should have received a copy of the GNU General Public License
!  along with this program.  If not, see <http://www.gnu.org/licenses/>.
!                       
!--------------------------------------------------------------------------------
    
!--------------------------------------------------------------------------------
!
! Module for global variables
!
! Note on grid allocations (2025-06-26, JB):
!     Any array holding grid data is allocated with dimensions (nopiy,nopix)
!     where nopiy and nopix are the number of pixels in the y and x directions.
!     This deviates from the Fortran standard, which requires that the first
!     dimension is the fastest varying one thus storing such an array on disk
!     would have y changing along each row of data. This is however changed
!     upon output to disk by the subroutine binary_out in mod_output.f90, which
!     applies a transpose before writing to file.
!     In principle there is no good reason to have y as the fastest varying
!     dimension, but this is how the code was originally written. There is no
!     issue with this, as long as all array operations take that into account.
!     Also the CUDA code works optimally. There is just the small inconvenience
!     that the arrays are internally accessed as (y,x).
!
!--------------------------------------------------------------------------------

module global_variables
    
    use m_precision, only: fp_kind
    
    implicit none
    
    integer(4) :: nt,nm,i_xtl                   !number of atom types, number of atoms, xtl file read flag
    integer(4) :: nopiy,nopix,npixels           !supercell
    integer(4) :: nopix_ucell,nopiy_ucell       !unit cell
    integer(4) :: ifactorx,ifactory             !unit cell tilings
    real(fp_kind) :: deltay,deltax              !real space spacing between pixels in x and y
    real(fp_kind) :: normalisation
	
    
    real(fp_kind), allocatable :: bwl_mat(:,:)            !bandwidth limiting matrix
                                          
    integer(4), allocatable    :: nat(:)                  !number of each atom type in the unit celll
    real(fp_kind), allocatable :: dz(:)					  !ionicity of each atom type
    real(fp_kind), allocatable :: tau(:,:,:)              !position of the atoms in the unit cell
    real(fp_kind), allocatable :: atf(:,:)                !atomic number, occupancy and DWF (urms)
    real(fp_kind), allocatable :: atomf(:,:),fx(:)        !electron scattering factor parameterisation from elsa
    real(fp_kind)  :: a0(3),deg(3),ekv,ss(7)              !a b c unit cell lengths, angle between a b c, accelerating voltage, tricyclinc info
    real(fp_kind)  :: thetad,surfn(3),orthog(3,3)
    real(fp_kind)  :: volts,ak                            !mean inner potential, wavevector (corrected for refraction)
    real(fp_kind)  :: ak1,relm                            !wavevector in freespace, relativistically corrected mass
    real(fp_kind),allocatable :: claue(:,:), Kz(:)        !Specimen tilt vector and z component of incident wave vector
    integer*4::n_tilts_total                              !Total number of specimen tilts
    real(fp_kind)  :: bvec(3)                             !Beam tilt vector
    
    !sample thickness and slicing variables
    real(fp_kind) :: thickness                        !sample thickness
    real(fp_kind),allocatable:: zarray(:)
    integer(4),allocatable::ncells(:)
    integer(4)::nz
    integer(4) :: n_cells  
    logical::even_slicing

    complex(fp_kind), allocatable :: fz(:,:,:)            !the scattering factors, in reciprocal space, calculated on the grid (supercell)
    complex(fp_kind), allocatable :: fz_DWF(:,:,:)        !the DWF smear_array, in reciprocal space, calculated on the grid (supercell)
    !complex(fp_kind), allocatable :: sinc(:,:)            !sinc function to correct for pixelation in the potential construction
    complex(fp_kind), allocatable :: inverse_sinc(:,:)    !1/sinc function to correct for pixelation in the potential construction

    real(fp_kind)  :: uvw1(3),uvw2(3)                     !real space scan vectors that are parallel 
    integer(4) :: ig1(3),ig2(3),izone(3)                  !to the reciprocal space vectors ig1,ig2.
        
    character*120 :: substance                            !label for the crystal substance
    character*10, allocatable :: substance_atom_types(:)

    !output variables
    integer(4) :: ndet,nseg ,nopiyout,nopixout            !number of integrating detectors and 4D STEM output
    real(fp_kind) :: seg_det_offset
    logical::segments
    real(fp_kind), allocatable :: outer(:),inner(:)       !detector ranges (in inverse angstrom)
    
    !interpolation variables
    integer(4) :: output_nopiy,output_nopix               !output number of pixels in the interpolated output image
    integer(4) :: tilex,tiley                             !Interpolated image tiling
    real(fp_kind)  ::  bwl_rad                            !band width limiting radius (in q space)
                                                
    !Constants data
    real(fp_kind),parameter :: pi = atan(1.0_fp_kind)*4.0_fp_kind 
    real(fp_kind),parameter :: tp = atan(1.0_fp_kind)*8.0_fp_kind 
    real(fp_kind),parameter :: hsq_on_twom = 150.4132_fp_kind  !h^2/2m
    complex(fp_kind) :: ci = cmplx(0.0_fp_kind,1.0_fp_kind)
    
    real(fp_kind),parameter :: const1 = 9.7846113e-07_fp_kind       ! const1 = 1/(2mc^2) in eV^(-1)
    real(fp_kind),parameter :: const2 = 12.263868_fp_kind           ! const2 = h  in eV.A
    real(fp_kind),parameter :: bohr = 0.529177_fp_kind              ! bohr radius
    real(fp_kind),parameter :: ryd = 13.60535_fp_kind               ! rydberg constant
    real(fp_kind),parameter :: fsc = 7.29735e-3_fp_kind             ! fsc = fine structure constant (dimensionless)
    real(fp_kind),parameter :: hbarc = 1973.26_fp_kind              ! hbarc = hbar * c in eV A units
    
    
    !logical types to pick inelastic calculations
    logical :: adf 
    logical :: EELS = .false.
    logical :: EDX = .false.
    logical :: SEI = .false.
    
    ! SE variables, 2025-06-26, JB
    real(fp_kind) :: se_det_scale = 1.0_fp_kind ! scale factor for the SE detector, initialised to accept full angular range
    real(fp_kind) :: se_imfp = 0.001_fp_kind ! inelastic mean free path (1/e decay) for secondary electrons, initialised to 0.001 A
    
    ! inelastic filter thresholds, 2025-07-29, JB
    real(fp_kind) :: Hn0_filter_threshold = 0.001_fp_kind ! threshold for Hn0 filter, initialised to 0.001
    real(fp_kind) :: iwf_filter_threshold = 0.01_fp_kind ! threshold for inelastic wavefunctions, initialised to 0.01
    
    logical :: qep,output_thermal,interpolation,fourdSTEM
    
    logical :: on_the_fly = .false.
    logical :: high_accuracy
    logical :: ionic = .false.
    logical :: tp_eels,istem
    logical :: single_channeling = .true. ! tp_eels single-channeling mode flag, 2025-07-29, JB
    logical :: plasmonmc = .false.
    logical :: linpoleels = .false.       ! linear eels range interpolation (2025-May-14, JB)
    
    ! arguments for debugging
    integer(4) :: arg_debug_wave = 0    ! 0 = no debug, 1 = print debug values of wave function, triggered by secret option mmouse_wave
    integer(4) :: arg_debug_intens = 0  ! 0 = no debug, 1 = print debug values of probe intensity, triggered by secret option mmouse_intens
    integer(4) :: arg_debug_stemdet = 0 ! 0 = no debug, 1 = print debug values of stem detectors, triggered by secret option mmouse_stemdet
    integer(4) :: arg_debug_dump = 0    ! 0 = no debug dumps, bits of this can trigger particular memory dumps
    
	contains
     
    
    function wavev(e)
      !  this function returns the wavevector in one over lambda, in a-1,
      !  for an input electron energy e in ev.
      
      use m_precision
      
      implicit none
      
      real(fp_kind) c1,c2,e
      real(fp_kind) wavev
      data c1, c2 / 9.78475598e-07_fp_kind, 12.2642596_fp_kind /
      
      wavev = sqrt( e + c1 *e ** 2.0_fp_kind ) / c2
      
      end function 
    
      subroutine constants()
      
          use m_precision
      
          implicit none
      
      
          real(fp_kind) con1,c2,delk
          data con1, c2 / 510.9989461_fp_kind, 1.956951198e-03_fp_kind /

          relm = ( con1 + ekv) * c2
          
    ! ak is the incident wave vector in the solid corrected for refraction (mean inner potential)
          ak   = wavev( ekv * 1000_fp_kind + volts)
          
    ! ak1 is the incident wave vector without corrections for refraction
          ak1   = wavev( ekv * 1000_fp_kind )
    !Initialise tilt to zero (on-axis)
          allocate(Kz(1),claue(3,1))
          n_tilts_total = 1
          Kz = ak1
          claue = 0_fp_kind
          
          delk = ak - ak1

	    write(*,*) 'Pertinent quantities for the crystal:'
          write(*,111) ekv,ak,'A',ak1,'A',delk,'A',relm
      111 format('          E = ',F12.3,' keV',/,&
         &' Crystal Ko = ',G20.10,1x,a1,'-1 (1/lambda)',/,&
         &'  Vacuum Ko = ',G20.10,1x,a1,'-1 (1/lambda)',/,&
         &'   Delta Ko = ',G20.10,1x,a1,'-1 (1/lambda)',/,&
         &'     m* / m = ',g12.5,&
         &' (relativistic mass increment)',/,/)

      end subroutine 
      

      end module global_variables
