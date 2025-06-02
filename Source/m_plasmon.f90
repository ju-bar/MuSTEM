!--------------------------------------------------------------------------------
!
!  Copyright (C) 2019  J. Barthel, L. J. Allen
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
!
!  Plasmon scattering Monte-Carlo implementing the approach of 
!    B. Mendis, Ultramic. 206 (2019) 112816.
!    doi:10.1016/j.ultramic.2019.112816
!
!  Contains corrections and alternative code based on discussions
!  between J. Barthel, B. Mendis, S. Findlay, and  L.J. Allen
!  (Aug.-Oct. 2019)
!
!  This is a re-implementation from original code written by J. Barthel as
!  published under
!  https://github.com/ju-bar/drprobe_clt/blob/master/msa/plasmon.f90
!
!--------------------------------------------------------------------------------
!
!  How to use this module:
!
!  call setup_plasmon_parameters to work out plasmon scattering parameters
!  from use input and use the code as follows below.
!
!  1) Initialization
!     call pl_init_plasmon for initializing bulk plasmon scattering with
!     automatically determined max. number of excitation levels.
!     or
!     call pl_init_lowloss for initializing with a direct setup of
!     mean free path, scattering angles and max. number of excitations.
!
!  2) Preparation before a multislice passes
!     call pl_reset to rest the excitation levels. Do this before you start
!     the next multislice run.
!
!  3) Prepare for scattering events
!     Store slice thickness, grid sampling and random number generator seed
!     in module variable pl_tslc, pl_grid_sqx, pl_grid_sqy, pl_idum before
!     the next slice iteration of the multislice sequence
!
!  4) Check for plasmon scattering
!     When the wave function is in Fourier-space representation, before
!     applying the free-space propagator, call pl_scatt_slc to determine
!     if an inelastic scattering event occurs.
!     The function returns the number of excitations occurring in a slice
!     and the amount of Fourier space pixel shift to be applied.
!
!  5) Populate the Monte-Carlo excitation statistics
!     Call pl_populate right after the multislice run is finished and before
!     the next run is going to happen with identical probe position to fill
!     the Monte-Carlo statistics.
!     pl_mc_num = number of Monte-Carlo runs
!     pl_mc_exc_pop(0:pl_npemax) = frequency of excitations 
!--------------------------------------------------------------------------------

module plasmon

	use m_precision, only: fp_kind

	! declaration of module members
	implicit none
	! - module constant
	real(fp_kind),public,parameter :: pl_erest = 510.9989461_fp_kind ! electron rest energy [keV]
	real(fp_kind),public,parameter :: pl_wthr = 0.01_fp_kind ! probability threshold for neglecting higher number of excitations (>0 and <1)
	! - module parameters (variables)
	integer(4),public :: pl_npemax ! max. number of inelastic transitions per probing electron
	real(fp_kind),public :: pl_lmfp ! mean free path [Angs] (inverse strength of inelastic interaction)
	real(fp_kind),public :: pl_theta_crit ! critical angle [rad] (cut-off angle of scattering)
	real(fp_kind),public :: pl_theta_char ! characteristic angle [rad] (decay angle of scattering)
	real(fp_kind),public :: pl_theta_crit2 ! critical angle squared
	real(fp_kind),public :: pl_theta_char2 ! characteristic angle squared
	
	data pl_npemax /0/ ! initialize with inactive scattering
	data pl_lmfp, pl_theta_char, pl_theta_crit / 100.0_fp_kind, 0.00005_fp_kind, 0.02_fp_kind /
	data pl_theta_char2, pl_theta_crit2 / 2.5E-9_fp_kind, 0.0004_fp_kind /
	
	! Monte-Carlo volatile handlers used by calling routines to store plasmon scattering relevant parameters
	! RNG seed
	integer(4),public :: pl_idum
	data pl_idum /0/
	! grid Fourier sampling rate [rad/pixel]
	real(fp_kind),public :: pl_grid_sqx,pl_grid_sqy
	data pl_grid_sqx,pl_grid_sqy / 0.0_fp_kind, 0.0_fp_kind /
	! slice thickness [Angs]
	real(fp_kind),public :: pl_tslc
	data pl_tslc /0.0_fp_kind /
	
	! results of the Monte-Carlo for one run
	! /!\ Allocations sizes on last dimension: (0:pl_npemax) /!\
	!     corresponds to the number of plasmon excitations.
	integer(4),public :: pl_exc_num ! number of excitations for the run <= pl_npemax
	data pl_exc_num /0/
	real(fp_kind),public,allocatable :: pl_exc_dq(:,:) ! excitation scattering angle [rad]
	real(fp_kind),public,allocatable :: pl_exc_dqtot(:,:) ! accumulated scattering angle [rad]
	
	! statistics of the whole Monte-Carlo
	integer(4),public :: pl_mc_num ! number of MC runs
	data pl_mc_num /0/
	integer(4),public,allocatable :: pl_mc_exc_pop(:) ! population numbers
  
	! infrastructure
	integer*4, public :: pl_num_err ! number of errors
	data pl_num_err /0/ 
	character(len=2048),public :: pl_msg_err ! last error message
	data pl_msg_err /""/
	
	contains
	

	subroutine setup_plasmon_parameters(ekin, tmax, doplasm)
    
		use m_user_input, only: get_input
		use m_string, only: to_string, command_line_title_box

		implicit none

		real(fp_kind),intent(in)::ekin,tmax
		logical,intent(inout)::doplasm
		
		real(fp_kind)::lmfp,prm1,prm2
		integer(4) :: i, in_mode, nep
		call command_line_title_box('Plasmon / low-loss scattering parameters')
		
10		in_mode = 0 ! preset plasmon input mode
		doplasm = .false.
		i = 0
		write(*,*) 'Simulate bulk plasmon / low-loss inelastic scattering?'
		write(*,*) '<1>  Setup bulk plasmons.'
		write(*,*) '<2>  Setup generic low-loss.'
		write(*,*) '<0>  No, continue.'
		call get_input('Plasmon setup: <1> bulk plasmons <2> generic low-loss <0> none.', i)
		write(*,*)
		if (i<1 .or. i>2) then
			doplasm = .false.
			return
		end if
		in_mode = i
		
		write(*,*) 'Bulk plasmon / low-loss inelastic scattering will be simulated'
		write(*,*) 'by a Monte-Carlo approach [Mendis, Ultramic. 206 (2019) 112816].'
		write(*,*) 'A large amount of QEP passes is recommended with this option.'
		write(*,*)
		
		select case (in_mode)
		case (1)
			write(*,*) 'Setup of bulk plasmon parameters:'
			write(*,*) 'Enter the mean free path for single scattering in Angstroms:'
			call get_input('Mean free path for single scattering in Angstroms', lmfp)
			write(*,*)
			write(*,*) 'Enter the bulk plasmon energy in eV:'
			call get_input('Bulk plasmon energy [eV]', prm1)
			write(*,*)
			call pl_init_plasmon(ekin, abs(prm1), abs(lmfp), tmax, i)
			if (i==0) then
				doplasm = .true.
			else
				write(6,101) 'during plasmon scattering setup', i
				write(*,*) '  ',trim(pl_msg_err)
				write(*,*) 'Plasmon scattering deactivated.'
				write(*,*)
				doplasm = .false.
			end if
		case (2)
			write(*,*) 'Setup of approximative low-loss scattering parameters'
			write(*,*) 'Enter the mean free path for single scattering in Angstroms:'
			call get_input('Mean free path for single scattering in Angstroms:', lmfp)
			write(*,*)
			write(*,*) 'Enter the characteristic angle in mrad:'
			call get_input('Characteristic scattering angle in mrad:', prm1)
			write(*,*)
			write(*,*) 'Enter the critical angle in mrad:'
			call get_input('Critical scattering angle in mrad:', prm2)
			write(*,*)
			write(*,*) 'Enter the max. number of allowed excitations:'
			call get_input('Max. number of excitations per probing electron:', nep)
			write(*,*)
			call pl_init_lowloss(lmfp, prm1*0.001_fp_kind, prm2*0.001_fp_kind, nep, i)
			if (i==0) then
				doplasm = .true.
			else
				write(6,101) 'during plasmon scattering setup', i
				write(*,*) '  ',trim(pl_msg_err)
				write(*,*) 'Plasmon scattering deactivated.'
				write(*,*)
				doplasm = .false.
			end if
		end select
		! cross check user input
		if (doplasm) then
			write(*,*) 'Plasmon scattering setup'
			write(*,*) '----------------------------------------------'
			write(6,102) pl_lmfp, 'A'
			write(6,103) pl_theta_char * 1000.0_fp_kind
			write(6,104) pl_theta_crit * 1000.0_fp_kind
			write(6,105) pl_npemax
			write(*,*) '----------------------------------------------'
			write(*,*) '<0> Proceed without changes.'
			write(*,*) '<1> Repeat plasmon setup.'
			write(*,*)
		else
			write(*,*) 'Setup of plasmon scattering Monte-Carlo failed.'
			write(*,*) 'Do you want to repeat the setup?'
			write(*,*) '<0> Proceed without changes.'
			write(*,*) '<1> Repeat plasmon setup.'
			write(*,*)
		end if
		call get_input('Repeat plasmon setup? <0> No. <1> Yes.', i)
		if (i==1) goto 10
		return
101		format('Error ', a, ' (code: ', i4.4, ')')
102		format('   mean free path             | ', f8.1, ' ', a1)
103		format('   characteristic angle       | ', f8.5, ' mrad')
104		format('   critical angle             | ', f8.5, ' mrad')
105		format('   max. number of excitations | ', i8)
    end subroutine


	!***************************************************************************
	!
	! pl_init_plasmon
	!
	! initializes the module for plasmon Monte-Carlo using plasmon parameters
	!
	! input:
	! 	ekin : kinetic energy of the probong electrons [eV]
	!   eplasm : plasmon energy [eV]
	!   lmfp : mean free path for single plasmon excitation [A]
	!   tmax : maximum sample thickness to calculate [A]
	!
	! output:
	!   nerr : error code (0 = success)
	!
	subroutine pl_init_plasmon(ekin,eplasm,lmfp,tmax,nerr)
		implicit none
		! interface
		real(fp_kind),intent(in) :: ekin, eplasm, lmfp, tmax
		integer(4),intent(inout) :: nerr
		! locals
		real(fp_kind) :: e0, tol, wt, facn
		integer(4) :: nalloc
		
		nerr = 0
		nalloc = 0
		e0 = pl_erest * 1000.0_fp_kind
		call pl_deinit()
		if (eplasm<=1.) goto 102
		if (ekin<eplasm) goto 103
		if (lmfp<=0.) goto 105
		if (tmax<0.) goto 106
		
		! calculate derived parameters
		! - characteristic angle [rad] (Edgerton (2011) Electron Energy-Loss Spectroscopy in the Electron Microscope)
		pl_theta_char = eplasm / ekin * (ekin + e0) / (ekin + 2.0_fp_kind*e0)
		pl_theta_char2 = pl_theta_char * pl_theta_char
		! - critical angle [rad] ! This parameter is critical.
		!   Egerton Ultramic. 107 (2007) 575-586.: approximation
		pl_theta_crit = sqrt(eplasm / ekin) 
		pl_theta_crit2 = pl_theta_crit * pl_theta_crit
		! - mean free path and t over lambda
		pl_lmfp = lmfp
		tol = tmax/lmfp
		! max number number of registered plasmon excitations
		pl_npemax = 0
		wt = 1. - exp(-tol)
		facn = 1.
		do while (wt > pl_wthr)	! include more excitation levels until remaining
								! total probability of all higher levels is below
								! the user defined threshold pl_wthr
			pl_npemax = pl_npemax + 1
			facn = facn * real(pl_npemax) ! n!
			wt = wt - (tol**pl_npemax * exp(-tol) / facn)
		end do
		!
		! allocations
		allocate(pl_exc_dq(1:2,0:pl_npemax), stat=nalloc)
		if (nalloc/=0) goto 107
		allocate(PL_exc_dqtot(1:2,0:PL_npemax), stat=nalloc)
		if (nalloc/=0) goto 107
		allocate(PL_mc_exc_pop(0:PL_npemax), stat=nalloc)
		if (nalloc/=0) goto 107
		PL_exc_dq = 0.0
		PL_exc_dqtot = 0.0
		PL_mc_exc_pop = 0
    
100		return
    
102		nerr = 2
		pl_num_err = pl_num_err + 1
		pl_msg_err = trim(pl_msg_err)//" (init): (eplasm<1) too small plasmon energy"
		goto 100
103		nerr = 3
		pl_num_err = pl_num_err + 1
		pl_msg_err = trim(pl_msg_err)//" (init): (ekin<eplasm) too small electron energy"
		goto 100
105		nerr = 5
		pl_num_err = pl_num_err + 1
		pl_msg_err = trim(pl_msg_err)//" (init): (lmfp<=0) mean-free path not positive"
		goto 100
106		nerr = 6
		pl_num_err = pl_num_err + 1
		pl_msg_err = trim(pl_msg_err)//" (init): (tmax<0) negative sample thickness"
		goto 100
107		nerr = 7
		pl_num_err = pl_num_err + 1
		pl_msg_err = trim(pl_msg_err)//" (init): memory allocation failed"
		goto 100
	end subroutine pl_init_plasmon
	
	
	!***************************************************************************
	!
	! pl_init_lowloss
	!
	! initializes the module for los-loss Monte-Carlo
	!
	! input:
	!   lmfp : mean free path for single plasmon excitation [A]
	!   theta_char : characteristic angle [rad]
	!   theta_crit : critical angle [rad]
	!   npemax : number of max. excitations per probing electron
	!
	! output:
	!   nerr : error code (0 = success)
	!
	subroutine pl_init_lowloss(lmfp,theta_char,theta_crit,npemax,nerr)
		implicit none
		! interface
		real(fp_kind),intent(in) :: lmfp,theta_char,theta_crit
		integer(4),intent(in) :: npemax
		integer(4),intent(inout) :: nerr
		! locals
		integer(4) :: nalloc
		real(fp_kind) :: tol
		
		nerr = 0
		nalloc = 0
		pl_num_err = 0
		pl_msg_err = ""
		pl_exc_num = 0
		pl_mc_num = 0
		if (allocated(pl_exc_dq)) deallocate(pl_exc_dq, stat=nalloc)
		if (allocated(pl_exc_dqtot)) deallocate(pl_exc_dqtot, stat=nalloc)
		if (allocated(pl_mc_exc_pop)) deallocate(pl_mc_exc_pop, stat=nalloc)
		if (theta_char<=0.) goto 102
		if (theta_crit<=theta_char) goto 103
		if (npemax<=0) goto 104
		if (lmfp<=0.) goto 105
		pl_lmfp = lmfp
		pl_npemax = npemax
		! calculate derived parameters
		pl_theta_char = theta_char
		pl_theta_char2 = pl_theta_char * pl_theta_char
		pl_theta_crit = theta_crit
		pl_theta_crit2 = pl_theta_crit * pl_theta_crit
		! allocations
		allocate(pl_exc_dq(1:2,0:pl_npemax), stat=nalloc)
		if (nalloc/=0) goto 107
		allocate(pl_exc_dqtot(1:2,0:pl_npemax), stat=nalloc)
		if (nalloc/=0) goto 107
		allocate(pl_mc_exc_pop(0:pl_npemax), stat=nalloc)
		if (nalloc/=0) goto 107
		pl_exc_dq = 0.0
		pl_exc_dqtot = 0.0
		pl_mc_exc_pop = 0
100		return
102		nerr = 2
		pl_num_err = pl_num_err + 1
		pl_msg_err = trim(pl_msg_err)//" (init): (theta_char<=0) too small characteristic angle"
		goto 100
103		nerr = 3
		pl_num_err = pl_num_err + 1
		pl_msg_err = trim(pl_msg_err)//" (init): (theta_crit<theta_char) too small critical angle"
		goto 100
104		nerr = 4
		pl_num_err = pl_num_err + 1
		pl_msg_err = trim(pl_msg_err)//" (init): (npemax<=0) too low number of allowed excitations"
		goto 100
105		nerr = 5
		pl_num_err = pl_num_err + 1
		pl_msg_err = trim(pl_msg_err)//" (init): (lmfp<=0) mean-free path not positive"
		goto 100
107		nerr = 7
		pl_num_err = pl_num_err + 1
		pl_msg_err = trim(pl_msg_err)//" (init): memory allocation failed"
		goto 100
	end subroutine pl_init_lowloss
	
	
	subroutine pl_deinit()
		implicit none
		integer(4) nalloc
		
		nalloc = 0
		pl_num_err = 0
		pl_msg_err = ""
		pl_npemax = 0
		pl_exc_num = 0
		pl_mc_num = 0
		if (allocated(pl_exc_dq)) deallocate(pl_exc_dq, stat=nalloc)
		if (allocated(pl_exc_dqtot)) deallocate(pl_exc_dqtot, stat=nalloc)
		if (allocated(pl_mc_exc_pop)) deallocate(pl_mc_exc_pop, stat=nalloc)
		return
	end subroutine pl_deinit
	
	
	!*******************************************************************!
	!
	! pl_mc_slc
	!
	! input:
	!     dt : slice thickness [A]
	!     iexc_max : max. allowed excitations
	!
	! output:
	!     iexc : excitation flag (0 = no, >0: number of excitations)
	!     dqx, dqy : scattering angle [rad]
	!
	! in/output:
	!     idum : rng seed transfer
	!
	! purpose: calculate whether inelastic scattering takes place
	!          in a slice, and if so, determines the scattering angle.
	!
	! results: None of the other module variables are changed. These
	!          changes need to be applied by the calling routine.
	!
	subroutine pl_mc_slc(dt, iexc_max, iexc, dqx, dqy, idum)
		use global_variables
		use m_numerical_tools
  		implicit none
		! interface
		real(fp_kind),intent(in) :: dt ! slice thickness [A]
		integer(4),intent(in) :: iexc_max ! max. allowed excitations
		integer(4),intent(out) :: iexc ! excitation flag
		real(fp_kind),intent(out) :: dqx, dqy ! scattering angle [rad]
		integer(4),intent(inout) :: idum
		! locals
		integer(4) :: i, ipr1
		real(fp_kind) :: pcur, qcur, r1, r2, qpowr
		
		iexc = 0
		dqx = 0.0_fp_kind
		dqy = 0.0_fp_kind
		! calculate a Poissonian random number (limited to max. allowed excitation level)
		ipr1 = poirand1(abs(dt)/pl_lmfp, idum)
		iexc = min(ipr1, iexc_max)
		if (iexc > 0) then ! excitation happens
			do i=1, iexc
				r1 = unirand(idum)
				r2 = unirand(idum)
				pcur = tp * r1 ! uniform RNG phi
				qpowr = (pl_theta_crit2 / pl_theta_char2 + 1.)**r2
				qcur = sqrt( pl_theta_char2 * qpowr - pl_theta_char2) ! Lorentzian Theta
				dqx = dqx + qcur * cos(pcur)
				dqy = dqy + qcur * sin(pcur)
			end do
		end if
100 	return
	end subroutine PL_mc_slc
	
	
	
	!*******************************************************************!
	!
	! pl_reset
	!
	! input:
	!     none
	!
	! output:
	!     none
	!
	! purpose: resets the housholder arrays keeping track of tilts and
	!          number of excitations per run through the specimen
	!
	subroutine pl_reset()
		implicit none
		pl_exc_num = 0
		pl_exc_dq = 0.0_fp_kind
		pl_exc_dqtot = 0.0_fp_kind
		return
	end subroutine pl_reset
	
	
	
	!*******************************************************************!
	!
	! pl_populate
	!
	! input:
	!     integer*4 :: num_exc ! number of plasmon excitations to register
	!                          ! for the finished run in the population
	!                          ! statistics
	!
	! output:
	!     none
	!
	! purpose: updates the population statistics in plasmon-loss channels
	!          call this when a multislice run is finished
	!
	subroutine pl_populate(num_exc)
		implicit none
		integer(4), intent(in) :: num_exc
		integer(4) :: n
		n = max(0, min( pl_npemax, num_exc) ) ! plasmon-loss channel 0 .. PL_npemax
		pl_mc_exc_pop(n) = PL_mc_exc_pop(n) + 1 ! increment excitation level count
		pl_mc_num = pl_mc_num + 1 ! increment number of MC passes
		return
	end subroutine pl_populate
	
	
	!*******************************************************************!
	!
	! pl_scatt_slc
	!
	! input:
	!     dt = current slice thickness [Angs]
	!     sqx, sqy = sampling rates of the wave function [rad/pixel]
	!
	! output:
	!     iexc = number of excitations triggered in the slice
	!     isx, isy = number of Fourier pixel shifts to apply
	!
	! in/output:
	!     idum = random number seed transfer
	!
	! purpose: Monte-Carlo plasmon scattering call for a slice
	!     returns number of excitations and possible Fourier shift
	!     to apply to the wave function (to do by the calling process)
	!
	subroutine pl_scatt_slc(dt, sqx, sqy, iexc, isx, isy, idum)
		implicit none
		! interface
		real(fp_kind), intent(in) :: dt ! slice thickness [Angs]
		real(fp_kind), intent(in) :: sqx, sqy ! sampling rates of the wave function [rad/pixel]
		integer(4), intent(out) :: iexc ! number of excitations happening (return)
		integer(4), intent(out) :: isx, isy ! number of pixels to shift the wave function
		integer(4), intent(inout) :: idum
		! locals
		integer(4) :: iexc_max, nexc ! number of excitations
		integer(4) :: isx0, isy0, isx1, isy1, isxd, isyd ! pixel shifts
		integer(4) :: i, j, i1, j1
		real(fp_kind) :: dqx, dqy

		iexc = 0
		iexc_max = pl_npemax - pl_exc_num ! max. further allowed excitations
		isx = 0
		isy = 0
		dqx = 0.0_fp_kind
		dqy = 0.0_fp_kind
		! decide if excitations happen in this slice
		call pl_mc_slc(dt, iexc_max, iexc, dqx, dqy, idum)
		if (iexc > 0) then ! excitations have happened
			! new excitation level
			nexc = pl_exc_num + iexc
			! set the change of scattering angle with this excitation
			pl_exc_dq(1, nexc ) = dqx
			pl_exc_dq(2, nexc ) = dqy
			! set the total scattering angle with this excitation
			do i=pl_exc_num+1, nexc ! loop in case that more than 1 excitation happened (iexc>=1)
				pl_exc_dqtot(1, i) = pl_exc_dqtot(1, pl_exc_num ) + dqx
				pl_exc_dqtot(2, i) = pl_exc_dqtot(2, pl_exc_num ) + dqy
			end do
			! get current pixel shift ( dqtot old -> pixels )
			isx0 = nint( pl_exc_dqtot(1, pl_exc_num) / sqx )
			isy0 = nint( pl_exc_dqtot(2, pl_exc_num) / sqy )
			! get new pixel shift ( dqtot new -> pixels )
			isx1 = nint( pl_exc_dqtot(1, nexc) / sqx )
			isy1 = nint( pl_exc_dqtot(2, nexc) / sqy )
			! get change of pixel shift (pixel new - pixel old)
			isx = isx1 - isx0
			isy = isy1 - isy0
			! update the number of excitations for this run
			pl_exc_num = pl_exc_num + iexc
		end if
		return
	end subroutine pl_scatt_slc
	

end module plasmon