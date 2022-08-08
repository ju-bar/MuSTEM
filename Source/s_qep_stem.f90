!--------------------------------------------------------------------------------
!
!  Copyright (C) 2017  L. J. Allen, H. G. Brown, A. J. D’Alfonso, S.D. Findlay, B. D. Forbes
!
!  modiefied:
!    2019-Dec-13, J. Barthel, added setup and use of plasmon scattering code
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

function qep_stem_GPU_memory(n_qep_grates, quick_shift, phase_ramp_shift) result(required_memory)
    
    use m_precision, only: fp_kind
    use global_variables, only: nopiy, nopix, ifactory, ifactorx, on_the_fly, ndet
    use m_lens, only: imaging
    use m_multislice
    
    implicit none
    
    logical,intent(in)::quick_shift,phase_ramp_shift
    integer*4,intent(in)::n_qep_grates
    
    real(fp_kind) :: required_memory
    
    real(fp_kind) :: array_size
    integer :: array_count
    
    array_size = 4.0_fp_kind * nopiy * nopix
    
    array_count = 2*(4 + n_slices + n_slices*n_qep_grates) + 2 + ndet
        
    if (on_the_fly.or.quick_shift.or.phase_ramp_shift) array_count = array_count + 2
    if (phase_ramp_shift) array_count = array_count + 2
   ! if (ionization) array_count = array_count + 1 + n_slices
    !if (eels) array_count = array_count + 1
    
    ! Temporary array used in cuda_stem_detector()
    array_count = array_count + 1
    
    required_memory = array_count * array_size
    
    if (phase_ramp_shift) required_memory = required_memory + 8.0_fp_kind*(nopiy*ifactory + nopix*ifactorx)
    
end function

    

subroutine qep_stem(STEM,ionization,PACBED)
    
    use global_variables
    use m_lens
    use m_user_input
    use m_precision
    use output
	use plasmon ! plasmon scattering module, only for main qep wave, not for double channeling, inserted JB-191213
#ifdef GPU	
    use cuda_array_library;use cudafor;use cuda_ms;use CUFFT;use cuda_potential;use cuda_setup
#endif
	use cufft_wrapper;use m_Hn0
	use m_crystallography;use m_tilt;use m_string;use m_multislice;use m_potential;use m_numerical_tools
    
    implicit none
    
	logical,intent(in)::STEM,ionization,PACBED
    !dummy variables
    integer(4) ::  i,j,l,m,i_qep_pass,iz,k,ii,jj,ntilt,count,shifty,shiftx,ny,nx,i_df,idet,total_slices,lengthimdf,starting_slice

    !random variables
    integer(4) :: idum,nsliceoutput,z_indx(1)
    !real(fp_kind) :: ran1
    complex(fp_kind) :: temp_transf(nopiy,nopix),prop(nopiy,nopix,n_slices)
       
    !probe variables
    complex(fp_kind),dimension(nopiy,nopix) :: psi,trans,psi_initial,psi_out,temp_image
     complex(fp_kind),dimension(nopiy,nopix,nz)::psi_elastic
     complex(fp_kind),dimension(nopiy,nopix,n_qep_grates,n_slices)::projected_potential,qep_grates
	 complex(fp_kind),allocatable::ctf(:,:,:)
    !output
    real(fp_kind),dimension(nopiy,nopix) :: image,temp
    
    !STEM image variables
    real(fp_kind) :: masks(nopiy,nopix,ndet),cbed(nopiy,nopix,nz)               !detector masks
    real(fp_kind),dimension(nysample,nxsample,probe_ndf,ndet,nz) :: stem_image,stem_elastic_image&
																  &,stem_inelastic_image
	real(fp_kind),dimension(nysample,nxsample,probe_ndf,nz) :: eels_correction_image
    real(fp_kind),allocatable :: ion_image(:,:,:),stem_ion_image(:,:,:,:,:),pacbed_pattern(:,:,:),pacbed_elastic(:,:,:),istem_image(:,:,:,:)
    
    !diagnostic variables
    real(fp_kind) :: intensity,intens
! intens added May 2019 to accommodate changes awhere it occurs below 
    real(fp_kind) :: t1, delta
    
    !output variables
    character(120) :: fnam, fnam_temp, fnam_det
#ifdef GPU
    !device variables
	integer :: plan
    complex(fp_kind),device,dimension(nopiy,nopix) :: psi_d,psi_out_d,psi_initial_d,trans_d
	complex(fp_kind),device,dimension(nopiy,nopix,nz)::psi_elastic_d
    real(fp_kind),device :: cbed_d(nopiy,nopix,nz),temp_d(nopiy,nopix)
	complex(fp_kind),device,allocatable :: prop_d(:,:,:),transf_d(:,:,:,:)
	complex(fp_kind),device,allocatable,dimension(:,:) ::shift_arrayx_d,shift_arrayy_d,shift_array_d
    real(fp_kind),device,allocatable,dimension(:,:,:) :: masks_d,ion_image_d,pacbed_pattern_d
    real(fp_kind),device,allocatable :: eels_correction_detector_d(:,:),ion_potential_d(:,:,:,:)
    
    !device variables for on the fly potentials
    complex(fp_kind),device,allocatable,dimension(:,:) :: bwl_mat_d,inverse_sinc_d,fz_dwf_d
    complex(fp_kind),device,allocatable,dimension(:,:,:) :: fz_d,fz_mu_d
    real(fp_kind),device,allocatable :: inelastic_potential_d(:,:)

	!Double channeling variables
	complex(fp_kind),device,allocatable,dimension(:,:) ::psi_inel_d,shiftarray,tmatrix_d,q_tmatrix_d
	complex(fp_kind),device,allocatable,dimension(:,:,:)::tmatrix_states_d,Hn0_shifty_coord_d,Hn0_shiftx_coord_d,ctf_d
	real(fp_kind),device,allocatable,dimension(:,:)::cbed_inel_dc_d
	real(fp_kind),device,allocatable,dimension(:,:,:)::Hn0_eels_detector_d
	real(fp_kind),device,allocatable,dimension(:,:,:,:)::efistem_image_d,istem_image_d
	real(fp_kind),allocatable,dimension(:,:,:,:,:)::Hn0_eels_dc
	real(fp_kind),allocatable,dimension(:,:,:)::tmatrix_states
	integer::i_target
#endif

    real(fp_kind),allocatable :: probe_intensity(:,:,:)
    character(1024) :: filename
    
    real(fp_kind)::qep_stem_GPU_memory,const
    logical::elfourd,manytilt,manyz,many_df
    integer*4::length,lengthdf
    
	
#ifdef GPU    
    call GPU_memory_message(qep_stem_GPU_memory(n_qep_grates, quick_shift, phase_ramp_shift), on_the_fly)
#endif
    
    
    call command_line_title_box('Pre-calculation setup')

    if(pacbed) then
        allocate(pacbed_pattern(nopiy,nopix,nz))
        pacbed_pattern=0
        if(output_thermal) then
            allocate(pacbed_elastic(nopiy,nopix,nz))
            pacbed_elastic = 0
        endif
	    elfourd = fourdSTEM.and.output_thermal
        length = ceiling(log10(maxval(zarray)))
    endif

#ifdef GPU
	if(double_channeling) then
        allocate(tmatrix_states_d(nopiy,nopix,nstates),psi_inel_d(nopiy,nopix),cbed_inel_dc_d(nopiy,nopix),tmatrix_states(nopiy,nopix,nstates))
		allocate(shiftarray(nopiy,nopix),tmatrix_d(nopiy,nopix),q_tmatrix_d(nopiy,nopix))
        tmatrix_states_d = setup_ms_hn0_tmatrices(nopiy,nopix,nstates)*alpha_n
        allocate(Hn0_shifty_coord_d(nopiy,maxval(natoms_slice_total),n_slices))
        allocate(Hn0_shiftx_coord_d(nopix,maxval(natoms_slice_total),n_slices))
        Hn0_shiftx_coord_d = Hn0_shiftx_coord
        Hn0_shifty_coord_d = Hn0_shifty_coord
		allocate(Hn0_eels_dc(nysample,nxsample,probe_ndf,nz,numeels))
		allocate(Hn0_eels_detector_d(nopiy,nopix,numeels),Hn0_eels_detector(nopiy,nopix,numeels))
		do l=1,numeels
			Hn0_eels_detector(:,:,l) = make_detector(nopiy,nopix,ifactory,ifactorx,ss,0.0_fp_kind,outerrad(l))
		enddo
		Hn0_eels_detector_d = Hn0_eels_detector
		Hn0_eels_dc = 0.0_fp_kind
		if(istem) then
			allocate(efistem_image_d(nopiy,nopix,imaging_ndf,nz))
			efistem_image_d=0
		endif
    endif
	if(istem) then
		allocate(istem_image_d(nopiy,nopix,imaging_ndf,nz))
		istem_image_d=0
	endif
#endif
	if(istem) then
		allocate(ctf(nopiy,nopix,imaging_ndf),istem_image(nopiy,nopix,imaging_ndf,nz))
		istem_image=0
		lengthimdf = calculate_padded_string_length(imaging_df,imaging_ndf)
		do i=1,imaging_ndf
			ctf(:,:,i) =  make_ctf([0.0_fp_kind,0.0_fp_kind,0.0_fp_kind],imaging_df(i),imaging_cutoff,imaging_aberrations,imaging_apodisation)
		enddo
#ifdef GPU
		allocate(ctf_d(nopiy,nopix,imaging_ndf),istem_image_d(nopiy,nopix,imaging_ndf,nz));ctf_d=ctf;istem_image_d=0
#endif
	endif
	manyz = nz>1
	manytilt = n_tilts_total>1
	many_df = probe_ndf .gt. 1
    length = calculate_padded_string_length(zarray,nz)
    lengthdf = calculate_padded_string_length(probe_df,probe_ndf)
    
! Precalculate the scattering factors on a grid
    call precalculate_scattering_factors()
    idum = seed_rng()
#ifdef GPU	
    if (on_the_fly) then
        call cuda_setup_many_phasegrate()               !setup the atom co-ordinate for on the fly potentials (the arguments are simply so that
	else
#else
    !Generally not practical for CPU calculations
    on_the_fly = .false.
    if (.not.on_the_fly) then
#endif
        projected_potential = make_qep_grates(idum)
        if(.not. load_grates) call load_save_add_grates(idum,projected_potential,nopiy,nopix,n_qep_grates,n_slices,nt,nat_slice)
        call make_local_inelastic_potentials(ionization)          !setup the REAL space inelastic potential (ionization and adf) for QUEP ADF is disabled
    endif

	if(stem) then
		do i=1,ndet/nseg
			do j=1,nseg
				if(nseg>1) masks(:,:,(i-1)*nseg+j) = make_detector(nopiy,nopix,ifactory,ifactorx,ss,inner(i),outer(i),2*pi*j/nseg-seg_det_offset,2*pi/nseg)
				if(nseg==1) masks(:,:,(i-1)*nseg+j) = make_detector(nopiy,nopix,ifactory,ifactorx,ss,inner(i),outer(i))
			enddo
		enddo
	endif

	t1 = secnds(0.0)    

    call command_line_title_box('Calculation running')

    if(stem) stem_image = 0.0_fp_kind
    if(stem) stem_elastic_image = 0.0_fp_kind
#ifdef GPU
	! Plan the fourier transforms
    if(fp_kind.eq.8)then
	    call cufftPlan(plan,nopix,nopiy,CUFFT_Z2Z)
	else
        call cufftPlan(plan,nopix,nopiy,CUFFT_C2C)
    endif
    ! Copy host arrays to the device
    if (stem.and.(ndet.gt.0)) then
        allocate(masks_d(nopiy,nopix,ndet))
        masks_d = masks
    endif

    if (ionization) then
		allocate(ion_image_d(nopiy,nopix,num_ionizations))
		if(.not.EDX) then
			allocate(eels_correction_detector_d(nopiy,nopix))
			eels_correction_detector_d=eels_correction_detector
		endif
    endif
    if(pacbed) then
		allocate(pacbed_pattern_d(nopiy,nopix,nz))
		pacbed_pattern_d=0
	endif
    allocate(prop_d(nopiy,nopix,n_slices))
    prop_d = prop
    if (on_the_fly) then
        allocate(bwl_mat_d(nopiy,nopix),inverse_sinc_d(nopiy,nopix),fz_d(nopiy,nopix,nt))
        fz_d = fz 
        if(ionization) then
            allocate(inelastic_potential_d(nopiy,nopix))
            allocate(fz_mu_d(nopiy,nopix,num_ionizations))
            fz_mu_d = ionization_mu
        endif
        inverse_sinc_d = inverse_sinc
        bwl_mat_d = bwl_mat
    else
        allocate(transf_d(nopiy,nopix,n_qep_grates,n_slices))
  	    transf_d = qep_grates
        if (ionization) then
            allocate(ion_potential_d(nopiy,nopix,num_ionizations,n_slices))
            ion_potential_d = ionization_potential
        endif
        if (qep_mode==3) then
            allocate(shift_array_d(nopiy,nopix))
            allocate(shift_arrayx_d(nopix,ifactorx))
            allocate(shift_arrayy_d(nopiy,ifactory))
            shift_arrayx_d = shift_arrayx
            shift_arrayy_d = shift_arrayy
        endif
    endif
#else
	if (ionization) allocate(ion_image(nopiy,nopix,num_ionizations)) 
#endif
    if (ionization) then
        allocate(stem_ion_image(nysample,nxsample,probe_ndf,nz,num_ionizations))
		stem_ion_image = 0.0_fp_kind
		if (.not.EDX) eels_correction_image = 0.0_fp_kind
        !call binary_out_unwrap(nopiy,nopix,eels_correction_detector,'eels_correction_detector')
	endif
    if (output_probe_intensity) allocate(probe_intensity(nopiy,nopix,size(output_thickness_list)))
    intensity = 1.0d0

    do ntilt=1,n_tilts_total
 	    do i = 1, n_slices
	        call make_propagator(nopiy,nopix,prop(:,:,i),prop_distance(i),Kz(1),ss,ig1,ig2,claue(:,1),ifactorx,ifactory)
	        prop(:,:,i) = prop(:,:,i) * bwl_mat
            qep_grates(:,:,:,i) = exp(ci*pi*a0_slice(3,i)/Kz(ntilt)*projected_potential(:,:,:,i))
			do j=1,n_qep_grates
			call fft2(nopiy,nopix,qep_grates(:,:,j,i),nopiy,psi,nopiy)
			if(qep_mode.eq.3) qep_grates(:,:,j,i)= psi*bwl_mat
			if(qep_mode.ne.3) call ifft2(nopiy,nopix,psi*bwl_mat,nopiy,qep_grates(:,:,j,i),nopiy)
		enddo
	enddo
#ifdef GPU
	prop_d=prop
	if(.not.on_the_fly) transf_d = qep_grates
#endif        
	lengthdf = ceiling(log10(maxval(abs(probe_df))))
	if(any(probe_df<0)) lengthdf = lengthdf+1
	do i_df = 1, probe_ndf
	do ny = 1, nysample ! code formatting relaxed due to too many loops
	do nx = 1, nxsample
#ifdef GPU
		write(6,903,ADVANCE='NO') achar(13), i_df, probe_ndf, ny, nysample, nx, nxsample, intensity
903	format(a1,' df:',i3,'/',i3,' y:',i3,'/',i3,' x:',i3,'/',i3,'  Intensity:', f6.3, ' (to monitor BWL)')	
#else
		write(6,900) i_df, probe_ndf, ny, nysample, nx, nxsample, intensity
900	format(1h+,1x,' df:',i3,'/',i3,' y:',i3,'/',i3,' x:',i3,'/',i3,'  Intensity:', f6.3, ' (to monitor BWL)')	
#endif        
!
!       Make STEM probe
!
!		psi_initial = make_ctf(probe_positions(:,ny,nx),probe_df(i_df),probe_cutoff,probe_aberrations,probe_apodisation)
!       call ifft2(nopiy, nopix, psi_initial, nopiy, psi_initial, nopiy)
!       psi_initial = psi_initial/sqrt(sum(abs(psi_initial)**2))
! The above three lines replaced by the following four to improve accuracy (May 2019) - Findlay suggestion
		psi_initial = make_ctf(probe_positions(:,ny,nx),probe_df(i_df),probe_cutoff,probe_aberrations,probe_apodisation)
		intens = sum(abs(psi_initial)**2)
		call ifft2(nopiy, nopix, psi_initial, nopiy, psi_initial, nopiy)
		psi_initial = psi_initial/sqrt(intens)

		call tilt_wave_function(psi_initial)
		if (output_probe_intensity) probe_intensity = 0.0_fp_kind
		cbed=0_fp_kind

		
#ifdef GPU
		cbed_d=0.0_fp_kind
		psi_elastic_d=0.0_fp_kind
		psi_initial_d = psi_initial
        do i_qep_pass = 1, n_qep_passes
		! Reset wavefunction
		psi_d = psi_initial_d
		if (ionization) ion_image_d=0.0_fp_kind

		! plasmon scattering reset, JB-191213
		if (plasmonmc) call pl_reset()

		do i = 1,maxval(ncells) ! main multislice loop of cells ! code formatting relaxed due to too many loops
		do j = 1, n_slices ! sub loop over slices per cell
			! Accumulate ionization cross section
			if(ionization) then
				do ii=1,num_ionizations
					call cuda_mod<<<blocks,threads>>>(psi_d,temp_d,1.0_fp_kind,nopiy,nopix) 
					if(on_the_fly) then
						call cuda_make_ion_potential(inelastic_potential_d,tau_slice(:,atm_indices(ii),:,j),nat_slice(atm_indices(ii),j),plan,&
														&fz_mu_d(:,:,ii),inverse_sinc_d,Volume_array(j))
						call cuda_multiplication<<<blocks,threads>>>(temp_d,inelastic_potential_d, temp_d,prop_distance(j),nopiy,nopix)
					else
						call cuda_multiplication<<<blocks,threads>>>(temp_d,ion_potential_d(:,:,ii,j), temp_d,prop_distance(j),nopiy,nopix)     !overlap
					endif
					call cuda_addition<<<blocks,threads>>>(temp_d,ion_image_d(:,:,ii),ion_image_d(:,:,ii),1.0_fp_kind,nopiy,nopix)                          !depth sum
				enddo ! ii
			endif

!Double channeling

			if(double_channeling) then
                    
                do i_target = 1, natoms_slice_total(j) ! Loop over targets
							
                    ! Calculate inelastic transmission matrix
					call cuda_make_shift_array<<<blocks,threads>>>(shiftarray,Hn0_shifty_coord_d(:,i_target,j),Hn0_shiftx_coord_d(:,i_target,j),nopiy,nopix)
					
					do k = 1, nstates
							
					call cuda_multiplication<<<blocks,threads>>>(tmatrix_states_d(:,:,k),shiftarray, q_tmatrix_d,1.0_fp_kind,nopiy,nopix)
                    call cufftExec(plan,q_tmatrix_d,tmatrix_d,CUFFT_INVERSE)
					call cuda_multiplication<<<blocks,threads>>>(psi_d,tmatrix_d,psi_inel_d,sqrt(normalisation),nopiy,nopix)
					starting_slice = j
					
					! Scatter the inelastic wave through the remaining cells
					do ii = i, n_cells
						do jj = starting_slice, n_slices
							! QEP multislice
							nran = floor(n_qep_grates*ran1(idum)) + 1
							shiftx = floor(ifactorx*ran1(idum));shifty = floor(ifactory*ran1(idum))
							if(on_the_fly) then
								call cuda_fph_make_potential(trans_d,ccd_slice_array(jj),tau_slice,nat_slice(:,jj),jj,prop_distance(jj),idum,plan,fz_d,inverse_sinc_d,bwl_mat_d)
								call cuda_multislice_iteration(psi_inel_d, trans_d, prop_d(:,:,jj), normalisation, nopiy, nopix,plan)
							elseif(qep_mode == 2) then
								call cuda_multislice_iteration(psi_inel_d, transf_d(:,:,nran,j), prop_d(:,:,j), normalisation, nopiy, nopix,shifty*nopiy_ucell,shiftx* nopix_ucell,plan)
							elseif(qep_mode == 3) then                       !randomly shift phase grate
								call cuda_phase_shift_from_1d_factor_arrays(transf_d(:,:,nran,jj),trans_d,shift_arrayy_d(:,shifty+1),shift_arrayx_d(:,shiftx+1),nopiy,nopix,plan)
								call cuda_multislice_iteration(psi_inel_d, trans_d, prop_d(:,:,jj), normalisation, nopiy, nopix,plan)
							else
								call cuda_multislice_iteration(psi_inel_d, transf_d(:,:,nran,jj), prop_d(:,:,jj), normalisation, nopiy, nopix,plan)
							endif
						enddo ! jj
						starting_slice=1
						if (any(ii==ncells)) then
							z_indx = minloc(abs(ncells-ii))
							call cufftExec(plan, psi_inel_d, psi_out_d, CUFFT_FORWARD)
							
							! Accumulate the EELS images
							do l=1,numeels
								Hn0_eels_dc(ny,nx,i_df,z_indx(1),l) = Hn0_eels_dc(ny,nx,i_df,z_indx(1),l)+cuda_stem_detector(psi_out_d, Hn0_eels_detector_d(:,:,l))
							enddo ! l
							! Accumulate EFISTEM images
							if (istem.and.i_df==1) then
								do l = 1, imaging_ndf
									call cuda_image(psi_out_d,ctf_d(:,:,l),temp_d,normalisation, nopiy, nopix,plan,.false.)
									call cuda_addition<<<blocks,threads>>>(efistem_image_d(:,:,l,z_indx(1)), temp_d, efistem_image_d(:,:,l,z_indx(1)), 1.0_fp_kind, nopiy, nopix)
								enddo
							endif
							!stop
						endif
					enddo ! ii
					enddo ! k
				enddo ! i_target
			endif  ! End loop over cells,targets and states and end double_channeling section

			! QEP multislice
			nran = floor(n_qep_grates*ran1(idum)) + 1
			shiftx = floor(ifactorx*ran1(idum));shifty = floor(ifactory*ran1(idum))
			! plasmon, inserted by JB-191213
			if (plasmonmc) then ! use plasmon code
				pl_tslc = prop_distance(j) ! store current slice thickness
				pl_grid_sqx = 1.0_fp_kind / a0_slice(1,j) / ak1 ! store reciprocal grid sampling along x in rad/pixel
				pl_grid_sqy = 1.0_fp_kind / a0_slice(2,j) / ak1 ! store reciprocal grid sampling along y in rad/pixel
				pl_idum = idum ! let plasmon module take over the random seed
			end if
			if(on_the_fly) then
				call cuda_fph_make_potential(trans_d,ccd_slice_array(j),tau_slice,nat_slice(:,j),j,prop_distance(j),idum,plan,fz_d,inverse_sinc_d,bwl_mat_d)
				call cuda_multislice_iteration(psi_d, trans_d, prop_d(:,:,j), normalisation, nopiy, nopix,plan)
			elseif(qep_mode == 2) then
				call cuda_multislice_iteration(psi_d, transf_d(:,:,nran,j), prop_d(:,:,j), normalisation, nopiy, nopix,shifty*nopiy_ucell,shiftx*nopix_ucell,plan)
			elseif(qep_mode == 3) then                       !randomly shift phase grate
				call cuda_phase_shift_from_1d_factor_arrays(transf_d(:,:,nran,j),trans_d,shift_arrayy_d(:,shifty+1),shift_arrayx_d(:,shiftx+1),nopiy,nopix,plan)
				call cuda_multislice_iteration(psi_d, trans_d, prop_d(:,:,j), normalisation, nopiy, nopix,plan)
			else
				call cuda_multislice_iteration(psi_d, transf_d(:,:,nran,j), prop_d(:,:,j), normalisation, nopiy, nopix,plan)
			endif

			if (output_probe_intensity) then
				k = (i-1)*n_slices+j
				if (output_cell_list(k)) then
					psi = psi_d
					probe_intensity(:,:,cell_map(k)) = probe_intensity(:,:,cell_map(k)) + abs(psi)**2
				endif
			endif
        enddo ! End loop j over slices
		!If this thickness corresponds to any of the output values then accumulate diffraction pattern
		if (any(i==ncells)) then
			
			call cufftExec(plan,psi_d,psi_out_d,CUFFT_FORWARD)
			call cuda_mod<<<blocks,threads>>>(psi_out_d,temp_d,normalisation,nopiy,nopix)
			z_indx = minloc(abs(ncells-i))
			
			! Accumulate elastic wave function
			call cuda_addition<<<blocks,threads>>>(psi_elastic_d(:,:,z_indx(1)),psi_d,psi_elastic_d(:,:,z_indx(1)),1.0_fp_kind,nopiy,nopix)

			! Accumulate diffaction pattern
			!call cufftExec(plan,psi_d,psi_out_d,CUFFT_FORWARD)
			call cuda_addition<<<blocks,threads>>>(cbed_d(:,:,z_indx(1)),temp_d,cbed_d(:,:,z_indx(1)),1.0_fp_kind,nopiy,nopix)
			
			if(ionization) then
			do ii=1,num_ionizations
				stem_ion_image(ny,nx,i_df,z_indx(1),ii) = stem_ion_image(ny,nx,i_df,z_indx(1),ii)+get_sum(ion_image_d(:,:,ii))
			enddo ! ii
			endif

			! Accumulate EFISTEM images
			if (istem.and.i_df==1) then;do l = 1, imaging_ndf
				call cuda_image(psi_out_d,ctf_d(:,:,l),temp_d,normalisation, nopiy, nopix,plan,.false.)
				call cuda_addition<<<blocks,threads>>>(istem_image_d(:,:,l,z_indx(1)), temp_d, istem_image_d(:,:,l,z_indx(1)), 1.0_fp_kind, nopiy, nopix)
			enddo;endif;
		endif
		enddo ! End loop i over cells
		
		! Plasmon scattering counting
		if (plasmonmc) call pl_populate(pl_exc_num)
			
		enddo ! End loop i_qep_pass over QEP passes
        
        intensity = get_sum(psi_d)
        cbed= cbed_d
		if (output_thermal.and.fourDSTEM) psi_elastic=psi_elastic_d
		do iz=1,nz

			! Integrate the diffraction pattern
			do idet = 1, ndet
				stem_image(ny,nx,i_df,idet,iz) = cuda_stem_detector(cbed_d(:,:,iz),masks_d(:,:,idet))
			enddo
			if (ionization.and.(.not.EDX)) eels_correction_image(ny,nx,i_df,iz) = cuda_stem_detector(cbed_d(:,:,iz),eels_correction_detector_d)
        
			! Integrate the elastic diffraction pattern
			call cufftExec(plan,psi_elastic_d(:,:,iz),psi_out_d,CUFFT_FORWARD)
			call cuda_mod<<<blocks,threads>>>(psi_out_d,temp_d,normalisation,nopiy,nopix) 
			do idet = 1, ndet
				stem_elastic_image(ny,nx,i_df,idet,iz)=cuda_stem_detector(temp_d,masks_d(:,:,idet))
			enddo
			psi_out = psi_out_d*sqrt(normalisation)
			
#else
		psi_elastic=0_fp_kind
			do i_qep_pass = 1, n_qep_passes 
			if (ionization) ion_image=0.0_fp_kind
            ! Reset wavefunction
            psi = psi_initial
            
			! Plasmon scattering reset, JB-191213
			if (plasmonmc) call pl_reset()

            do i = 1,maxval(ncells)
	            do j = 1, n_slices
                    ! Accumulate ionization cross section
					if(ionization) then
                        do ii=1,num_ionizations
					        temp = abs(psi)**2 * ionization_potential(:,:,ii,j) * prop_distance(j)
					        ion_image(:,:,ii) = temp+ion_image(:,:,ii)
                        enddo
					endif
					
					! plasmon, inserted by JB-191213
					if (plasmonmc) then ! use plasmon code
						pl_tslc = prop_distance(j) ! store current slice thickness
						pl_grid_sqx = 1.0_fp_kind / a0_slice(1,j) / ak1 ! store reciprocal grid sampling along x in rad/pixel
						pl_grid_sqy = 1.0_fp_kind / a0_slice(2,j) / ak1 ! store reciprocal grid sampling along y in rad/pixel
						pl_idum = idum ! let plasmon module take over the random seed
					end if
                    
					call qep_multislice_iteration(psi,prop(:,:,j),qep_grates(:,:,:,j),nopiy,nopix,ifactory,ifactorx,idum,n_qep_grates,qep_mode,shift_arrayy,shift_arrayx)
                    
					if (output_probe_intensity) then
						k = (i-1)*n_slices+j
						if (output_cell_list(k)) then
							probe_intensity(:,:,cell_map(k)) = probe_intensity(:,:,cell_map(k)) + abs(psi)**2
						endif
					endif
		        enddo ! End loop over slices
				
                !If this thickness corresponds to any of the output values then accumulate diffration pattern
				if (any(i==ncells)) then
					z_indx = minloc(abs(ncells-i))
                    
                    ! Accumulate elastic wave function - this will be Fourier transformed later
					psi_elastic(:,:,z_indx(1)) = psi_elastic(:,:,z_indx(1)) + psi
                    
					!Transform into diffraction space
					call fft2(nopiy,nopix,psi,nopiy,psi_out,nopiy)
                    
					! Accumulate diffaction pattern
					temp = abs(psi_out)**2
					cbed(:,:,z_indx(1)) = cbed(:,:,z_indx(1)) + temp

					if(ionization) stem_ion_image(ny,nx,i_df,z_indx(1),:) = stem_ion_image(ny,nx,i_df,z_indx(1),:)+ sum(sum(ion_image,dim=2),dim=1)
                endif
                
            enddo ! End loop over cells
			
			! Plasmon scattering counting
			if (plasmonmc) call pl_populate(pl_exc_num)
			
		enddo !End loop over QEP passes
		intensity = sum(abs(psi)**2)
        
		do iz=1,nz
			
        
			! Integrate the elastic diffraction pattern
			call fft2(nopiy,nopix,psi_elastic(:,:,iz),nopiy,psi_out,nopiy)
            if(stem) stem_image(ny,nx,i_df,1:ndet,iz) = sum(sum(spread(cbed(:,:,iz),dim=3,ncopies=ndet)*masks(:,:,:),dim=1),dim=1)
            if(stem) stem_elastic_image(ny,nx,i_df,1:ndet,iz)= sum(sum(spread(abs(psi_out)**2,dim=3,ncopies =ndet)*masks(:,:,1:ndet),dim=1),dim=1)
			if(ionization.and.(.not.EDX)) eels_correction_image(ny,nx,i_df,z_indx(1)) = sum(cbed(:,:,iz)*eels_correction_detector)
#endif            
            if(pacbed) then
			    !Output 4D STEM diffraction pattern
			    if(fourDSTEM) then
					    !Output total (elastic and inelastic) diffraction pattern
					    filename = trim(adjustl(output_prefix))
						if (probe_ndf>1) filename = trim(adjustl(filename))//defocus_string(probe_df(i_df),lengthdf)
					    if (nz>1) filename = trim(adjustl(filename))//'_z='//to_string(int(zarray(iz)))//'_A'
					    call binary_out_unwrap(nopiy, nopix, cbed(:,:,iz)/n_qep_passes, trim(adjustl(filename)) //'_pp_'//&
										      &to_string(nx)//'_'//to_string(ny)//'_Diffraction_pattern',write_to_screen=.false.&
                                              &,nopiyout=nopiyout,nopixout=nopixout)
			    endif

			    if(i_df==1) pacbed_pattern(:,:,iz) = pacbed_pattern(:,:,iz) + cbed(:,:,iz)/n_qep_passes
                if(output_thermal) pacbed_elastic(:,:,iz) =pacbed_elastic(:,:,iz) + abs(psi_out)**2/n_qep_passes**2

			    if (output_thermal.and.fourDSTEM) then 
					    !Output elastic only diffraction pattern
                        filename = trim(adjustl(output_prefix))
						if (probe_ndf>1) filename = trim(adjustl(filename))//defocus_string(probe_df(i_df),lengthdf)
					    if (nz>1) filename = trim(adjustl(filename))//'_z='//to_string(int(zarray(iz)))//'_A'
					    filename = trim(adjustl(filename))//'_pp_'//to_string(nx)//'_'//to_string(ny)//'_Elastic_Diffraction_pattern'
					    call binary_out_unwrap(nopiy, nopix, abs(psi_out)**2/n_qep_passes**2, filename,write_to_screen=.false.&
                                              &,nopiyout=nopiyout,nopixout=nopixout)
                endif
            endif
		enddo
			

                
		if (output_probe_intensity) call probe_intensity_to_file(probe_intensity,i_df,ny,nx,n_qep_passes,probe_ndf,nysample,nxsample)
	enddo ! End loop over x probe positions
	enddo ! End loop over y probe positions
	enddo ! End loop over defocus series
	
    ! QEP normalisation
    
    if(ionization) then
		stem_ion_image = stem_ion_image/float(n_qep_passes)
		if(.not.EDX) eels_correction_image = eels_correction_image/float(n_qep_passes)
	endif

    if (ndet.gt.0) then
	    stem_image = stem_image/float(n_qep_passes)
        stem_elastic_image = stem_elastic_image/float(n_qep_passes*n_qep_passes)
    endif

    delta = secnds(t1)
    
    write(*,*) 
    write(*,*) 
    
    write(*,*) 'Calculation is finished.'
    write(*,*) 
    write(*,*) 'Time elapsed ', delta, ' seconds.'
    write(*,*)    
    
	if(timing) then
		open(unit=9834, file=trim(adjustl(output_prefix))//'_timing.txt', access='append')
		write(9834, '(a, g, a, /)') 'The multislice calculation took ', delta, 'seconds.'
		close(9834)
    endif

    if (fp_kind.eq.8) then
        write(*,*) 'The following files were outputted (as 64-bit big-endian floating point):'
	else
        write(*,*) 'The following files were outputted (as 32-bit big-endian floating point):'
	endif
    
    write(*,*)
    fnam = trim(adjustl(output_prefix))
    if (n_tilts_total>1) fnam = trim(adjustl(output_prefix))//tilt_description(claue(:,ntilt),ak1,ss,ig1,ig2)
    stem_inelastic_image = stem_image - stem_elastic_image
    do idet = 1, ndet
        
        if(output_thermal) then
			fnam_temp = trim(adjustl(fnam)) // '_DiffPlaneElastic_Detector'
			call add_zero_padded_int(fnam_temp, fnam_det, idet, 2)
			call output_stem_image(stem_elastic_image(:,:,:,idet,:),fnam_det,probe_df)
        
			
			fnam_temp = trim(adjustl(fnam)) // '_DiffPlaneTDS_Detector'
			call add_zero_padded_int(fnam_temp, fnam_det, idet, 2)
			call output_stem_image(stem_inelastic_image(:,:,:,idet,:),fnam_det,probe_df)
        
			fnam_temp = trim(adjustl(fnam)) // '_DiffPlaneTotal_Detector'
		else
			fnam_temp = trim(adjustl(fnam)) // '_DiffPlane_Detector'
		endif

        call add_zero_padded_int(fnam_temp, fnam_det, idet, 2)
        call output_stem_image(stem_image(:,:,:,idet,:),fnam_det,probe_df)
    enddo

    !ionization
    if(ionization) then
        do ii=1,num_ionizations
        if(EDX) then
            filename = trim(adjustl(fnam)) // '_'//trim(adjustl(substance_atom_types(atm_indices(ii))))//'_'//trim(adjustl(Ion_description(ii)))//'_shell_EDX'
            call output_stem_image(stem_ion_image(:,:,:,:,ii), filename,probe_df)
        else
            filename = trim(adjustl(fnam))// '_'//trim(adjustl(substance_atom_types(atm_indices(ii))))//'_'//trim(adjustl(Ion_description(ii)))//'_orbital_EELS'
            call output_stem_image(stem_ion_image(:,:,:,:,ii), filename,probe_df)
            stem_ion_image(:,:,:,:,ii) = stem_ion_image(:,:,:,:,ii)*eels_correction_image

            filename =  trim(adjustl(fnam)) //'_'//trim(adjustl(substance_atom_types(atm_indices(ii))))//'_'//trim(adjustl(Ion_description(ii)))//'_orbital_EELS_Corrected'            
            call output_stem_image(stem_ion_image(:,:,:,:,ii), filename,probe_df)
        endif
        enddo
        if(.not.EDX) then
            filename = trim(adjustl(fnam)) // '_EELS_CorrectionMap' 
            call output_stem_image(eels_correction_image, filename,probe_df)
        endif
       
    endif
    
    if(pacbed) then
            const = float(nysample*nxsample)
        	do i=1,nz
            
		    filename = trim(adjustl(output_prefix))
            if (n_tilts_total>1) fnam = trim(adjustl(output_prefix))//tilt_description(claue(:,ntilt),ak1,ss,ig1,ig2)
		    if(nz>1) filename=trim(adjustl(filename))//'_z='//zero_padded_int(int(zarray(i)),length)//'_A_'
		    call binary_out_unwrap(nopiy,nopix,pacbed_pattern(:,:,i)/const,trim(adjustl(filename))//'_PACBED_Pattern',nopiyout=nopiyout,nopixout=nopixout)
        
            if(output_thermal) then
                call binary_out_unwrap(nopiy,nopix,PACBED_elastic(:,:,i)/const,trim(adjustl(filename))//'_elastic_PACBED_Pattern',nopiyout=nopiyout,nopixout=nopixout)    
                call binary_out_unwrap(nopiy,nopix,pacbed_pattern(:,:,i)/const-PACBED_elastic(:,:,i)/const,trim(adjustl(filename))//'_thermal_PACBED_Pattern',nopiyout=nopiyout,nopixout=nopixout)
            endif
	    enddo
    endif

#ifdef GPU
				if(double_channeling.and.istem) then
				do i=1,nz;do l=1,imaging_ndf
				temp = efistem_image_d(:,:,l,i)
				call output_TEM_result(output_prefix,tile_out_image(temp,ifactory,ifactorx)/n_qep_passes/nysample/nxsample,'energy_filtered_ISTEM',nopiy,nopix,manyz,imaging_ndf>1,manytilt,z=zarray(i)&
									&,lengthz=length,lengthdf=lengthimdf,tiltstring = tilt_description(claue(:,ntilt),ak1,ss,ig1,ig2),df = imaging_df(l))
			enddo;enddo;endif	
		
	if(double_channeling) then
		do l=1,numeels
			filename =  trim(adjustl(fnam))//'_double_channeling_EELS_'//zero_padded_int(l,2)
			call output_stem_image(Hn0_eels_dc(:,:,:,:,l), filename,probe_df)
		enddo
		endif

		if(istem) then
				do i=1,nz;do l=1,imaging_ndf
				temp = istem_image_d(:,:,l,i)
				call output_TEM_result(output_prefix,tile_out_image(temp,ifactory,ifactorx)/n_qep_passes/nysample/nxsample,'ISTEM',nopiy,nopix,manyz,imaging_ndf>1,manytilt,z=zarray(i)&
									&,lengthz=length,lengthdf=lengthimdf,tiltstring = tilt_description(claue(:,ntilt),ak1,ss,ig1,ig2),df = imaging_df(l))
			enddo;enddo;endif	
#endif
	enddo
	
end subroutine qep_stem
