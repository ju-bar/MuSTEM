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
  use m_conv2dr ! convolution module, used for applying transmission point spread for SE imaging
    
  implicit none
    
  logical,intent(in)::STEM,ionization,PACBED
  !dummy variables
  integer(4) :: i,j,l,m,i_qep_pass,iz,k,ii,jj,ntilt,cntscan,numscan
  integer(4) :: shifty,shiftx,ny,nx,i_df,idet,total_slices
  integer(4) :: lengthimdf,starting_slice

  !random variables
  integer(4) :: idum,nsliceoutput,z_indx(1)
  !real(fp_kind) :: ran1
  complex(fp_kind) :: prop(nopiy,nopix,n_slices)!,temp_transf(nopiy,nopix)
       
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
  real(fp_kind),allocatable :: ion_image(:,:,:),stem_ion_image(:,:,:,:,:),pacbed_pattern(:,:,:),&
                        &pacbed_elastic(:,:,:),istem_image(:,:,:,:)
    
  !diagnostic variables
  real(fp_kind) :: intensity, intens
  ! intens added May 2019 to accommodate changes awhere it occurs below 
  real(fp_kind) :: t1, dnow, dnext, delta, dest, dest_last
    
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
  real(fp_kind),device,allocatable :: se_transmission_d(:,:,:) ! per slice secondary electron transmission function, JB-2025-06-26
  real(fp_kind),device,allocatable :: se_acctransm_d(:,:) ! accumulated SE transmission function, JB-2025-06-26
    
  !device variables for on the fly potentials
  complex(fp_kind),device,allocatable,dimension(:,:) :: bwl_mat_d,inverse_sinc_d,fz_dwf_d
  complex(fp_kind),device,allocatable,dimension(:,:,:) :: fz_d,fz_mu_d
  real(fp_kind),device,allocatable :: inelastic_potential_d(:,:)

  !Double channeling variables
  integer,allocatable :: inel_wf_valid(:,:) ! validity flags for inelastic wave functions (1: valid, 0: invalid), result of strength filter, (target atom, state)
  real(fp_kind),allocatable :: inel_wf_intens(:,:) ! intensity of inelastic wave functions, input for strength filter, (target atom, state)
  complex(fp_kind),device,allocatable,dimension(:,:) ::psi_inel_d,shiftarray,tmatrix_d,q_tmatrix_d
  complex(fp_kind),device,allocatable,dimension(:,:,:)::tmatrix_states_d,ctf_d
  complex(fp_kind),device,allocatable,dimension(:,:,:)::Hn0_shifty_coord_d,Hn0_shiftx_coord_d
  real(fp_kind),device,allocatable,dimension(:,:)::cbed_inel_dc_d
  real(fp_kind),device,allocatable,dimension(:,:,:)::Hn0_eels_detector_d
  real(fp_kind),device,allocatable,dimension(:,:,:,:)::efistem_image_d,istem_image_d
  real(fp_kind),allocatable,dimension(:,:,:,:,:)::Hn0_eels_dc
  real(fp_kind),allocatable,dimension(:,:,:)::tmatrix_states
  !real(fp_kind),allocatable,dimension(:,:,:,:,:) :: inel_intens, inel_intens_valid ! 2025-07-28 JB: (debug) inelastic intensities per atom, transition, z, y, x
  integer::i_target
#endif

  real(fp_kind),allocatable :: probe_intensity(:,:,:)
  character(1024) :: filename
  real(fp_kind)::qep_stem_GPU_memory,const,rnorm
  logical::elfourd,manytilt,manyz,many_df,has_eta
  integer*4::length,lengthdf
  real(fp_kind),allocatable :: se_acctransm(:,:) ! accumulated SE transmission function, JB-2025-06-26
  
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
  if(tp_eels) then
    allocate(tmatrix_states_d(nopiy,nopix,nstates),psi_inel_d(nopiy,nopix),cbed_inel_dc_d(nopiy,nopix),tmatrix_states(nopiy,nopix,nstates))
    allocate(shiftarray(nopiy,nopix),tmatrix_d(nopiy,nopix),q_tmatrix_d(nopiy,nopix))
    tmatrix_states_d = setup_ms_hn0_tmatrices(nopiy,nopix,nstates)*alpha_n
    allocate(Hn0_shifty_coord_d(nopiy,maxval(natoms_slice_total),n_slices))
    allocate(Hn0_shiftx_coord_d(nopix,maxval(natoms_slice_total),n_slices))
    Hn0_shiftx_coord_d = Hn0_shiftx_coord
    Hn0_shifty_coord_d = Hn0_shifty_coord
    !allocate(inel_intens(maxval(natoms_slice_total), nstates, n_slices*maxval(ncells), nysample, nxsample))
    !allocate(inel_intens_valid(maxval(natoms_slice_total), nstates, n_slices*maxval(ncells), nysample, nxsample))
    !inel_intens = 0.0_fp_kind
    !inel_intens_valid = 0.0_fp_kind
    allocate(inel_wf_valid(maxval(natoms_slice_total), nstates)) ! validity flags for inelastic wave functions, 2025-07-30 JB
    allocate(inel_wf_intens(maxval(natoms_slice_total), nstates)) ! intensity of inelastic wave functions, 2025-07-30 JB
    allocate(Hn0_eels_dc(nysample,nxsample,probe_ndf,nz,dc_numeels))
    allocate(Hn0_eels_detector_d(nopiy,nopix,dc_numeels),Hn0_eels_detector(nopiy,nopix,dc_numeels))
    !do l=1,numeels ! 2025-07-14, JB: replaced by annular and segmented EELS options in DC, see below
    !  Hn0_eels_detector(:,:,l) = make_detector(nopiy,nopix,ifactory,ifactorx,ss,0.0_fp_kind,outerrad(l))
    !enddo
    do i=1,dc_numeels/dc_eels_nseg; do j=1,dc_eels_nseg
        l=(i-1)*dc_eels_nseg+j
		if(dc_eels_nseg>1) Hn0_eels_detector(:,:,l) = make_detector(nopiy,nopix,ifactory,ifactorx, &
                &   ss,dc_eels_inner(i),dc_eels_outer(i),2*pi*j/dc_eels_nseg-dc_eels_seg_offset,2*pi/dc_eels_nseg)
		if(dc_eels_nseg==1) Hn0_eels_detector(:,:,l) = make_detector(nopiy,nopix,ifactory,ifactorx, &
                &   ss,dc_eels_inner(i),dc_eels_outer(i))
	enddo; enddo
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
      ctf(:,:,i) =  make_ctf([0.0_fp_kind,0.0_fp_kind,0.0_fp_kind],imaging_df(i),imaging_cutoff,&
                          &imaging_aberrations,imaging_apodisation)
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
    
! Pre-calculate the scattering factors on a grid
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
    if(.not. load_grates) then
      call load_save_add_grates(idum,projected_potential,nopiy,nopix,n_qep_grates,n_slices,nt,nat_slice)
    endif
    call make_local_inelastic_potentials(ionization)  !setup the REAL space inelastic potential (ionization and adf) for QUEP ADF is disabled
  endif

  if (SEI) then ! 2025-06-26, JB: Secondary electron transmission function (always not on-the-fly)
    call make_se_transmission_grates() ! calculate the transmission function for SE per slice
    allocate(se_acctransm(nopiy,nopix)) ! allocate the accumulated SE transmission function
    se_acctransm = 1.0_fp_kind ! initialize to 1 (full transmission)
    ! setup convolution module for applying the PSF to the SE transmission functions
    call conv2dr_init(n_se_psf) ! initialize with number of PSFs
    do i=1, n_se_psf
      call conv2dr_set_kernel(i, se_psf(:,:,i), nopiy, nopix) ! set the kernel data
    enddo
  endif

  if(stem) then
    do i=1,ndet/nseg
      do j=1,nseg
        if(nseg>1) then
          masks(:,:,(i-1)*nseg+j) = make_detector(nopiy,nopix,ifactory,ifactorx,ss,inner(i),outer(i),&
                            &2*pi*j/nseg-seg_det_offset,2*pi/nseg)
        endif
        if(nseg==1) then
          masks(:,:,(i-1)*nseg+j) = make_detector(nopiy,nopix,ifactory,ifactorx,ss,inner(i),outer(i))
        endif
      enddo
    enddo
  endif
  

  t1 = secnds(0.0)
  rnorm = 1.0_fp_kind
  intensity = 1.0d0
  numscan = nysample * nxsample * probe_ndf * n_tilts_total
  cntscan = 0
  dnow = 0.0_fp_kind
  dnext = 15.0_fp_kind
  has_eta = .false.
  dest = 0.0_fp_kind
  dest_last = 0.0_fp_kind
  
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
    if(EELS) then ! 2025-06-07, JB: EELS flag replaces previous .not.EDX conditions
      allocate(eels_correction_detector_d(nopiy,nopix))
      eels_correction_detector_d=eels_correction_detector
    endif
    if(SEI) then ! 2025-06-26, JB: SE transmission functions on device
      allocate(se_transmission_d(nopiy,nopix,n_slices),se_acctransm_d(nopiy,nopix))
      se_transmission_d = se_transmission ! copy from host to device
      se_acctransm_d = se_acctransm ! initialize the accumulated SE transmission function
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
    if (EELS) eels_correction_image = 0.0_fp_kind ! 2025-06-07, JB: EELS flag replaces previous .not.EDX conditions
    !call binary_out_unwrap(nopiy,nopix,eels_correction_detector,'eels_correction_detector')
  endif
  if (output_probe_intensity) allocate(probe_intensity(nopiy,nopix,size(output_thickness_list)))
  
  
  ! loop over the specimen tilts
  do ntilt=1,n_tilts_total
    do i = 1, n_slices
      call make_propagator(nopiy,nopix,prop(:,:,i),prop_distance(i),Kz(1),ss,ig1,ig2,claue(:,1),ifactorx,ifactory)
      prop(:,:,i) = prop(:,:,i) * bwl_mat
      qep_grates(:,:,:,i) = exp(ci*pi*a0_slice(3,i)/Kz(ntilt)*projected_potential(:,:,:,i))
      do j=1,n_qep_grates
        call fft2(nopiy,nopix,qep_grates(:,:,j,i),psi)
        if(qep_mode.eq.3) qep_grates(:,:,j,i)= psi*bwl_mat
        if(qep_mode.ne.3) call ifft2(nopiy,nopix,psi*bwl_mat,qep_grates(:,:,j,i))
      enddo ! n_qep_grates
    enddo ! n_slices
#ifdef GPU
    prop_d=prop
    if(.not.on_the_fly) transf_d = qep_grates
#endif        
    lengthdf = ceiling(log10(maxval(abs(probe_df))))
    if(any(probe_df<0)) lengthdf = lengthdf+1

    
! --------------------------------------------------------------------
! >> SCANNING LOOPS (code formatting relaxed due to too many loops) <<
    
    do i_df = 1, probe_ndf
    do ny = 1, nysample
    do nx = 1, nxsample
        
    if (cntscan > 0 .and. cntscan < numscan-1) then ! valid progress for eta
        if (dnow > dnext) then
            dest = dnow / cntscan * (real(numscan-cntscan,fp_kind))
            dnext = dnext + 3.0_fp_kind
            dest_last = dest
            has_eta = .true.
        else
            dest = dest_last
        end if
    end if
        
#ifdef GPU
    if (has_eta .AND. cntscan < numscan-1) then
        write(6,902,advance='no') achar(13), i_df, probe_ndf, ny, nysample, &
                    & nx, nxsample, intensity, int(dest)
    else
        write(6,901,advance='no') achar(13), i_df, probe_ndf, ny, nysample, &
                    & nx, nxsample, intensity
    end if
901 format(a1,' df:',i3,'/',i3,' y:',i3,'/',i3,' x:',i3,'/',i3, &
                    & '  Intensity:',f6.3,' (to monitor BWL)      ')
902 format(a1,' df:',i3,'/',i3,' y:',i3,'/',i3,' x:',i3,'/',i3, &
                    & '  Intensity:',f6.3,' (ETA: ',i0,' s)       ')
#else
    if (has_eta .AND. cntscan < numscan-1) then
        write(6,902) i_df, probe_ndf, ny, nysample, nx, nxsample, intensity, int(dest)
    else
        write(6,901) i_df, probe_ndf, ny, nysample, nx, nxsample, intensity
    end if
901 format(1h+,1x,' df:',i3,'/',i3,' y:',i3,'/',i3,' x:',i3,'/',i3, &
                    &'  Intensity:',f6.3,' (to monitor BWL)      ')
902 format(1h+,1x,' df:',i3,'/',i3,' y:',i3,'/',i3,' x:',i3,'/',i3, &
                    &'  Intensity:',f6.3,' (ETA: ',i0,' s)       ')
#endif      
    flush(6)
!
!       Make STEM probe
!
!		psi_initial = make_ctf(probe_positions(:,ny,nx),probe_df(i_df),probe_cutoff,probe_aberrations,probe_apodisation)
!       call ifft2(nopiy, nopix, psi_initial, nopiy, psi_initial, nopiy)
!       psi_initial = psi_initial/sqrt(sum(abs(psi_initial)**2))
! The above three lines replaced by the following four to improve accuracy (May 2019) - Findlay suggestion
    psi_initial = make_ctf(probe_positions(:,ny,nx),probe_df(i_df),probe_cutoff,probe_aberrations,probe_apodisation)
    intens = sum(abs(psi_initial)**2)
    call ifft2(nopiy, nopix, psi_initial, psi_initial)
    psi_initial = psi_initial/sqrt(intens)
    
    if (arg_debug_wave>0) write(*,*) 'Initial psi(0,0) = ', psi_initial(1,1)
    
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
      if (SEI) se_acctransm_d = se_acctransm ! reset attenuation for secondary electrons, 2025-06-26, JB
                                    ! note: se_acctransm remains as initialized, i.e. 1.0_fp_kind for GPU code
      do i = 1,maxval(ncells) ! main multislice loop of cells ! code formatting relaxed due to too many loops
      do j = 1, n_slices ! sub loop over slices per cell
        ! Accumulate ionization cross section
        if(ionization) then
          do ii=1,num_ionizations
            call cuda_mod<<<blocks,threads>>>(psi_d,temp_d,1.0_fp_kind,nopiy,nopix)
            ! overlap
            if(on_the_fly) then
              call cuda_make_ion_potential(inelastic_potential_d,tau_slice(:,atm_indices(ii),:,j),&
                            &nat_slice(atm_indices(ii),j),plan,&
                            &fz_mu_d(:,:,ii),inverse_sinc_d,Volume_array(j))
              call cuda_multiplication<<<blocks,threads>>>(temp_d,inelastic_potential_d,&
                            &temp_d,prop_distance(j),nopiy,nopix)
            else
              call cuda_multiplication<<<blocks,threads>>>(temp_d,ion_potential_d(:,:,ii,j),&
                            &temp_d,prop_distance(j),nopiy,nopix)
            endif
            if (SEI) then ! 2025-06-26, JB: Apply SE transmission function
              ! Multiply the current accumulated SE transmission function to
              ! the ionization cross section overlap of the current layer.
              call cuda_multiplication<<<blocks,threads>>>(temp_d,se_acctransm_d,&
                            &temp_d,1.0_fp_kind,nopiy,nopix)
            endif
            
            ! depth sum, 2025-06-26, JB
            call cuda_addition<<<blocks,threads>>>(ion_image_d(:,:,ii),temp_d,ion_image_d(:,:,ii),&
                            &1.0_fp_kind,nopiy,nopix)  
          enddo ! ii
        endif

        !Double channeling, i.e., propagation of inelastic wave functions through the specimen
        if(tp_eels) then
            
          inel_wf_intens = 0.0_fp_kind ! init all iwfs intensities with zero, 2025-07-30 JB
          inel_wf_valid = 1 ! init all iwfs as valid
          if (.NOT.single_channeling .AND. iwf_filter_threshold > 1.0e-6_fp_kind) then ! prepare data for inel-wf filter only when real double channeling
          ! This only runs when double-channeling is enabled, i.e., single_channeling is false,
          ! because the overhead of calculating the intensities of inelastic wavefunctions is
          ! equal to the cost of calculating single channeling results. This means we would not save
          ! any computation but actually do more work.
          ! We may want to skip this section also when exciting close to the exit plane, because
          ! then double channeling doesn't cost a lot of computation and the filter is usually also
          ! less effective.
          do i_target = 1, natoms_slice_total(j) ! Loop over targets
            ! calculate shift matrix for current target atom
            call cuda_make_shift_array<<<blocks,threads>>>(shiftarray,Hn0_shifty_coord_d(:,i_target,j),&
                            &Hn0_shiftx_coord_d(:,i_target,j),nopiy,nopix)
            do k = 1, nstates
              if (state_use(k) == 0) cycle ! skip unused states, Hn0 strength filter, 2025-07-28 JB
              ! shift current Hn0 to atomic position
              call cuda_multiplication<<<blocks,threads>>>(tmatrix_states_d(:,:,k),shiftarray, &
                            &q_tmatrix_d,1.0_fp_kind,nopiy,nopix)
              ! transform to real space
              call cufftExec(plan,q_tmatrix_d,tmatrix_d,CUFFT_INVERSE)
              ! multiply with the wave function psi_d * tmatrix_d --> psi_inel_d
              call cuda_multiplication<<<blocks,threads>>>(psi_d,tmatrix_d,psi_inel_d,&
                            &sqrt(normalisation),nopiy,nopix)
              ! store the inelastic wave function intensity
              inel_wf_intens(i_target,k) = get_sum(psi_inel_d) ! 2025-07-30 JB: store the intensity of the inelastic wave function
              ! (debug)calculate the total intensity of the inelastic wave function
              !ii = j + (i-1)*n_slices ! i.e., the current slice index in the whole specimen
              !inel_intens(i_target,k,ii,ny,nx) = get_sum(psi_inel_d)
              !inel_intens_valid(i_target,k,ii,ny,nx) = 1.0_fp_kind
            enddo
          enddo
          ! inelastic wave function strength filter, 2025-07-30 JB
          call strength_filter(inel_wf_intens, inel_wf_valid, maxval(natoms_slice_total)*nstates, &
                            &  iwf_filter_threshold)
          endif
          
          do i_target = 1, natoms_slice_total(j) ! Loop over targets
            ! calculate shift matrix for current target atom
            call cuda_make_shift_array<<<blocks,threads>>>(shiftarray,Hn0_shifty_coord_d(:,i_target,j),&
                            &Hn0_shiftx_coord_d(:,i_target,j),nopiy,nopix)
            do k = 1, nstates
              if (state_use(k) == 0) cycle ! skip unused states, Hn0 strength filter, 2025-07-28 JB
              if (inel_wf_valid(i_target,k) == 0) cycle ! skip invalidated (too weak) inelastic wave functions, 2025-07-30 JB
              ! shift current Hn0 to atomic position
              call cuda_multiplication<<<blocks,threads>>>(tmatrix_states_d(:,:,k),shiftarray, &
                            &q_tmatrix_d,1.0_fp_kind,nopiy,nopix)
              ! transform to real space
              call cufftExec(plan,q_tmatrix_d,tmatrix_d,CUFFT_INVERSE)
              ! multiply with the wave function psi_d * tmatrix_d --> psi_inel_d
              call cuda_multiplication<<<blocks,threads>>>(psi_d,tmatrix_d,psi_inel_d,&
                            &sqrt(normalisation),nopiy,nopix)
              
              starting_slice = j
              
              do ii = i, n_cells ! Scatter the inelastic wave through the remaining cells
                ! if the following loop is skipped, then we would calculate single channeling
                !
                if (.NOT.single_channeling) then ! 2025-07-28 JB: skipping actual double channeling
                do jj = starting_slice, n_slices
                  ! QEP multislice
                  nran = floor(n_qep_grates*ran1(idum)) + 1
                  shiftx = floor(ifactorx*ran1(idum));shifty = floor(ifactory*ran1(idum))
                  if(on_the_fly) then
                    call cuda_fph_make_potential(trans_d,ccd_slice_array(jj),tau_slice,nat_slice(:,jj),&
                            &jj,prop_distance(jj),idum,plan,fz_d,inverse_sinc_d,bwl_mat_d)
                    call cuda_multislice_iteration(psi_inel_d, trans_d, prop_d(:,:,jj), normalisation,&
                            &nopiy, nopix,plan)
                  elseif(qep_mode == 2) then
                    call cuda_multislice_iteration(psi_inel_d, transf_d(:,:,nran,j), prop_d(:,:,j),&
                            &normalisation, nopiy, nopix,shifty*nopiy_ucell,shiftx* nopix_ucell,plan)
                  elseif(qep_mode == 3) then                       !randomly shift phase grate
                    call cuda_phase_shift_from_1d_factor_arrays(transf_d(:,:,nran,jj),trans_d,&
                            &shift_arrayy_d(:,shifty+1),shift_arrayx_d(:,shiftx+1),nopiy,nopix,plan)
                    call cuda_multislice_iteration(psi_inel_d, trans_d, prop_d(:,:,jj), normalisation,&
                            &nopiy, nopix,plan)
                  else
                    call cuda_multislice_iteration(psi_inel_d, transf_d(:,:,nran,jj), prop_d(:,:,jj),&
                            &normalisation, nopiy, nopix,plan)
                  endif
                enddo ! jj
                endif ! skip double channeling propagation 
                starting_slice=1
                if (any(ii==ncells)) then
                  z_indx = minloc(abs(ncells-ii))
                  call cufftExec(plan, psi_inel_d, psi_out_d, CUFFT_FORWARD)
							
                  ! Accumulate the EELS images
                  do l=1,dc_numeels
                    Hn0_eels_dc(ny,nx,i_df,z_indx(1),l) = Hn0_eels_dc(ny,nx,i_df,z_indx(1),l)+&
                            &cuda_stem_detector(psi_out_d, Hn0_eels_detector_d(:,:,l))
                  enddo ! l
                  ! Accumulate EFISTEM images
                  if (istem.and.i_df==1) then
                    do l = 1, imaging_ndf
                      call cuda_image(psi_out_d,ctf_d(:,:,l),temp_d,normalisation, nopiy, nopix,plan,.false.)
                      call cuda_addition<<<blocks,threads>>>(efistem_image_d(:,:,l,z_indx(1)), temp_d,&
                            &efistem_image_d(:,:,l,z_indx(1)), 1.0_fp_kind, nopiy, nopix)
                    enddo
                  endif
							    !stop
                endif
              enddo ! ii
            enddo ! k
          enddo ! i_target
        endif  ! End loop over cells,targets and states and end tp_eels section

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
          call cuda_fph_make_potential(trans_d,ccd_slice_array(j),tau_slice,nat_slice(:,j),j,&
                            &prop_distance(j),idum,plan,fz_d,inverse_sinc_d,bwl_mat_d)
          call cuda_multislice_iteration(psi_d, trans_d, prop_d(:,:,j), normalisation, nopiy, nopix,plan)
        elseif(qep_mode == 2) then
          call cuda_multislice_iteration(psi_d, transf_d(:,:,nran,j), prop_d(:,:,j), normalisation,&
                            &nopiy, nopix,shifty*nopiy_ucell,shiftx*nopix_ucell,plan)
        elseif(qep_mode == 3) then                       !randomly shift phase grate
          call cuda_phase_shift_from_1d_factor_arrays(transf_d(:,:,nran,j),trans_d,shift_arrayy_d(:,shifty+1),&
                            &shift_arrayx_d(:,shiftx+1),nopiy,nopix,plan)
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
        
        if (SEI) then ! 2025-06-26, JB: Update accumulated SE transmission function
          ! multiply the TF of the current slice to the accumulated SE transmission function
          call cuda_multiplication<<<blocks,threads>>>(se_acctransm_d,se_transmission_d(:,:,j),&
                            &se_acctransm_d,1.0_fp_kind,nopiy,nopix)
          ! 2025-08-21, JB: added simulation of the secondary electron isotropic momentum taking
          !                 effect in the propagation to the detector and on the absorption
          ! convolute with the PSF of the current slice to model the point spread of
          ! the secondary electrons propagating with perpendicular momenta and backwards to the detector
          call conv2dr_apply(1+MODULO(j-1,n_se_psf), se_acctransm_d, se_acctransm_d, nopiy, nopix)
        endif
        
      enddo ! End loop j over slices

      !If this thickness corresponds to any of the output values then accumulate diffraction pattern
      if (any(i==ncells)) then
        z_indx = minloc(abs(ncells-i))
        
        if (arg_debug_wave>0) then
            psi = psi_d ! download psi
            write(*,*) 'psi(0,0) = ', psi(1,1), ' at z=',z_indx(1)
        endif
        
        ! Accumulate elastic wave function
        call cuda_addition<<<blocks,threads>>>(psi_elastic_d(:,:,z_indx(1)),psi_d,&
                            &psi_elastic_d(:,:,z_indx(1)),1.0_fp_kind,nopiy,nopix)
        
        ! Calculate the diffraction pattern
        call cufftExec(plan,psi_d,psi_out_d,CUFFT_FORWARD)
        call cuda_mod<<<blocks,threads>>>(psi_out_d,temp_d,normalisation,nopiy,nopix)
        
        if (arg_debug_intens>0) then
            temp = temp_d ! download psi
            write(*,*) 'cbed(0,0) = ', temp(1,1), ' at z=',z_indx(1)
        endif

        ! Accumulate diffaction pattern
        !call cufftExec(plan,psi_d,psi_out_d,CUFFT_FORWARD)
        call cuda_addition<<<blocks,threads>>>(cbed_d(:,:,z_indx(1)),temp_d,cbed_d(:,:,z_indx(1)),&
                            &1.0_fp_kind,nopiy,nopix)
			
        if(ionization) then
          do ii=1,num_ionizations
            temp = ion_image_d(:,:,ii) ! download ion_image from device to host
            stem_ion_image(ny,nx,i_df,z_indx(1),ii) = stem_ion_image(ny,nx,i_df,z_indx(1),ii) +&
                            &get_sum(ion_image_d(:,:,ii))
          enddo ! ii
        endif

        ! Accumulate EFISTEM images
        if (istem.and.i_df==1) then
          do l = 1, imaging_ndf
            call cuda_image(psi_out_d,ctf_d(:,:,l),temp_d,normalisation, nopiy, nopix,plan,.false.)
            call cuda_addition<<<blocks,threads>>>(istem_image_d(:,:,l,z_indx(1)), temp_d, &
                            &istem_image_d(:,:,l,z_indx(1)), 1.0_fp_kind, nopiy, nopix)
          enddo
        endif
      endif
      enddo ! End loop i over cells ! relaxed indentation
		
      ! Plasmon scattering counting
      if (plasmonmc) call pl_populate(pl_exc_num)
			
    enddo ! End loop i_qep_pass over QEP passes ! relaxed indentation
        
    intensity = get_sum(psi_d)
    cbed= cbed_d
    if (output_thermal.and.fourDSTEM) psi_elastic=psi_elastic_d
    do iz=1,nz

      ! Integrate the diffraction pattern
      do idet = 1, ndet
        stem_image(ny,nx,i_df,idet,iz) = cuda_stem_detector(cbed_d(:,:,iz),masks_d(:,:,idet))
      enddo
      if (arg_debug_stemdet>0) then
          cbed(:,:,iz) = cbed_d(:,:,iz) ! download cbeds
          write(*,*) 'cbed(0,0) = ', cbed(1,1,iz), ' at iz=',iz
          write(*,*) 'stem_image_1 = ', stem_image(ny,nx,i_df,1,iz), ' at iz=',iz
      end if
      if (EELS) eels_correction_image(ny,nx,i_df,iz) = & ! 2025-06-07, JB: EELS flag replaces previous .not.EDX conditions
                   & cuda_stem_detector(cbed_d(:,:,iz),eels_correction_detector_d)
        
      ! Integrate the elastic diffraction pattern
      call cufftExec(plan,psi_elastic_d(:,:,iz),psi_out_d,CUFFT_FORWARD)
      call cuda_mod<<<blocks,threads>>>(psi_out_d,temp_d,normalisation,nopiy,nopix) 
      do idet = 1, ndet
        stem_elastic_image(ny,nx,i_df,idet,iz)=cuda_stem_detector(temp_d,masks_d(:,:,idet))
      enddo
      psi_out = psi_out_d*sqrt(normalisation)
			
#else ! code indentation discontinued -> CPU code
    
    psi_elastic=0_fp_kind
    do i_qep_pass = 1, n_qep_passes 
      if (ionization) ion_image=0.0_fp_kind
      ! Reset wavefunction
      psi = psi_initial
            
      ! Plasmon scattering reset, JB-191213
      if (plasmonmc) call pl_reset()
    
      if (SEI) se_acctransm = 1.0_fp_kind ! reset accumulated SE transmission function, 2025-06-26, JB
      do i = 1,maxval(ncells)
        do j = 1, n_slices
          ! Accumulate ionization cross section
          if(ionization) then
            do ii=1,num_ionizations
              temp = abs(psi)**2 * ionization_potential(:,:,ii,j) * prop_distance(j)
              if (SEI) temp = temp * se_acctransm ! 2025-06-26, JB: Apply SE transmission function
              ion_image(:,:,ii) = temp + ion_image(:,:,ii) ! depth sum
            enddo
          endif
					
          ! plasmon, inserted by JB-191213
          if (plasmonmc) then ! use plasmon code
            pl_tslc = prop_distance(j) ! store current slice thickness
            pl_grid_sqx = 1.0_fp_kind / a0_slice(1,j) / ak1 ! store reciprocal grid sampling along x in rad/pixel
            pl_grid_sqy = 1.0_fp_kind / a0_slice(2,j) / ak1 ! store reciprocal grid sampling along y in rad/pixel
            pl_idum = idum ! let plasmon module take over the random seed
          end if
                    
          call qep_multislice_iteration(psi,prop(:,:,j),qep_grates(:,:,:,j),nopiy,nopix,&
                          &ifactory,ifactorx,idum,n_qep_grates,qep_mode,shift_arrayy,shift_arrayx)
                    
          if (output_probe_intensity) then
            k = (i-1)*n_slices+j
            if (output_cell_list(k)) then
              probe_intensity(:,:,cell_map(k)) = probe_intensity(:,:,cell_map(k)) + abs(psi)**2
            endif
          endif
          
          ! 2025-06-26, JB: Update accumulated SE transmission function
          if (SEI) then
              se_acctransm = se_acctransm * se_transmission(:,:,j)
              ! 2025-08-21, JB: added simulation of the secondary electron isotropic momentum taking
              !                 effect in the propagation to the detector and on the absorption
              ! convolute with the PSF of the current slice to model the point spread of
              ! the secondary electrons propagating with perpendicular momenta and backwards to the detector
              call conv2dr_apply(1+MODULO(j-1,n_se_psf), se_acctransm, se_acctransm, nopiy, nopix)
          endif
        enddo ! End loop j over slices
				
        !If this thickness corresponds to any of the output values then accumulate diffraction pattern
        if (any(i==ncells)) then
          z_indx = minloc(abs(ncells-i))
          
          if (arg_debug_wave>0) write(*,*) 'psi(0,0) = ', psi(1,1), ' at z=',z_indx(1)
                    
          ! Accumulate elastic wave function - this will be Fourier transformed later
          psi_elastic(:,:,z_indx(1)) = psi_elastic(:,:,z_indx(1)) + psi
                    
          !Transform into diffraction space
          call fft2(nopiy,nopix,psi,psi_out)
                    
          ! Accumulate diffaction pattern
          temp = abs(psi_out)**2
          if (arg_debug_intens>0) write(*,*) 'cbed(0,0) = ', temp(1,1), ' at z=',z_indx(1)
          cbed(:,:,z_indx(1)) = cbed(:,:,z_indx(1)) + temp

          if (ionization) then
            stem_ion_image(ny,nx,i_df,z_indx(1),:) = stem_ion_image(ny,nx,i_df,z_indx(1),:) +&
                          &sum(sum(ion_image,dim=2),dim=1)
          endif
            
        endif
                
      enddo ! End loop i over cells
			
      ! Plasmon scattering counting
      if (plasmonmc) call pl_populate(pl_exc_num)
			
    enddo !End loop over QEP passes
      
    intensity = sum(abs(psi)**2)
        
    do iz=1,nz
        
      ! Integrate the elastic diffraction pattern
      call fft2(nopiy,nopix,psi_elastic(:,:,iz),psi_out)
      if (stem) then
        stem_image(ny,nx,i_df,1:ndet,iz) = sum(sum(spread(cbed(:,:,iz),dim=3,ncopies=ndet)*&
                          &masks(:,:,:),dim=1),dim=1)
        if (arg_debug_stemdet>0) then
            write(*,*) 'cbed(0,0) = ', cbed(1,1,iz), ' at iz=',iz
            write(*,*) 'stem_image_1 = ', stem_image(ny,nx,i_df,1,iz), ' at iz=',iz
        end if
        stem_elastic_image(ny,nx,i_df,1:ndet,iz)= sum(sum(spread(abs(psi_out)**2,dim=3,ncopies =ndet)*&
                          &masks(:,:,1:ndet),dim=1),dim=1)
      endif
      if (EELS) then ! 2025-06-07, JB: EELS flag replaces previous .not.EDX conditions
        eels_correction_image(ny,nx,i_df,z_indx(1)) = sum(cbed(:,:,iz)*eels_correction_detector)
      endif

#endif            

      if (pacbed) then
        ! Output 4D STEM diffraction pattern
        if(fourDSTEM) then
        ! Output total (elastic and inelastic) diffraction pattern
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
          ! Output elastic only diffraction pattern
          filename = trim(adjustl(output_prefix))
          if (probe_ndf>1) filename = trim(adjustl(filename))//defocus_string(probe_df(i_df),lengthdf)
          if (nz>1) filename = trim(adjustl(filename))//'_z='//to_string(int(zarray(iz)))//'_A'
          filename = trim(adjustl(filename))//'_pp_'//to_string(nx)//'_'//to_string(ny)//'_Elastic_Diffraction_pattern'
          call binary_out_unwrap(nopiy, nopix, abs(psi_out)**2/n_qep_passes**2, filename,write_to_screen=.false.&
                          &,nopiyout=nopiyout,nopixout=nopixout)
        endif
      endif
    enddo ! iz
			
    if (output_probe_intensity) then
      call probe_intensity_to_file(probe_intensity,i_df,ny,nx,n_qep_passes,probe_ndf,nysample,nxsample)
    endif
    
    cntscan = cntscan + 1
    dnow = secnds(t1)
      
    enddo ! End loop over x probe positions
    enddo ! End loop over y probe positions
    enddo ! End loop over defocus series
    
! >> SCANNING LOOPS (code formatting relaxed due to too many loops) <<
! --------------------------------------------------------------------
    
    
	
    ! QEP normalisation
    if(ionization) then
      stem_ion_image = stem_ion_image/float(n_qep_passes)
      if(EELS) eels_correction_image = eels_correction_image/float(n_qep_passes) ! 2025-06-07, JB: EELS flag replaces previous .not.EDX conditions
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
    
!#ifdef GPU
!    ! debug output of inelastic intensities
!    if (ALLOCATED(inel_intens)) then
!        write(*,*) 'Creating dump of gathered inelastic intensities...'
!        open(unit=456,file='dump_inel_intens.bin',form='binary',status='replace',convert='big_endian')
!        write(456) inel_intens
!        close(456)
!        write(*,*) 'Dump of inelastic intensities written to dump_inel_intens.bin'
!        write(*,*) '  (target_atom, transition, slice, scan_y, scan_x)'
!        write(*,*) 'Shape: ', shape(inel_intens)
!        write(*,*)
!    endif
!    if (ALLOCATED(inel_intens_valid)) then
!        open(unit=457,file='dump_inel_intens_valid.bin',form='binary',status='replace',convert='big_endian')
!        write(457) inel_intens_valid
!        close(457)
!    endif
!#endif

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
          filename = trim(adjustl(fnam)) // '_'//trim(adjustl(substance_atom_types(atm_indices(ii))))//&
                            &'_'//trim(adjustl(Ion_description(ii)))//'_shell_EDX'
          call output_stem_image(stem_ion_image(:,:,:,:,ii), filename,probe_df)
        elseif(SEI) then
          filename = trim(adjustl(fnam)) // '_'//trim(adjustl(substance_atom_types(atm_indices(ii))))//&
                            &'_'//trim(adjustl(Ion_description(ii)))//'_shell_SE'
          call output_stem_image(stem_ion_image(:,:,:,:,ii)*se_det_scale, filename,probe_df)
        else ! EELS
          filename = trim(adjustl(fnam))// '_'//trim(adjustl(substance_atom_types(atm_indices(ii))))//&
                            &'_'//trim(adjustl(Ion_description(ii)))//'_orbital_EELS'
          call output_stem_image(stem_ion_image(:,:,:,:,ii), filename,probe_df)
          stem_ion_image(:,:,:,:,ii) = stem_ion_image(:,:,:,:,ii)*eels_correction_image
          filename =  trim(adjustl(fnam)) //'_'//trim(adjustl(substance_atom_types(atm_indices(ii))))//&
                            &'_'//trim(adjustl(Ion_description(ii)))//'_orbital_EELS_Corrected'            
          call output_stem_image(stem_ion_image(:,:,:,:,ii), filename,probe_df)
        endif
      enddo
      if (SEI) then ! 2025-06-07, combined SE image output
          filename = trim(adjustl(fnam)) // '_SE'
          call output_stem_image(sum(stem_ion_image(:,:,:,:,:),dim=5)*se_det_scale, filename, probe_df)
      endif
      if(EELS) then ! 2025-06-07, JB: EELS flag replaces previous .not.EDX conditions
        filename = trim(adjustl(fnam)) // '_EELS_CorrectionMap' 
        call output_stem_image(eels_correction_image, filename,probe_df)
      endif
    endif
    
    if(pacbed) then
      const = float(nysample*nxsample)
      do i=1,nz
        filename = trim(adjustl(output_prefix))
        if (n_tilts_total>1) fnam = trim(adjustl(output_prefix))//tilt_description(claue(:,ntilt),ak1,ss,ig1,ig2)
        if(nz>1) filename=trim(adjustl(filename))//'_z='//zero_padded_int(int(zarray(i)),length)//'_A'
        call binary_out_unwrap(nopiy,nopix,pacbed_pattern(:,:,i)/const,trim(adjustl(filename))//&
                            &'_PACBED_Pattern',nopiyout=nopiyout,nopixout=nopixout)
        if(output_thermal) then
          call binary_out_unwrap(nopiy,nopix,PACBED_elastic(:,:,i)/const,trim(adjustl(filename))//&
                            &'_elastic_PACBED_Pattern',nopiyout=nopiyout,nopixout=nopixout)    
          call binary_out_unwrap(nopiy,nopix,pacbed_pattern(:,:,i)/const-PACBED_elastic(:,:,i)/const,&
                            &trim(adjustl(filename))//'_thermal_PACBED_Pattern',nopiyout=nopiyout,nopixout=nopixout)
        endif
      enddo
    endif

#ifdef GPU

        
    if(tp_eels.and.istem) then ! EFISTEM
      rnorm = 1.0 / float(n_qep_passes*nysample*nxsample)
      do i=1,nz
        do l=1,imaging_ndf
          temp = efistem_image_d(:,:,l,i)
          call output_TEM_result(output_prefix,tile_out_image(temp,ifactory,ifactorx)*rnorm,&
                            &'energy_filtered_ISTEM',nopiy,nopix,manyz,imaging_ndf>1,manytilt,z=zarray(i),&
                            &lengthz=length,lengthdf=lengthimdf,tiltstring = tilt_description(claue(:,ntilt),&
                            &ak1,ss,ig1,ig2),df = imaging_df(l))
        enddo
      enddo
    endif	
		
    if(tp_eels) then ! DC-EELS
      rnorm = 1.0 / float(n_qep_passes)
      do l=1,dc_numeels
        if (single_channeling) then
            filename = trim(adjustl(fnam))//'_single_channeling_EELS_'//zero_padded_int(l,2)
        else
            filename = trim(adjustl(fnam))//'_double_channeling_EELS_'//zero_padded_int(l,2)
        endif
        call output_stem_image(Hn0_eels_dc(:,:,:,:,l)*rnorm, filename,probe_df) ! 2025-07-14, JB: QEP normalization added
      enddo
    endif

    if(istem) then ! ISTEM
      rnorm = 1.0 / float(n_qep_passes*nysample*nxsample)
      do i=1,nz
        do l=1,imaging_ndf
          temp = istem_image_d(:,:,l,i)
          call output_TEM_result(output_prefix,tile_out_image(temp,ifactory,ifactorx)*rnorm,&
                            &'ISTEM',nopiy,nopix,manyz,imaging_ndf>1,manytilt,z=zarray(i),&
                            &lengthz=length,lengthdf=lengthimdf,tiltstring = tilt_description(claue(:,ntilt),&
                            &ak1,ss,ig1,ig2),df = imaging_df(l))
        enddo
      enddo
    endif	
#endif
  enddo ! n_tilts_total
    
  call conv2dr_uninit() ! uninitialize the convolution module

end subroutine qep_stem
