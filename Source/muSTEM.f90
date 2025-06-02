!--------------------------------------------------------------------------------
!   Program: MU_STEM (GPU VERSION by GPU preprocessor flag with PGI fortran)
!                    (CPU VERSION by CPU preprocessor flag with Intel fortran)
!    
!   Description:    Calculate (S)TEM images and diffraction patterns using the
!                   multislice approach.
!                   Includes the contributions of TDS using either the Hall & Hirsch
!                   absorptive model or the Quantum Excitation of Phonons model.
!                   Both plane wave and convergent beam illumination may be used.
!                   STEM EDX/EELS images can be calculated within the local approximation.
!
!                   Outputs are big-endian floating point numbers in either
!                   32-bit or 64-bit precision, depending on which precision
!                   is chosen in mod_precision.f90.
!    
!   Maintainer:     Hamish Brown
!   Email:          hamish.brown@monash.edu
!   Date:           August 2017
!   Requirements:   PGI Fortran
!
!   version:        6.1  
!  
!  Copyright (C) 2025  L. J. Allen, H. G. Brown, A. J. D’Alfonso, S.D. Findlay, B. D. Forbes, J. Barthel
!
!  v6.1 includes modifications by J. Barthel and L. J. Allen (2019 - 2025)
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

    program MU_STEM
    
        use m_precision
        use m_user_input
        use global_variables
        use m_lens
#ifdef GPU        
        use cuda_setup, only: setup_GPU
        use cuda_array_library, only: set_blocks
#endif
        use output, only: setup_output_prefix,timing
        use m_tilt, only: prompt_tilt
        use m_absorption
        use m_potential
        use m_multislice
        use m_string
        use m_electron
        use m_Hn0
        
        implicit none
        
        integer :: i_illum, i_tds_model, i_cb_menu, i_cb_calc_type,ifile,nfiles,i_arg,idum,i,ieftem
        integer :: arg_num_threads, nerr
        
        logical :: nopause = .false.,there,ionization,stem,pacbed
        character(512)::command_argument
        character(120)::fnam, pnam
        character(2):: symb
        
        nerr                     = 0
        arg_num_threads          = 0
        pnam                     = "single"
        if(fp_kind==Double) pnam = "double"
108     write(6,109) trim(pnam)
        109     format(&
       &1x,'|----------------------------------------------------------------------------|',/,&
       &1x,'|              Melbourne University (scanning) transmission electron         |',/,&
       &1x,'|                            microscopy computing suite                      |',/,&
       &1x,'|      __       __  __    __   ______  ________  ________  __       __       |',/,&
       &1x,'|     |  \     /  \|  \  |  \ /      \|        \|        \|  \     /  \      |',/,&
       &1x,'|     | $$\   /  $$| $$  | $$|  $$$$$$\\$$$$$$$$| $$$$$$$$| $$\   /  $$      |',/,&
       &1x,'|     | $$$\ /  $$$| $$  | $$| $$___\$$  | $$   | $$__    | $$$\ /  $$$      |',/,&
       &1x,'|     | $$$$\  $$$$| $$  | $$ \$$    \   | $$   | $$  \   | $$$$\  $$$$      |',/,&
       &1x,'|     | $$\$$ $$ $$| $$  | $$ _\$$$$$$\  | $$   | $$$$$   | $$\$$ $$ $$      |',/,&
       &1x,'|     | $$ \$$$| $$| $$__/ $$|  \__| $$  | $$   | $$_____ | $$ \$$$| $$      |',/,&
       &1x,'|     | $$  \$ | $$ \$$    $$ \$$    $$  | $$   | $$     \| $$  \$ | $$      |',/,&
       &1x,'|      \$$      \$$  \$$$$$$   \$$$$$$    \$$    \$$$$$$$$ \$$      \$$      |',/,&
       &1x,'|                                                                            |',/,&
       &1x,"|       Copyright (C) 2025 L.J. Allen, H.G. Brown, A.J. D'Alfonso,           |",/,&
       &1x,'|              S.D. Findlay, B.D. Forbes, J. Barthel                         |',/,&
	   &1x,'|       email: hgbrown@unimelb.edu.au                                        |',/,&
       &1x,'|              ju.barthel@fz-juelich.de (for this version)                   |',/,&
       &1x,'|       This program comes with ABSOLUTELY NO WARRANTY.                      |',/,&
       &1x,'|                                                                            |',/,&
       &1x,'|       This program is licensed to you under the terms of the GNU           |',/,&
	   &1x,'|       General Public License Version 3 as published by the Free            |',/,&
	   &1x,'|       Software Foundation.                                                 |',/,&
       &1x,'|                                                                            |',/,&
#ifdef GPU
       &1x,'|       GPU Version 6.1 (branch https://github.com/ju-bar 2025-06-02)        |',/,&
  	   &1x,'|           ! absorptive model calculations are known to crash occasionally  |',/,&
#else
       &1x,'|       CPU only Version 6.1 (branch https://github.com/ju-bar 2025-06-02)   |',/,&
#endif
       &1x,'|           (',a6,' precision compile)                                       |',/,&
       &1x,'|                                                                            |',/,&
       &1x,'|       Optional arguments:  see muSTEM.exe options                          |',/,&
       &1x,'|                                                                            |',/,&
       &1x,'|----------------------------------------------------------------------------|',/)
       

        ! Process command line arguments
        do i_arg = 1, command_argument_count()
            call get_command_argument(i_arg, command_argument)
            idum = index(command_argument, 'omp_num_threads')
            if (idum > 0) then
                read(unit=command_argument(idum+16:),fmt='(i)',iostat=nerr) arg_num_threads
                cycle
            end if
            select case (trim(adjustl(command_argument)))
            case ('options')
                    write(*,*) 
                    write(*,*) '  List of options for muSTEM:    e.g.  muSTEM nopause'
                    write(*,*) '    nopause      - avoids program pauses'
                    write(*,*) '    ionic        - applies ionic form factors (requires charge input via xtl)'
                    write(*,*) '    linpoleels   - applies linear interpolation on EELS energy windows'
                    write(*,*) '    omp_num_threads={n} sets number of OpenMP threads, e.g. omp_num_threads=3'
                    write(*,*)
                    stop
                case ('nopause')
                    nopause = .true.
			    case ('timing')
				    timing = .true.
			    case ('ionic')
                    ionic = .true.
                case ('linpoleels')
                    linpoleels = .true.
                !case ('mmouse_wave')
                !    arg_debug_wave = 1
                !case ('mmouse_intens')
                !    arg_debug_intens = 1
                !case ('mmouse_stemdet')
                !    arg_debug_stemdet = 1
            end select
        enddo
        
        if (.not. nopause) then
            write(*,*) ' Press enter to continue.'
            read(*,*)
        endif
       
        
        
        ! Set up user input routines, nfiles is the number
        ! of user input files to play if "play all" is inputted
        nfiles = init_input()
       
        do ifile=1,nfiles
            !If play or play all open relevant user input file.
            if(input_file_number.ne.5) then
                fnam = get_driver_file(ifile)
                inquire(file=fnam,exist = there)
                if (there) then
                    open(unit=in_file_number, file=fnam, status='old')
                else
                    write(*,*) "Couldn't find user input file: ",trim(adjustl(fnam))
                    cycle
                endif
            endif
        ! Set up CPU multithreading
        call setup_threading(arg_num_threads)
#ifdef GPU
        ! Set up GPU
        call setup_GPU
		double_channeling =.true. !Double channeling possible with GPU calculation
#else
        open(6,carriagecontrol ='fortran')
		double_channeling =.false. !Double channeling not possible with CPU calculation
#endif
        
        call command_line_title_box('Dataset output')
        ! Ask for the prefix for dataset output
        call setup_output_prefix
        
        call command_line_title_box('Specimen details')
        ! Read in the xtl file
        call set_xtl_global_params 
        call validate_xtl
        
      
        ! Calculate the mean inner potential
        call set_volts(nt, atf, nat, atomf, volts, ss)       
        
        ! Set electron quantities
        call constants           
        
        ! Ask for slicing scheme
        call setup_slicing_depths
        
        ! Ask for thickness
        call setup_specimen_thickness
       
        ! Set up the potential calculation method
        call prompt_high_accuracy    
               
        ! Set the unit cell tiling and grid size
        call set_tiling_grid
        
        i_illum = -1
        do while(i_illum<1.or.i_illum>2)
            call command_line_title_box('Illumination')
            write(*,*) '<1> Plane wave      (including HRTEM images and diffraction patterns)'
            write(*,*) '<2> Convergent-beam (including STEM images and CBED patterns)'
                                                                                               
            call get_input('Calculation type', i_illum) 
            write(*,*)
        enddo
        
        ! Set illumination flags
        pw_illum = (i_illum == 1)
        cb_illum = .not. pw_illum
        

        
        ! For convergent-beam, set up the probe forming lens
        if (cb_illum) call setup_lens_parameters('Probe',probe_aberrations,probe_cutoff)

        
        i_tds_model=-1
        do while(i_tds_model<1.or.i_tds_model>2)
            call command_line_title_box('Thermal scattering model')
            write(*,*) '<1> Quantum Excitation of Phonons (QEP) model'
            write(*,*) '    (Accurately accounts for inelastic scattering'
            write(*,*) '     due to phonon excitation)'
            write(*,*) '<2> Absorptive model'
            write(*,*) '    (Calculates faster than QEP but only approximates'
            write(*,*) '     inelastic scattering due to phonon excitation)'

            call get_input('<1> QEP <2> ABS', i_tds_model)
            write(*,*)
        enddo
        
        ! Set TDS model flags
        complex_absorption = (i_tds_model == 2)
        qep = .not. complex_absorption
        
        ! Prompt for including absorptive potential
        if (complex_absorption) call prompt_include_absorption
        
        if((complex_absorption.and.include_absorption).or.qep) then
        write(*,*) ' Options for output of inelastically and elastically scattered components'
        write(*,*) ' <1> Only output total signal (ie that measured in experiment)'
        write(*,*) ' <2> Separately output elastic, inelastic and total signal'
        write(*,*) 'Note: option <2> does not apply to EELS or EDX output'
        if(complex_absorption.and.include_absorption) write(*,*) 'Note: option <2> only applies to STEM imaging'
        write(*,*)
        call get_input('Elastic and inelastic scattering output choice', i_tds_model)
        write(*,*)
        output_thermal = i_tds_model==2
        endif
         

        ! Prompt user for a tilt for either kind of illumination
        call prompt_tilt        
#ifdef GPU        
        ! Setup the CUDA thread hierachy for nopiy, nopix arrays
        call set_blocks(nopiy, nopix)
#endif
        ! Calculate the slices of the supercell
        call calculate_slices
        
        ! Calculate the bandwidth limiting matrix
        call make_bwl_mat                           
        
        ! Ask for QEP parameters
        if (qep) call setup_qep_parameters(n_qep_grates,n_qep_passes,nran,quick_shift,ifactory,ifactorx)
		
		! Ask plasmon scattering parameters
		if (qep) call setup_plasmon_parameters(ekv*1000.0_fp_kind, thickness, plasmonmc)
        
        ! Save/load transmission functions
        call prompt_save_load_grates
        
        ! Set up the imaging lens
        if (pw_illum) then
			 call setup_lens_parameters('Image',imaging_aberrations,imaging_cutoff)
#ifdef GPU
			write(*,*) 'Calculate an energy-filtered TEM (EFTEM) image? <1> Yes <2> No'
			call get_input('Do EFTEM? <1> Yes <2> No',ieftem)
			double_channeling = ieftem == 1
			if(double_channeling) call ionizer_init(.false.)
#endif
		else
			imaging_ndf = 1
		endif
        
        ! Choose the convergent-beam calculation type
        if (cb_illum) then
            
            call command_line_title_box('Calculation type')
            write(*,*) 'Choose a calculation type:'
115         write(*,*) '<1> CBED pattern'
            write(*,*) '<2> STEM (BF/ABF/ADF/EELS/EDX/PACBED/4D-STEM)'
            
            call get_input('<1> CBED <2> STEM/PACBED', i_cb_calc_type)
            write(*,*)
            
            select case (i_cb_calc_type)
                case (1)
                    ! CBED pattern
                    call place_probe(probe_initial_position)
					          double_channeling=.false.
                case (2)
                    ! STEM images
                    call STEM_options(STEM,ionization,PACBED,istem,double_channeling)
                    fourdSTEM= .false.
                    if(pacbed) call fourD_STEM_options(fourdSTEM,nopiyout,nopixout,nopiy,nopix)
                    call setup_probe_scan(PACBED.and.(.not.(ionization.or.STEM)))
                    call prompt_output_probe_intensity
                    if(STEM) call setup_integration_measurements
                    adf = STEM.and.(.not.qep)
					          if(istem) call setup_lens_parameters('Image',imaging_aberrations,imaging_cutoff)
                    
                    ! Precalculate the scattering factors on a grid
                    call precalculate_scattering_factors()
                    if(ionization) call setup_inelastic_ionization_types
					          if(double_channeling) call ionizer_init(.true.)
                    
                case default
                    goto 115
                    
            end select
            
            
        endif
        
        ! Start simulation
        if (pw_illum .and. qep) then
            ! Plane wave QEP
            call qep_tem

        elseif (pw_illum .and. complex_absorption) then
            ! Plane wave absorptive
            call absorptive_tem
           
        elseif (cb_illum .and. qep .and. i_cb_calc_type.eq.1) then
            ! Convergent-beam QEP CBED
            call qep_tem

        elseif (cb_illum .and. qep .and. i_cb_calc_type.eq.2) then
            ! Convergent-beam QEP STEM
            call qep_stem(STEM,ionization,PACBED)
            
        elseif (cb_illum .and. complex_absorption .and. i_cb_calc_type.eq.1) then
            ! Convergent-beam absorptive CBED
            call absorptive_tem

        elseif (cb_illum .and. complex_absorption .and. i_cb_calc_type.eq.2) then
            ! Convergent-beam absorptive STEM
            call absorptive_stem(STEM,ionization,PACBED)
            
        endif
         
        write(*,*)
        close(in_file_number)
        call reset_allocatable_variables
        
        enddo
        if (.not. nopause) then
            write(*,*) ' Press enter to exit.'
            read(*,*) 
        endif
    end program Mu_STEM

   
   
    subroutine setup_threading(nthreads)
    
        use m_string, only: to_string,command_line_title_box
    
        implicit none
      
        integer, intent(in) :: nthreads
        integer*4 :: num_cores, num_threads
        integer*4 :: omp_get_max_threads, omp_get_num_procs
    
        num_cores = omp_get_num_procs()
        if (nthreads > 0) then
            num_threads = min(nthreads, num_cores-1) ! Use the user-specified number of threads, capped to one less than the number of cores
        else
            num_threads = num_cores/2 ! Default to half the number of processors.
        end if
        num_threads = max(1, num_threads) ! Ensure at least one thread is used.
        
        call omp_set_num_threads(num_threads)
#ifdef GPU
#else
        call dfftw_init_threads()
        call dfftw_plan_with_nthreads(num_threads)
#endif
        call command_line_title_box('CPU multithreading')
        write(*,*) 'The number of available logical cores is: ' // to_string(num_cores)
        write(*,*) 'The number of threads being used on the CPU is: ' // to_string(num_threads)
        write(*,*)
    
    end subroutine
  
	  subroutine reset_allocatable_variables()
      
		!The downside of using global variables... :(
        use global_variables
		use m_lens
        use m_potential
        use m_multislice
#ifdef GPU
		use cuda_potential
#endif
        if(allocated(Kz)                      ) deallocate(Kz)
        if(allocated(claue)                   ) deallocate(claue)
		if(allocated(nat)                     ) deallocate(nat)    
		if(allocated(tau)                     ) deallocate(tau)    
		if(allocated(atf)                     ) deallocate(atf)    
        if(allocated(atomf)                   ) deallocate(atomf)
		if(allocated(fx)                      ) deallocate(fx)
		if(allocated(dz)                      ) deallocate(dz)   
		if(allocated(zarray)                  ) deallocate(zarray)   
		if(allocated(ncells)                  ) deallocate(ncells)   
		if(allocated(fz    )                  ) deallocate(fz)   
		if(allocated(fz_DWF)                  ) deallocate(fz_DWF)   
		!if(allocated(sinc)                    ) deallocate(sinc)   
		if(allocated(inverse_sinc)            ) deallocate(inverse_sinc)   
		if(allocated(substance_atom_types    )) deallocate(substance_atom_types)   
		if(allocated(outer                   )) deallocate(outer)   
		if(allocated(inner                   )) deallocate(inner)   
		if(allocated(a0_slice                )) deallocate(a0_slice          )    
		if(allocated(nat_slice               )) deallocate(nat_slice         )    
		if(allocated(nat_slice_unitcell      )) deallocate(nat_slice_unitcell)    
		if(allocated(tau_slice               )) deallocate(tau_slice         )    
		if(allocated(tau_slice_unitcell      )) deallocate(tau_slice_unitcell)    
		if(allocated(prop_distance           )) deallocate(prop_distance     )    
		if(allocated(depths                  )) deallocate(depths            )    
		if(allocated(ss_slice                )) deallocate(ss_slice          )  
		if(allocated(probe_df                )) deallocate(probe_df          )  
		if(allocated(imaging_df              )) deallocate(imaging_df          )  
		if(allocated(atm_indices             )) deallocate(atm_indices          )  
		if(allocated(ion_description         )) deallocate(ion_description         )  
		if(allocated(ionization_mu           )) deallocate(ionization_mu           )  
		if(allocated(fz_adf                  )) deallocate(fz_adf                  )  
		if(allocated(adf_potential           )) deallocate(adf_potential           )  
		if(allocated(ionization_potential    )) deallocate(ionization_potential    )  
		if(allocated(eels_correction_detector)) deallocate(eels_correction_detector)
#ifdef GPU
		if(allocated(ccd_slice_array		 )) deallocate(ccd_slice_array)
		if(allocated(Volume_array			 )) deallocate(Volume_array)
#endif
		end subroutine