    ! A collection of routines pertinent to calculating the Hn0's
    !
    ! To do:
    !      
    !gspline and gsplint are repeated functions that are ALMOST redundant
    !The function calls are slightly different to the spline and splint subroutiunes
    !in numerical_recipies.f90
    !I intend to replace the appropriate calls so to remove redundancy
    !-------------------------------------------------------------------------------- 
    
    
module m_Hn0

    use global_variables
    use m_precision
    use m_user_input
    use m_lens
	use m_string
    
    implicit none
    
    integer(4) :: target_atom
    integer(4) :: numwve = 100
    integer(4) :: nsplpts = 300 !continuum wavefunction para
    integer(4) :: nstates                        !number of atomic states used in a calculation
    integer(4), allocatable :: state_vector(:,:) !lists the states as, ml,lpr,ml_pr
    integer(4), allocatable :: state_use(:) !lists the states used in the calculation, 1 if used, 0 if not used
    !real(fp_kind)  ::  bwl_rad  !band width limiting radius (in q space)
    real(fp_kind) :: VA(10000)  !potential of the atom seen by the ejected electron
    real(fp_kind) :: Pnl(2000)  !Bound state wavefunction of the target electron to be ionized
    integer(4) :: lorb             !angular momentum quantum number of the bound state
    real(fp_kind) :: Z          !atomic number of the target atom
    real(fp_kind) :: eval       !Energy of the ejected electron
    integer(4) :: nrad       !number of points that the radial wavefunctions are computed on
    real(fp_kind) :: et         !threshold energy
    real(fp_kind), save  :: FAK(100)  ! factorials necessary for the cleb3j code, save state needed for OMP thread safety
                                                  ! must be run 
    complex(fp_kind) :: alpha_n  !constant for the current conversion factor
    real(fp_kind),allocatable :: qmin(:),gint(:,:,:,:),dgint(:,:,:,:),qpos(:,:)
    logical :: dc_eels_segments ! EELS segment flag
    integer*4 :: dc_numeels, dc_eels_nseg ! number of EELS detectors and segments
    real(fp_kind) :: dc_eels_seg_offset !azimthal offset for the EELS detector segments
    real(fp_kind),allocatable :: dc_eels_outer(:), dc_eels_inner(:) !EELS detector ranges
    real(fp_kind),allocatable :: Hn0_eels_detector(:,:,:) ! EELS detector data for the Hn0 calculation
    !complex(fp_kind) :: tmatrix(nopiy,nopix)
    
    !the code constructs transition matrices at the origin and shifts uses the 1D factorisation
    ! shift arrays for the atoms in the unit cell
    integer(4),allocatable :: natoms_slice_total(:) !the number of target atoms per slice
    complex(fp_kind),allocatable :: Hn0_shiftx_coord(:,:,:),Hn0_shifty_coord(:,:,:)
    !integer*4::ycells,xcells
    
    logical :: all_atoms
    real(fp_kind) :: ionize_fraction
    real(fp_kind) :: tpb                                  ! constant 2pi*bohr
	real(fp_kind),allocatable :: ion_tau(:,:)
    

	logical::efistem,eftem,hn0_calc
    contains 
      
    
    !--------------------------------------------------------------------------------
    !This subroutine initialises the Hn0 ionization calculation
    !It also reads in the atomic structure file 
    !(*_1s.orb,etc) from RCN and stores the bound state wavefunction
    !in array Pnl(2000) and the potential from which the continuum 
    !wavefunctions will be calculated in array VA(10000).
    !Both are defined on regular grids with a 0.0015 Bohr unit step
    !size.  The potential is in Rydberg units.
      
    subroutine ionizer_init(stem)
		
        use m_multislice, only:get_cbed_detector_info, make_detector
        use output
        implicit none
    
        integer(4) :: i,itype,j,order,comaindex,nd,idc
		logical::stem,detectors,outputdetectors
        real(fp_kind) :: str_filter
        character*100::dstring
        !select atom to be ionized
        !calculate the factorials required for the cleb code
        call fakred()
		
		call command_line_title_box('Ionization transition-potential menu')

		itype=-1
		do while(itype<1.or.itype>nt) 
	 40  write(*,610)
    610 format (' Enter the index for the type of atom for which ionizations occur:',/)
     45 write(*,612) (substance_atom_types(i),i,i=1,nt)
    612 format(1x,a6,3x,'index =',i3,/)

        call get_input("Index of atom to ionize ",itype)
        write(*,*)

		if(itype<1.or.itype>nt) write(*,615) itype
		enddo
	615 format('Index ',i3,' invalid, please try again') 
		
        
        target_atom=itype
        !atomic no of atom follows
        Z = ATF(1,itype)

     66 call numset()   !read in from orb file
        write(*,*) 'Enter the energy (in eV) above threshold for the ionization calculation.'
		write(*,*) 'Note: This calculates the differential cross-section at one energy, i.e.,'
		write(*,*) 'there is no integration over an energy window.'
		write(*,*)
        call get_input('Ionization energy',eval)
        write(*,*)
		
        !ycells=1
        !xcells=1
		!If STEM EELS get detectors
		if(stem) then
            dc_numeels=0
            dc_eels_nseg=1
            dc_eels_segments=.false.
            outputdetectors=.false.
            dc_eels_seg_offset=0.0_fp_kind
            
            ! single or double channeling
            write(*,*) "Calculating STEM-EELS with transition potentials can include the effect"
            write(*,*) "of further elastic and thermal-diffuse scattering after the ionization,"
            write(*,*) "aka. double channeling. Double-channeling effects become more important"
            write(*,*) "for thicker specimens, smaller EELS detectors and strong dynamical"
            write(*,*) "diffraction conditions. However, double-channeling requires much more"
            write(*,*) "additional computation."
            idc = -1
            do while (idc<0.OR.idc>1)
              write(*,*) "Do you want to do double-channeling or run a faster but approximate"
              write(*,*) "single-channeling calculation?"
              write(*,*) "Select  <0> for single-channeling  or  <1> for double-channeling."
              call get_input('<0> single-channeling or <1> double-channeling', idc)
              write(*,*)
            enddo
            single_channeling = (idc == 0) ! set single channeling switch
            if (single_channeling) iwf_filter_threshold = 0.0_fp_kind ! do not filter inelastic wave functions with single-channeling
            
            
            do while(dc_numeels<1) ! ensure at least one detector
                
			write(*,*) 'muSTEM allows you to calculate STEM-EELS images for EELS detectors of'
			write(*,*) 'different sizes simultaneously.'
            write(*,*) 'Enter the number of detectors in the diffraction plane:'
		    write(*,*) "To segment detectors input a comma ',' and then the number"
		    write(*,*) "of angular segments (eg. for 4 rings and 4 quadrants input '4,4')."
		    write(*,*) 'To output the detectors for inspection end the input with a '
		    write(*,*) "question mark ('?'), ie '4,4?'."
            write(*,*)
            call get_input('Number of EELS detectors', dstring)
            write(*,*)
            
            !Check if there is a ',' which indicates the user
		    !wants segmented detectors
		    comaindex = index(dstring,',')
		    dc_eels_segments = (comaindex.ne.0)
            
            !Check if there is a ?, which indicates the user would
		    !like to output the detectors
		    outputdetectors = (index(dstring,'?').ne.0)
		    !Remove the ? so that it doesn't cause problems later
		    if (outputdetectors) dstring = dstring(:index(dstring,'?')-1)
            
            !Read detector string
            if(dc_eels_segments) then ! there are segments
			    read(dstring(:comaindex),*) dc_numeels
			    read(dstring(comaindex+1:),*) dc_eels_nseg
                dc_numeels = dc_numeels*dc_eels_nseg
			    if(dc_eels_nseg>1) then 
				    write(*,*) 'Please input orientation offset for segmented detectors in degrees:'
				    call get_input('Segment orientation offset (degrees)', dc_eels_seg_offset)
				    dc_eels_seg_offset = dc_eels_seg_offset / 180_fp_kind * pi !Convert from degrees to mrad
                else
                    dc_eels_segments = .false. ! no segments or one, turn option off
			    endif
            else ! no segments
			    read(dstring,*) dc_numeels
			    dc_eels_nseg = 1
            endif
            
            if (dc_numeels<1) then
                write(*,*) 'Error: Number of EELS detectors must be at least 1.'
                write(*,*) 'Please try again.'
                write(*,*)
            endif
            
            enddo ! while(dc_numeels<1)
            
			if(allocated(dc_eels_outer)) deallocate(dc_eels_outer)
            if(allocated(dc_eels_inner)) deallocate(dc_eels_inner)
            nd = int(dc_numeels/dc_eels_nseg)
			allocate(dc_eels_outer(nd), dc_eels_inner(nd))
            
            call get_cbed_detector_info(nd,ak1,dc_eels_outer,dc_eels_inner,.true.)

		    if(outputdetectors) then
			    do i=1,nd
			    do j=1,dc_eels_nseg
				    if(dc_eels_nseg>1) call binary_out_unwrap(nopiy,nopix, &
                            & make_detector(nopiy,nopix,ifactory,ifactorx,ss,dc_eels_inner(i), &
                            &   dc_eels_outer(i),2*pi*j/dc_eels_nseg-dc_eels_seg_offset,2*pi/dc_eels_nseg), &
                            &   'EELSDetector'//to_string((i-1)*dc_eels_nseg+j))
				    if(dc_eels_nseg==1) call binary_out_unwrap(nopiy,nopix, &
                            & make_detector(nopiy,nopix,ifactory,ifactorx,ss,dc_eels_inner(i), &
                            &   dc_eels_outer(i)),&
                            &   'EELSDetector'//to_string((i-1)*dc_eels_nseg+j))
			    enddo
			    enddo
            endif
            
            if (.NOT.single_channeling) then ! iwf filter setup
                ! insert setup of Hn0-strength filter here
                write(*,*) "To speed up STEM-EELS double channeling, muSTEM can omit weak inelastic"
                write(*,*) "wave functions. The quality of the simulation may be reduced when this"
                write(*,*) "filter is too strong. Using zero filter strength completely deactivates it."
                write(*,*) "This filter analyses the overlaps between the probe wave function and the"
                write(*,*) "transition potentials of all atoms in a slice."
                str_filter = -1.0_fp_kind
                do while(str_filter < 0.0_fp_kind .or. str_filter >= 1.0_fp_kind)
                    write(*,*) "Please, enter the fraction of total overlap intensity that is removed"
                    write(*,*) "by omitting the weakest contributions (typical range: 0.01 to 0.1)."
                    call get_input("inelastic wavefunction filter strength",str_filter)
                    write(*,*)
                enddo
                iwf_filter_threshold = str_filter
            endif
        
		
!			if(ifactory>1.and.ifactorx>1) then ! deprecated with new filters: remove this section and the related code using ycells, xcells, all_atoms
!				write(*,140)
!                140 format('In STEM the probe is sufficiently localised such that only a fraction of atoms',/,&
!				         &'in the simulation supercell have significant probability of ionization.',/,&
!						 &'To speed up calculation with large supercells, you have the option of',/,&
!						 &'choosing the fraction of the supercell within which atoms will be ionized.',/&
!                         &'(Works only well when the number of tiles is big and the unit cell is small.)',/,&
!						 &'Please input the fraction of the supercell that will be ionized.')    
!				call get_input('Fraction of supercell to ionize',ionize_fraction)
!				write(*,*)
!            
!				!Round this number up so that an odd integer number of unit cells are ionized
!				ycells = ceiling(ionize_fraction*ifactory)
!				ycells = ycells + mod(ycells+1,2)
!				xcells = ceiling(ionize_fraction*ifactorx)
!				xcells = xcells + mod(xcells+1,2)
!                if (ycells>= ifactory .or. xcells>=ifactorx) then ! one of the two directions covers the supercell extent
!                    ! fall back to ionization of all atoms
!                    all_atoms = .true.
!                    ycells = ifactory
!                    xcells = ifactorx
!                    write(unit=*,fmt=161)
!161                 format('All atoms in the supercell will be ionised.')    
!                else
!                    all_atoms = .false.
!				    write(unit=*,fmt=160) xcells,ycells,ifactorx,ifactory
!160                 format(i0,'x',i0,' unit cells will be ionised out of ',i0,'x',i0,' total cells in the supercell.')    
!                endif
!				write(*,*)
!            endif
        endif
        
        !call make_state_vector()
        write(*,*) "Please input the maximum order of the transition potential to be calculated,"
        write(*,*) "with order = abs(l - l') the absolute difference between initial and final"
        write(*,*) "angular momentum quantum numbers. As a rule of thumb, transitions of order=1"
        write(*,*) "are sufficient for qualitative agreement with experiment, but inclusion of"
        write(*,*) "higher order quantum numbers, such as transitions of order=2 or more are"
		write(*,*) "necessary for quantitative agreement with experiment."
		call get_input("maximum order of transitions",order)
        write(*,*)
        
        !call fill_state_vector(order,lorb) ! old version omitting delta_l=0 transitions for order=1
        call fill_state_vector_version2(order,lorb) ! 2025-07-17, JB: new version including delta_l=0 transitions for all orders
        write(*,*)
        
        ! setup of Hn0-strength filter
        write(*,*) "To speed up calculations with transition potentials, muSTEM can omit weak"
        write(*,*) "transitions. The quality of the simulation may be reduced when this"
        write(*,*) "filter is too strong. Using zero filter strength completely deactivates it."
        str_filter = -1.0_fp_kind
        do while(str_filter < 0.0_fp_kind .or. str_filter >= 1.0_fp_kind)
            write(*,*) "Please, enter the fraction of total transition strength that is removed"
            write(*,*) "by omitting the weakest transitions (typical range: 0.001 to 0.1)."
            call get_input("transition filter strength",str_filter)
            write(*,*)
        enddo
        Hn0_filter_threshold = str_filter
        
        !allocate arrays for the hn0 calculations
        if(allocated(qmin)) deallocate(qmin)
        if(allocated(gint)) deallocate(gint)
        if(allocated(dgint)) deallocate(dgint)
        if(allocated(qpos)) deallocate(qpos)
        allocate(qmin(nstates))
        allocate(gint(0:numwve,nsplpts,-1:1,nstates))      
        allocate(dgint(0:numwve,nsplpts,-1:1,nstates))
        allocate(qpos(300,nstates))
    end subroutine
    
    
    
    !--------------------------------------------------------------------------------
    !   subroutine numset
    !   gets the potential and the radial wavefunction from the orb files
    !   as calculated by Cowan

    subroutine numset
    
        implicit none

        integer(4) :: nradtest,i
	    real(fp_kind) :: pmax
	    real(fp_kind) :: dpnl
	

        character(120) :: infile
        ! 2025-07-16, JB: Added AOM check to match implementation limits
        lorb=-1 ! initialise the angular momentum quantum number to be invalid
        do while(lorb.lt.0.or.lorb.gt.2) ! ensure a valid angular momentum quantum number is entered
        write(*,300)
    300 format (' Enter Input File   file.orb')
        call get_input("orb file",infile)
        write(*,*) 
        open(66,file=trim(adjustl(infile)),status='old')
        pmax=0.0_fp_kind
        nrad=0
        nradtest=0
        !read in threshold energy and angular momentum of bound state
        read(66,*) et
        read(66,*) lorb
        write(*,100) et,lorb
    100 format(' Ionization Threshold Energy: ',g15.8,' eV',/,' Angular momentum of bound state: ',I3)
        write(*,*) 'Reading bound state wave function...'
        if (lorb.lt.0.OR.lorb.gt.2) then ! handle invalid angular momentum quantum number
            write(*,*) 'Error: Invalid angular momentum quantum number.'
            write(*,*) 'muSTEM currently supports l=0 (s-orbital) up to l=2 (d-orbitals).'
            write(*,*) 'Please try again.'
            write(*,*)
            close(66) ! close the file before retrying
        endif
        enddo ! loop entering valid orb file
        do i=1,2000
            Pnl(i)=0.0_fp_kind
            read(66,*) dPnl
            if(abs(dPnl).lt.1.0e-30_fp_kind) then
                Pnl(i)=0.0_fp_kind
            else
                Pnl(i)=dpnl
            endif
            pmax=max(pmax,abs(Pnl(i)))
        enddo

        ! locate tail
        do i=2000,1,-1
            nrad=i
            if(abs(Pnl(i)).gt.0.002_fp_kind*pmax) exit
        enddo
    
        ! ensure nrad is even for simpson integration
    101 if(2.0_fp_kind*float(int(float(nrad)/2.0_fp_kind)).ne.float(nrad)) then
            nrad=nrad+1
            if(nrad.ge.2000) nrad=2000
            goto 101
        endif
        write(*,*) 'Reading potential...'
        write(*,*) 

        ! read in potentials
        do I=1,10000
            VA(i)=0.0_fp_kind
            read(66,*) VA(i)
        enddo

        close(66)

    end subroutine
    
    
    
    !-------------------------------------------------------------------------------   
    !  subroutine to establish the vector ml,lpr,ml_pr used for the final state
 
    subroutine make_state_vector
    
        use global_variables
        use m_user_input
    
        implicit none
    
        integer*4   order,switch,singlestate
        !integer*4   state_counter
     
        write(*,*) ' Do you want to ionize to a <1> single state <2> otherwise'
        call get_input('single state <1>',singlestate)
        write(*,*)
        
        if(singlestate.eq.1) then
              nstates=1               !set the number of states to one
    100       if(allocated(state_vector)) deallocate(state_vector)
              allocate(state_vector(nstates,3))
              if(lorb.gt.0) then
    101             write(*,*) 'Enter a value for an initial ml quantum number (-l,l)'
                    call get_input('input ml quantum number',state_vector(1,1))
                    if(iabs(state_vector(1,1)).gt. iabs(lorb)) then
                          write(*,*) 'Error unphysical quantum number please reenter' 
                          goto 101
                    endif
              elseif(lorb.eq.0) then
                    state_vector(1,1)=0
              endif
              write(*,*) "Enter the final state quantum numbers l', m_l'"
              call get_input("Final state l'",state_vector(1,2))
              call get_input("Final state m_l'",state_vector(1,3))
              write(*,*) 'Initial l:',lorb,' ml:',state_vector(1,1),'l_pr',state_vector(1,2),'ml_pr',state_vector(1,3)
              write(*,*) '<1> Continue'
              write(*,*) '<2> Choose again'
              call get_input('<1> accept quantum numbers',switch)
              if(switch.ne.1) goto 100
              
        elseif(singlestate.ne.1) then 
            
            write(*,*) "Please input the maximum order of the transition potential to be calculated,"
            write(*,*) "with order = abs(l - l') the absolute difference between initial and final"
            write(*,*) "angular momentum quantum numbers. As a rule of thumb, dipole transitions"
            write(*,*) "(order=1) are sufficient for qualitative agreement with experiment,"
            write(*,*) "but inclusion of higher order quantum numbers such as quadrupole transitions"
		    write(*,*) "(order=2) are necessary for quantitative agreement with experiment."
		    call get_input("maximum order of transitions",order)
            write(*,*)
        
            call fill_state_vector(order,lorb) 
            write(*,*)
        
        endif 
    
    
    end subroutine
    
    
    
!-------------------------------------------------------------------------------
!     Calculates the number of states for a given input l
!     and a given order for transitions
!-------------------------------------------------------------------------------
    ! This version does not exclude transitions with l=l_p for order=1.
    subroutine fill_state_vector_version2(order, l_input)
    
        implicit none

        integer(4) ::   ml,l_p,ml_p,l_input
        integer(4) ::   l_p_max,ml_max,ml_p_max
        integer(4) ::   order,i,j
        integer(4) ::   state_vector_temp(10000,3)
    
        ml_max=l_input ! this the initial state angular momentum quantum number, sets max initial ml
        l_p_max=l_input+order ! maximum l' quantum number for the final state, order the transition order (delta l)
    
    
        ml=-l_input !initialise the first ml of the initial state
        j=1   !start the counter on the state_vector
        do while(ml.le.ml_max)                       ! limit sum over ml
            l_p=l_input-order                        ! initialise l' quantum number to l - order
            if(l_p.lt.0) l_p=0                       ! limit l' to 0
            do while(l_p.le.l_p_max)                 ! loop over l'
                ml_p=-l_p                            ! initialise ml' loop counter
                ml_p_max=l_p                         ! set max ml' for loop counter
                do while(ml_p.le.ml_p_max)           ! loop over ml'
                    state_vector_temp(j,:)=[ml,l_p,ml_p]
                    j=j+1                            ! increment state_vector counter
                    ml_p=ml_p+1                      ! increment ml' loop counter
                    if(j.gt.9999) PAUSE 'Too many states (out of vector bounds), program will HALT'
                enddo
                l_p=l_p+1                            ! increment l' loop counter
            enddo
            ml=ml+1                                  ! increment ml to next state
        enddo
        nstates=j-1                                  ! set the number of states
        
        if (nstates==0) then
            write(*,*) 'Error: No states found for the given input parameters.'
            write(*,*) 'Please check the input values for l and order.'
            stop 1
        endif
    
        if(allocated(state_vector)) deallocate(state_vector)
        allocate(state_vector(nstates,3))
        do i = 1,nstates
              state_vector(i,:) = state_vector_temp(i,:)
			  write(*,200) i,nstates,lorb,state_vector(i,:)
        enddo
		write(*,*)
    200 format('Transition ',i3,'/',i3,': l = ',i3,' ml = ',i3," l' = ",i3," ml' = ",i3,"  ")
    end subroutine
    
    ! This version does exclude transitions with l=l_p for order=1 but not for order>1.
    subroutine fill_state_vector(order, l_input)
    
        implicit none

        integer(4) ::   ml,l_p,ml_p,l_input
        integer(4) ::   l_p_max,ml_max,ml_p_max
        integer(4) ::   order,i,j
        integer(4) ::   state_vector_temp(10000,3)

    
        ml_max=l_input
        l_p_max=lorb+order ! don't understand why lorb and not l_input is used here
    
    
        ml=-l_input !initialise the first ml
        j=1   !start the counter on the state_vector
        if(order.gt.1) then                                   !Quadrupole and higher transitions 
              do while(ml.le.ml_max)                          !limit sum over ml
                     l_p=l_input-order
                     if(l_p.lt.0) l_p=0                       !reset l_p for next ml      
                     do while(l_p.le.l_p_max)                 !loop over lp
                          ml_p=-l_p                           !initialise ml_p loop counter
                          ml_p_max=l_p                        !set max ml_p for loop counter
                          do while(ml_p.le.ml_p_max)          !loop over ml_p
                                state_vector_temp(j,:)=[ml,l_p,ml_p]
                                j=j+1
                                ml_p=ml_p+1
								
                          enddo
                          l_p=l_p+1
						  
                     enddo
                     ml=ml+1
              enddo
		    nstates=j-1
        elseif(order.eq.1) then                               !Dipole transition case
              do while(ml.le.ml_max)
                     l_p=lorb-order
                     if(l_p.lt.0) l_p=0
                     do while(l_p.le.l_p_max)                 
                     !-----------------------------------------
                     !Special case for dipole transitions
                          if(l_p.eq.lorb) then
                                l_p=l_p+1    
                               cycle
                          endif
                     !-----------------------------------------
                          ml_p=-l_p
                          ml_p_max=l_p
                          do while(ml_p.le.ml_p_max)
                                state_vector_temp(j,:)=[ml,l_p,ml_p]
								!write(*,*) state_vector_temp(j,:)
                                j=j+1
                                ml_p=ml_p+1
                                if(j.gt.9999) PAUSE 'Too many states will crash now'
                          enddo
                          l_p=l_p+1
                     enddo
                     ml=ml+1
              enddo
		    nstates = j-1
        endif 
    
        ! nstates=j-1
    
        if(allocated(state_vector)) deallocate(state_vector)
        allocate(state_vector(nstates,3))
        do i = 1,nstates
              state_vector(i,:) = state_vector_temp(i,:)
			  write(*,200) i,nstates,lorb,state_vector(i,:)
        enddo
		write(*,*)
    200 format('Transition ',i3,'/',i3,': l = ',i3,' ml = ',i3," l' = ",i3," ml' = ",i3,"  ")
    end subroutine 
    

      
    !--------------------------------------------------------------------------------
    !
    subroutine fakred
        !  Calculate factorials required for CLEB code 
        !  JB: FAK is declared in this module
        !      Declaration is using save state
        !      to achieve thread-safety in OpenMP
        !
        !      This calculates FAK(N) = GAMMA(N)*10**(-N+1)
        !                               GAMMA(N) = (N-1)!
    
        implicit none
        
        real(fp_kind) :: AI
        integer(4) :: I

        FAK(1) = 1.0_fp_kind
        do I = 1, 99
            AI = I
            FAK(I+1) = FAK(I)*0.1_fp_kind*AI
        enddo
        
    end subroutine
    
    

    real(fp_kind) FUNCTION CLEB(A1,A2,A3,B2,B3)
        !     The Clebsch--Gordan coefficients 
        !
        !              Cleb = <a1, a2, b1, b2 | a3, b3>
        !     
        !     Uses FAK(N) = GAMMA(N)*10**(-N+1)  formed by
        !     a call to routine `FAKRED" that must be inserted
        !     in the main program of the code top use these CG's
        !           FAK is defined at the start of this file
    
        implicit none
    
        real(fp_kind) :: del
        real(fp_kind) :: A1,A2,A3
        real(fp_kind) :: B1,B2,B3
        real(fp_kind) :: TEST1,TEST2,TEST3
        real(fp_kind) :: ss_cleb
        integer(4) :: K,K1,K2,K3,K4,K5,K6
        integer(4) :: J01,J02,J03,J04
        integer(4) :: J1,J2,J3,J4,J5,J6,J7,J8
        integer(4) :: M,MO,MU

        real(fp_kind) :: SU
    
        CLEB = 0.0_fp_kind                                
        B1 = B3 - B2           
        IF ( ABS ( B1 ) .GT. A1 .OR. ABS ( B2 ) .GT. A2 .OR. ABS ( B3 ) .GT. A3      ) RETURN     
        TEST1 = A1 + B1 + 100.0 - float(INT ( A1 + B1 + 100.0 ))                 
        TEST2 = A2 + B2 + 100.0 - float(INT ( A2 + B2 + 100.0 ))      
        TEST3 = A1 + A2 + A3    - float(INT ( A1 + A2 + A3    ))          
        IF ( ( TEST1 .GT. 0.1_fp_kind ) .OR. ( TEST2 .GT. 0.1_fp_kind) .OR. ( TEST3 .GT. 0.1_fp_kind )      ) RETURN 
        K1 = A3 - A2 + B1 + 100.1_fp_kind                              
        K1 = K1 - 100                                   
        K2 = A3 - A1 - B2 + 100.1_fp_kind                     
        K2 = K2 - 100           
        K3 = A1 + A2 - A3 + 2.001_fp_kind      
        K4 = A1 - B1 + 2.001_fp_kind
        K5 = A2 + B2 + 2.001_fp_kind
        K6 = A1 + B1 + 0.001_fp_kind
        J01 = K3 - 1                                              
        J02 = K2 + K5 - 1                                     
        J03 = K1 + K4 - 1                             
        J04 = J02 + K4 + K6 - 1               
        IF ( MIN(J01,J02,J03) .LE. 0 ) RETURN       
        DEL = 0.1_fp_kind * FAK(J01) / FAK(J04) * FAK(J02) * FAK(J03)       
        MO = MIN(K3,K4,K5) - 1                     
        MU = MIN(K1,K2)        
        M = 1   
        IF ( MU ) 1, 2, 2        
      1 M = -MU + 1
      2 SU = 0.0_fp_kind
        IF ( M .GT. MO ) RETURN                             
        DO K = M, MO                             
           J1 = K1 + K                          
           J2 = K2 + K                 
           J3 = K3 - K                
           J4 = K4 - K           
           J5 = K5 - K          
           SU = -SU + 1. / ( FAK(K) * FAK(J1) * FAK(J2) * FAK(J3)* FAK(J4) * FAK(J5) )  
        END DO
        J6 = K2 + K3 - 1                                             
        J7 = K2 + K4 - 1                                          
        J8 = K1 + K5 - 1                           
        SS_cleb = -1.0_fp_kind                                         
        IF ( MOD(MO,2) .EQ. 1 ) SS_cleb = 1.0_fp_kind
        CLEB = SS_cleb * SU * SQRT ( DEL * FAK(K6+1) * FAK(K4-1) * FAK(K5-1)* FAK(J6) * FAK(J7) * FAK(J8) * ( 2.0_fp_kind * A3 + 1.0_fp_kind ) )
        IF ( ABS ( CLEB ) .LT. 1.0e-4_fp_kind ) CLEB = 0.0_fp_kind
        
    END function CLEB
    
    
    
    real(fp_kind) function plgndr(l_temp,m,x)
    
        implicit none
    
        integer(4) :: l_temp,m
        integer(4) :: i,ll
        real(fp_kind) :: x
        real(fp_kind) :: fact,pll,pmm,pmmp1,somx2
    
        if(m.lt.0.or.m.gt.l_temp.or.abs(x).gt.1.0_fp_kind)pause ' ERROR: Bad arguments in plgndr'
    
        pmm=1.0_fp_kind
        if(m.gt.0) then
            somx2=sqrt((1.0_fp_kind-x)*(1.0_fp_kind+x))
            fact=1.0_fp_kind
            do i=1,m
                pmm=-pmm*fact*somx2
                fact=fact+2.0_fp_kind
            enddo
        endif
        if(l_temp.eq.m) then
            plgndr=pmm
        else
            pmmp1=x*(2*m+1)*pmm
            if(l_temp.eq.m+1) then
                plgndr=pmmp1
            else
                do ll=m+2,l_temp
                    pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m)
                    pmm=pmmp1
                    pmmp1=pll
                enddo
                plgndr=pll
            endif
        endif
    
    end function
      
    
    
!--------------------------------------------------------------------------------

    complex(fp_kind) function sph_harm(l_val,m,theta,phi)
    implicit none
    
    integer(4) :: l_val,m,m1,flag
    real(fp_kind) :: theta,phi
    real(fp_kind) :: const
    
    if (m.lt.0) then
       flag = 1
       m1 = -m
    else
       m1 = m
    endif
    
    const = float(2*l_val+1) / (4.0_fp_kind*pi)
    const = const * fac(l_val-m1) / fac(l_val+m1)
    const = sqrt(const)
    const = const * plgndr(l_val,m1,cos(theta))
    sph_harm = const * exp(ci*float(m1)*phi)
    
    if (flag.eq.1) then
        !sph_harm = conjg(sph_harm) * (-1**mod(m,2))
        !BDF FIXED 2016-07-22 
        sph_harm = conjg(sph_harm) * (-1)**mod(m,2)
    endif
    
    return
    end function
    
    complex(fp_kind) function sph_harm2(l_val,m,theta,phi) result(sph_harm)
    implicit none
    
    integer(4) :: l_val,m
    real(fp_kind) :: theta,phi
    real(fp_kind) :: const
    
    const = float(2*l_val+1) / (4.0_fp_kind*pi)
    const = const * fac(l_val-abs(m)) / fac(l_val+abs(m))
    const = sqrt(const)
    const = const * plgndr(l_val,abs(m),cos(theta))
    sph_harm = const * exp(ci*abs(m)*phi)
    
    if (m.lt.0) then
        sph_harm = conjg(sph_harm) * (-1)**abs(m)
    endif
    
    return
    end function
    
    !--------------------------------------------------------------------------------
    !     Numerical recipies factorial function
    !
    real(fp_kind) function fac(n)
    implicit none
    integer(4) :: n
    integer(4) :: j,ntop
    real(fp_kind) :: a(33)
    SAVE ntop,a
    DATA ntop,a(1)/0,1.0_fp_kind/
    
    if (n.lt.0) then
          pause 'negative factorial in factrl'
    else if (n.le.ntop) then
          fac=a(n+1)
    else if (n.le.32) then
          do j=ntop+1,n
                a(j+1)=j*a(j)
          enddo    
          ntop=n
          fac=a(n+1)
    else
          fac=exp(gammln(float(n)+1.0_fp_kind))
    endif
    return
    END function

    !-------------------------------------------------------------------------------
    real(fp_kind) FUNCTION gammln(xx)
    implicit none
    INTEGER(4) :: j
    real(fp_kind) :: xx
    real(fp_kind) :: ser,stp,tmp,x,y,cof(6)
    SAVE cof,stp
    DATA cof,stp/76.18009172947146_fp_kind,-86.50532032941677_fp_kind,24.01409824083091_fp_kind,-1.231739572450155_fp_kind,.1208650973866179e-2_fp_kind,-.5395239384953e-5_fp_kind,2.5066282746310005_fp_kind/
    x=xx
    y=x
    tmp=x+5.5_fp_kind
    tmp=(x+0.5_fp_kind)*log(tmp)-tmp
    ser=1.000000000190015_fp_kind
    do 11 j=1,6
      y=y+1.0_fp_kind
      ser=ser+cof(j)/y
11  continue
    gammln=tmp+log(stp*ser/x)

    return
    END function
    
    
    
    
    
    
    
    !---------------------------------------------------------------------
    !   subroutine calc_gsplint calculates the overlap integral 
    !   between the continuum and the bound state radial wavefunctions
    !   
    !   input:  l,ml    initial state
    !           lpr     final state angular momentum
    !                   *note that the m_lpr is used when filling the tmatrix
    !           eval    energy of the ejected electron
    !
    !   output: 
    !subroutine calc_gsplint(VA,Pnl,nrad,Z,l,ml,lpr,mlpr,qpos,qmin,gint,dgint)
    subroutine calc_gsplint(ml,lpr,mlpr,qpos_local,qmin_local,gint_local,dgint_local)
    
    use m_electron
    use m_crystallography, only: trimr
        
    implicit none
    
    integer(4) :: ml,lpr,mlpr   !l is in the module
    integer(4) :: lammax,lprmax!,lammin
    integer(4),parameter :: mono = 1 !flag to kill the monopole contribution (for test purposes, only mono = 1 gives physical results)
        ! Why is the monopole contribution killed?
        ! The reason is that we are calculating with partially nonorthogonalized
        ! wavefunctions. The monopole term (lam=0) can give rise to unphysical
        ! contributions in this case. By removing it here, we ensure that only
        ! physically meaningful multipole contributions are considered.
        ! The constant monopole term is expected to contribute zero in case of orthogonal
        ! sets of wavefunctios. <phi_f|1|psi_i> = 0 if f != i
    
    !real(fp_kind) :: Z             !atomic number
    real(fp_kind) :: ztmp   = 1     !flag for wavfuncsolo subroutine (1 for non hydrogenic potential)
    real(fp_kind) :: kappa          !wavevector and energy of the ejected electron (energy above threshold)
    real(fp_kind) :: eloss          !Total energy lost by the fast electron
    real(fp_kind) :: eka,akd        !Total energy (in kev) and wavevec of the scattered fast electron respectively
    real(fp_kind) :: qmin_local     !momentum transfer between ejected and incident electron
    real(fp_kind) :: qminr          !qmin*two pi times the born radius
    real(fp_kind) :: qmax           !maximum reciprocal space vector (i.e. q_min*nopix/2) magnitude
    real(fp_kind) :: qmaxr          !qmax+qmin*two pi times the born radius (i.e. maximum momentum transfer)
    real(fp_kind) :: qrange         !range of the momentum transfers
    real(fp_kind) :: delq           !deltaq, qrange/3150 works out the q spacing to evaluate the sbessel's on
                                    !of the sbessels
                                    !for each block 50 in the overlap integral the spacing doubles for evaluation
    real(fp_kind) :: qpos_local(300) !keeps a record of the qpositions that the sbessels are evaluated on.   
    real(fp_kind) :: gint_local(0:numwve,nsplpts,-1:1) !interpolation data for the g integrals
    real(fp_kind) :: dgint_local(0:numwve,nsplpts,-1:1) !derivative calculated by gspline
    real(fp_kind) :: g(3)           !test for the largest gvector 
 !   real(fp_kind) :: VA(10000),Pnl(2000)    !Potential and target wavefunction
    real(fp_kind) :: sbess(2000,0:numwve) !spherical bessels
    real(fp_kind) :: Rel(0:2000)          !continuum wavefunction calculated by wavfunc_solo
    !integer(4)   nrad
    
    !counters 
    real(fp_kind) :: qtmp            !qvalue for the bessel evaluations
    integer(4) :: iqblk,iq,iqstep !
    integer(4) :: lpr_temp
      
    !calculate the factorials required for the cleb code
    !call fakred()
    
    tpb = tp * bohr !two pi * bohr radius
    sbess = 0.0_fp_kind ! init bessel functions
    ! determine maximum reciprocal lattice vector
    g(1:3) = float(ig1(1:3))*float(nopix)/float(2*ifactorx) &
            & + float(ig2(1:3))*float(nopiy)/float(2*ifactory)
    qmax = trimr(g,SS) ! maximum reciprocal space vector in 1/A
    
    lammax = lorb + lpr ! maximum order included
    !lammin = abs(lorb - lpr) ! minimum order required (not used, we always start from lam=0)
    lprmax = lammax ! maximum order copy
    kappa = wavev(eval)                 ! wave number of the ejected electron 1/A
    
    call wavfunc_solo((eval/RYD),lpr,Ztmp,Rel) ! calculate continuum wave function Ztmp= 1 for non hydrogenic potential
    Rel = Rel * sqrt(bohr*kappa/(2.0_fp_kind*eval)) ! Wavefunction in eV 
140 eloss = eval + et                   ! energy loss in eV = energy of the ejected electron + threshold energy
    eka = 1000.0_fp_kind * ekv - eloss  ! energy of the scattered electron in eV
    if (abs(eka).le.0.01_fp_kind) then 
        akd = 0.0_fp_kind               ! ionization onset case
    else
        akd = wavev(eka)                ! wavenumber of the scattered electron 1/A
    endif
    qmin_local = ak - akd               ! qmin = momentum transfer q=k-k' between 
                                        !        the incident and scattered electrons 
    ! q***r variables below are in units of k * 2 pi * bohr radius (no dimension)
    qminr = qmin_local * tpb            ! determine maximum and minimum values for spline of G integrals
    if(qminr.lt.0.001_fp_kind) qminr = 0.001_fp_kind ! limit qmin
    qmaxr = (qmin_local+qmax) * tpb     ! adjust qmaxr
    qrange = qmaxr - qminr              ! set qr range
    delq = qrange / 3150.0_fp_kind      ! set qr step of range / (1000 pi)
    qtmp = qminr - delq                 ! initialize qr for calculation of sbessels
    do iqblk=1, 6                       ! run in 6 blocks of ...
        do iq=1, 50                     !  .. 50 q steps (in total 300 q steps)
            qtmp = qtmp + delq          ! increment qtmp
            call sbessel(qtmp,lprmax,sbess) ! Calculate the Sbessels at qtmp
            if(mono.eq.1) sbess(:,0)=sbess(:,0)-1.0_fp_kind ! kill the spurious monopole contribution
            ! Why is the monopole contribution killed here?
            ! The reason is that we are calculating with partially nonorthogonalized
            ! wavefunctions. The monopole term (lam=0) can give rise to unphysical
            ! contributions in this case. By removing it here, we ensure that only
            ! physically meaningful multipole contributions are considered.
            ! The constant monopole term is expected to contribute zero in case of orthogonal
            ! sets of wavefunctios. <phi_f|1|psi_i> = 0 if f != i
            iqstep = (iqblk-1)*50 + iq 
            qpos_local(iqstep) = qtmp   ! store value of qtmp
            call ovlapint_solo_new(lpr,iqstep,Pnl,Rel,gint_local,sbess) ! calculate overlap integrals
        enddo
        delq=delq+delq
    enddo
    ! call spline to set up cubic spline (these use the triange relation to determine
    !   the values of lambda that are non zero, an good approximation for spherically
    !   symmetric scattering, ie. we do not consider directional bonding here.
    !   the triangle relations are from the angular part)
    !   * |l-l'| <= lam <= l+l' for the integral to be non zero, 
    !        note that we start from lam=0 here, we do too much work? Yes, but only for one case
    !   * some lam values will fail parity check and are also excluded here
    !   * handles all cases well up to l=2 (d orbitals)
    do lpr_temp=0,lprmax ! loop over lam from 0 to l + l'
        ! S orbital calculation
        if(lorb.eq.0) then ! only lam = lpr
            call gspline(qpos_local,nsplpts,lpr_temp,0,gint_local,dgint_local)
        endif
        ! P orbital calculation
        if(lorb.eq.1) then ! only lam = lpr-1, lpr+1 (lam=lpr fails parity since lorb is odd)
            call gspline(qpos_local,nsplpts,lpr_temp,1,gint_local,dgint_local)
            call gspline(qpos_local,nsplpts,lpr_temp,-1,gint_local,dgint_local)
        endif
        ! D orbital calculation
        if(lorb.eq.2) then ! only lam = lpr-2, lpr, lpr+2 (other lam fail parity)
            call gspline(qpos_local,nsplpts,lpr_temp,1,gint_local,dgint_local)
            call gspline(qpos_local,nsplpts,lpr_temp,0,gint_local,dgint_local)
            call gspline(qpos_local,nsplpts,lpr_temp,-1,gint_local,dgint_local)
        endif
        ! F orbitals are not supported by this code. It would require more than 3
        ! slots in the multipole expansion.
    enddo   !end loop over lpr
    return
    end subroutine
    
    ! --------------------------------------------------------------------------------
    ! This function calculates the integral of the product of the arrays Rel, sbess_col and Pnl
    ! using Simpson's rule with a grid spacing of h.
    ! This is used for the overlap integral calculation.
    ! The number of radial points nrad is assumed to be even and a module variable.
    function simpsonrad(Rel, sbess_col, Pnl, h) result(integral)
        real(fp_kind), intent(in) :: Rel(0:), sbess_col(:), Pnl(:)
        real(fp_kind), intent(in) :: h
        real(fp_kind) :: integral

        integer :: I
        real(8) :: sumod, sumev, sumd
        integer :: nstpev, nstpod

        sumod = 0.0d0
        sumev = 0.0d0
        nstpev = 2
        nstpod = 1

        ! Even-indexed points
        do I = nrad - 2, nstpev, -2
            sumd = real(Rel(I) * sbess_col(I) * Pnl(I), kind=8) ! Use double precision for better accuracy
            sumev = sumev + sumd
        end do

        ! Odd-indexed points
        do I = nrad - 1, nstpod, -2
            sumd = real(Rel(I) * sbess_col(I) * Pnl(I), kind=8) ! Use double precision for better accuracy
            sumod = sumod + sumd
        end do

        integral = h * real((4.0d0 * sumod + 2.0d0 * sumev) / 3.0d0, fp_kind) ! Convert back to single precision
    end function

    
    
    !--------------------------------------------------------------------------------
    ! This subroutine calculates the overlap integral for a given lpr and iqstp.
    ! This version is shorter but does the same as the original.
    ! The simpsonrad function is used to calculate the integral.
    subroutine ovlapint_solo_new(lpr, iqstp, Pnl, Rel, gint_local, sbess)
      implicit none
    
      integer(4), parameter :: numwve = 100, nsplpts = 300
      integer(4) :: lpr, iqstp
      real(fp_kind), parameter :: h = 0.0015_fp_kind
      real(fp_kind) :: Pnl(2000), Rel(0:2000)
      real(fp_kind) :: gint_local(0:numwve, nsplpts, -1:1)
      real(fp_kind) :: sbess(2000, 0:numwve)
      integer(4) :: midx, bessel_idx
    
      do midx = -1, 1
        bessel_idx = lpr + midx * lorb ! read as lpr + dlam with dlam = midx * lorb
        ! Exclude cases where the bessel_idx is out of bounds
        if (bessel_idx < 0 .or. bessel_idx > numwve) cycle
    
        ! For lorb=0, only dlam=0 allowed
        if (lorb == 0 .and. midx /= 0) cycle
    
        ! For lorb=1, dlam=0 typically omitted (fails parity check)
        if (lorb == 1 .and. midx == 0) cycle
    
        gint_local(lpr, iqstp, midx) = simpsonrad(Rel, sbess(:, bessel_idx), Pnl, h)
      end do
    
    end subroutine ovlapint_solo_new


    ! Old version of the overlap integral calculation with inline simpson integration.
    ! changed summation to run with double precision
    !--------------------------------------------------------------------------------
    ! This subroutine calculates the integral of the product of the arrays Pel and
    ! Pnl with a spherical bessel funtion sbess using simpsons rule and a grid 
    ! spacing of h=0.0015 a.u.
    subroutine ovlapint_solo(lpr,iqstp,Pnl,Rel,gint_local,sbess)
    
    integer(4) numwve,nsplpts
    
    parameter(numwve=100,nsplpts=300)
    
    integer(4) :: lpr     !angular momentum of continuum state
    !integer(4) :: l	!angular momentum of bound state
    integer(4) :: iqstp
    real(fp_kind) :: h    !grid spacing
    real(8) :: sumod(-1:1),sumev(-1:1) 
    real(8) :: sumd
    integer(4) I
    integer(4) nstpev,nstpod
    !integer(4) nrad
    real(fp_kind) :: Pnl(2000)
    
    real(fp_kind) :: Rel(0:2000)
    real(fp_kind) :: gint_local(0:numwve,nsplpts,-1:1)
    real(fp_kind) :: sbess(2000,0:numwve) !spherical bessel functions
    
    h=0.0015_fp_kind
    
    sumev = 0.0d0
    sumod = 0.0d0
    
    nstpev=2
    nstpod=1
    
    ! sum over even terms
    do I=nrad-2,nstpev,-2
    
          if(lorb.eq.0) then
                sumd=real((Rel(I))*(sbess(I,lpr))*(Pnl(I)), kind=8)
                sumev(0)=sumev(0)+(sumd)
          endif
    
          if(lorb.eq.1) then
                if(lpr.ne.0) then
                      sumd=real((Rel(I))*(sbess(I,lpr-1))*(Pnl(I)), kind=8)
                else
                      sumd=0.0d0
                endif
                sumev(-1)=sumev(-1)+(sumd)
                sumd=real((Rel(I))*(sbess(I,lpr+1))*(Pnl(I)), kind=8)
                sumev(1)=sumev(1)+(sumd)
          endif
    
          if(lorb.eq.2) then
                if(lpr.gt.1) then
                      sumd=real((Rel(I))*(sbess(I,lpr-2))*(Pnl(I)), kind=8)
                else
                      sumd=0.0d0
                endif
                sumev(-1)=sumev(-1)+(sumd)
    
                sumd=real((Rel(I))*(sbess(I,lpr))*(Pnl(I)), kind=8)
                sumev(0)=sumev(0)+(sumd)
                sumd=real((Rel(I))*(sbess(I,lpr+2))*(Pnl(I)), kind=8)
                sumev(1)=sumev(1)+(sumd)
          endif
    enddo
      
    ! sum over odd terms

    do I=nrad-1,nstpod,-2
    
          if(lorb.eq.0) then
                sumd=real((Rel(I))*(sbess(I,lpr))*(Pnl(I)), kind=8)
                sumod(0)=sumod(0)+(sumd)
          endif
    
          if(lorb.eq.1) then
                if(lpr.ne.0) then
                      sumd=real((Rel(I))*(sbess(I,lpr-1))*(Pnl(I)), kind=8)
                else
                      sumd=0.0d0
                endif
                sumod(-1)=sumod(-1)+(sumd)
                sumd=real((Rel(I))*(sbess(I,lpr+1))*(Pnl(I)), kind=8)
                sumod(1)=sumod(1)+(sumd)
          endif
    
          if(lorb.eq.2) then
                if(lpr.gt.1) then
                      sumd=real((Rel(I))*(sbess(I,lpr-2))*(Pnl(I)), kind=8)
                else
                      sumd=0.0d0
                endif
                sumod(-1)=sumod(-1)+(sumd)
                sumd=real((Rel(I))*(sbess(I,lpr))*(Pnl(I)), kind=8)
                sumod(0)=sumod(0)+(sumd)
                sumd=real((Rel(I))*(sbess(I,lpr+2))*(Pnl(I)), kind=8)
                sumod(1)=sumod(1)+(sumd)
          endif
    enddo
    
    ! calculate integral
    if(lorb.eq.0) then
        gint_local(lpr,iqstp, 0)=h*real((4.0d0*sumod( 0)+2.0d0*sumev( 0))/3.0d0,kind=fp_kind)
    endif
    
    if(lorb.eq.1) then
        gint_local(lpr,iqstp,-1)=h*real((4.0d0*sumod(-1)+2.0d0*sumev(-1))/3.0d0,kind=fp_kind)
        gint_local(lpr,iqstp, 1)=h*real((4.0d0*sumod( 1)+2.0d0*sumev( 1))/3.0d0,kind=fp_kind)
    endif
    
    if(lorb.eq.2) then
        gint_local(lpr,iqstp,-1)=h*real((4.0d0*sumod(-1)+2.0d0*sumev(-1))/3.0d0,kind=fp_kind)
        gint_local(lpr,iqstp, 0)=h*real((4.0d0*sumod( 0)+2.0d0*sumev( 0))/3.0d0,kind=fp_kind)
        gint_local(lpr,iqstp, 1)=h*real((4.0d0*sumod( 1)+2.0d0*sumev( 1))/3.0d0,kind=fp_kind)
    endif
301 continue

    return
    end subroutine
    
    
    !--------------------------------------------------------------------------------
    ! This subroutine calculates continuum wavefunctions on a 0.0015 a.u. grid 
    ! to a radius of  15 a.u.  Normalisation is carried out by matching to Coulomb 
    ! functions at this radius.
    !--------------------------------------------------------------------------------
    !subroutine wavfunc_solo(E,lpr,Zatom,VA,nrad,Rel)
    subroutine wavfunc_solo(E,lpr,Zatom,Rel)
    implicit none
    integer(4) :: lpr,nl    !angular momentum
    integer(4) :: I,J       !step number
    integer(4) :: NSTEP     !number of steps
    integer(4) :: nstepx,icount

    ! output
    real(fp_kind) :: Pel(0:10000)    !continuum wave function
    real(fp_kind) :: Rel(0:2000)  !output wavefucntion
    ! input
    real(fp_kind) :: Zatom        !atomic number (1 for non hydrogenic potential)
    real(fp_kind) :: E            !energy of the free elecron  (Ryd.)
    ! local variables
    real(fp_kind) :: R,RP,RP2     !radial grid variables
    real(fp_kind) :: V,VP         !potential terms
    real(fp_kind) :: k,kP,kM,C,CP,CM,CON  !variables for Numerov
    real(fp_kind) :: kappa        !wavevector of ejected electron
    real(fp_kind) :: h,h2         !step size and (step size)*2
    real(fp_kind) :: alpha2       !fine structure constant squared
    real(fp_kind) :: DPel         !derivative of Pel
    real(fp_kind) :: f,fdot,g,gdot,Bkl,Ckl,delta,norm,Pout!normalisation variables
    real(fp_kind) :: testmax,beta
    real(8) :: dmax_val
    real(4) :: smax_val

    ! set initial values
    h=0.0015_fp_kind
    h2=2.0_fp_kind*h
    alpha2=1.0_fp_kind/(137.036_fp_kind)**2
    Pel(0)=0.0_fp_kind
    nl=lpr+1
    if(nl.gt.80) nl=80
    Pel(1)=h**nl
    CON=h**2/12.0_fp_kind
    NSTEP=10000
    kappa=sqrt(E)
    
    ! calculate initial step Pel(2)
30  R=h         !radius at I=1
    RP=R+h      !radius at I=2
    RP2=RP+h    !radius at I=3
    !calculate integration variables 
    V=VA(1)    !potential at I=1
    VP=VA(2)  !potential at I=2

    k=-(float(lpr*(lpr+1))/R**2 + V-E)
    kP=-(float(lpr*(lpr+1))/RP**2 + VP-E)

    !calculate Pel(2)
    Pel(2)=2.0_fp_kind*Pel(1)*(1.0_fp_kind-5.0_fp_kind*CON*k)/(1.0_fp_kind+kP*CON)

    !Start integration
    do I=2,NSTEP-1   !calculation of Pel(I+1)
          R=RP
          RP=RP2
          RP2=h*float(I+2)
          V=VA(I)    !potential at I
          VP=VA(I+1)  !potential at I+1

          kM=k
          k=kP
          kP=-(float(lpr*(lpr+1))/RP**2 + VP-E)

          ! Numerov Code
          CM=1.0_fp_kind+kM*CON
          C=1.0_fp_kind-5.0_fp_kind*k*CON
          CP=1.0_fp_kind+kP*CON
          Pel(I+1)=(2.0_fp_kind*Pel(I)*C-Pel(I-1)*CM)/CP
          !If Pel is too large reduce to avoid overflow
          !adding clauses for single and double precision
          if(fp_kind.eq.Double) dmax_val = 1.0d280
          if(fp_kind.eq.Single) smax_val = 1.0e37
          if(fp_kind.eq.Double.and.abs(Pel(I+1)).gt.dmax_val.or.fp_kind.eq.Single.and.abs(Pel(I+1)).gt.smax_val) then
              do J=1,I+1
                  
                  !if (Pel(J).le.1.0e-255_fp_kind) Pel(J)=0.0_fp_kind
                  !Pel(J)=Pel(J)/1.0e50_fp_kind
                  if (Pel(J).le.1.0e-37_fp_kind) Pel(J)=0.0_fp_kind
                  Pel(J)=Pel(J)/1.0e20_fp_kind
                  
              enddo
          endif
    enddo

    !Determine maximum value of Pel(NSTEPX-2) for best normalisation
    testmax=0.0_fp_kind
    do i = NSTEP-2,NSTEP-7,-1
          testmax=max(testmax,abs(Pel(I)))
          if(testmax.eq.abs(Pel(I))) NSTEPX=i+2
    enddo
    !Calculate derivative of Pel at NSTEPX-2
    DPel=(2.0_fp_kind*Pel(NSTEPX-4)-16.0_fp_kind*Pel(NSTEPX-3)+16.0_fp_kind*Pel(NSTEPX-1)-2.0_fp_kind*Pel(NSTEPX))/(24.0_fp_kind*h)

    ! normalise Rel
    if(Pel(NSTEP-2).eq.0.0_fp_kind) goto 201
    beta=DPel/Pel(NSTEPX-2)      !calculate logarithmic derivative
    !Calculate F,F',G,G'
    R=float(NSTEPX-2)*h

    call FJN(lpr,R*kappa,-Zatom/kappa,F,FDOT,G,GDOT)

    delta=ATAN(-(kappa*fdot-beta*f)/(kappa*gdot-beta*g))

    Bkl=2.0_fp_kind*cos(delta)
    Ckl=2.0_fp_kind*sin(delta)
    Pout=Bkl*f+Ckl*g
    norm=Pout/Pel(NSTEPX-2)
    icount=0
    do I=nrad,1,-1
          Rel(I)=real(norm*Pel(I))
    enddo
201 continue
    return
    end subroutine
    
    !--------------------------------------------------------------------------------
    ! This subroutine calculates spherical bessel funtions of orders 0 to 20 for
    ! a range of arguments from 0.0001 to 1200. For compatability with ABSION2 a
    ! radial grid of up to 2000 points is assumed and a step size of 0.0015 a.u.
    ! is used. Based on Abramowitz + Stegun p452 and recurrence relation 10.1.19
    ! due to J.C.P. Miller (Abramowitz + Stegun ref [9.20]).
    ! orders 21-50 are also calculated to a lower accuracy.
    !
    ! A.J.D 2014    There is a problem with sbtmp(4) at single precision, I have coded 
    !               the reccurance seeds manually using double precision and then convert
    !               to single precision.
    !--------------------------------------------------------------------------------
    subroutine sbessel(q,lprmax,sbess)
    implicit none
    
    integer(4) :: maxit,numwve,nradmax
    real(fp_kind) :: h
    
    parameter(h=0.0015_fp_kind,nradmax=2000,maxit=1200,numwve=100)
    ! external inputs
	real(fp_kind) :: q        !momentum transfer (a.u.)**(-1)
    integer(4) :: lprmax
    ! local variables
    real(fp_kind) :: F(0:maxit)       !dummy function
    real(8) :: sbtmp(0:4)       !temporary storage for analytic functions
    real(8) :: x                !argument
    real(fp_kind) :: bessmax          !maximum of first five orders
    real(fp_kind) :: norm             !normalisation constant
    integer(4) :: i,j,jj        !counters
    integer(4) :: Imax          !order of bessmax
    integer(4) :: nit           !number of iterations
    ! return values passed through common statement
    real(fp_kind) :: sbess(2000,0:numwve) !spherical bessel functions
    
    
    ! check to see if nrad is exceeded
    if(nrad.gt.nradmax) then
          write(*,21) nradmax
21        format(' nrad is larger than maximum allowed!',/,'  nrad will be reset to ',I5)
          pause 'press enter to continue or Ctrl+C to exit'
          nrad=nradmax
    endif
    
    ! calculate first five spherical bessel functions analytically to
    ! determine normalisation
    do i=nrad,1,-1
          x=dfloat(i)*h*q
          sbtmp(0)=dsin(x)/x
          sbtmp(1)=dsin(x)/(x*x)-dcos(x)/(x)
          sbtmp(2)=(3.0d0/(x*X*x)-1.0d0/x)*dsin(x)-3.0d0*dcos(x)/(x*x)
          sbtmp(3)=(15.0d0/(x*X*x*x)-6.0d0/(x*x))*dsin(x)+(1.0d0/x-15.0d0/(x*x*x))*dcos(x)
          sbtmp(4)=(1.05d2/(x*x*x*x*x)-45.0d0/(x*x*x)+1.0d0/x)*dsin(x)+(10.0d0/(x*x)-1.05d2/(x*x*x*x))*dcos(x)
          ! fill first five spherical bessel functions
          sbess(i,0)=sngl(sbtmp(0))
          sbess(i,1)=sngl(sbtmp(1))
          sbess(i,2)=sngl(sbtmp(2))
          sbess(i,3)=sngl(sbtmp(3))
          sbess(i,4)=sngl(sbtmp(4))
          ! select appropriate number of iterrations
          nit=int(x*1.1_fp_kind)
          if(nit.lt.lprmax*2) nit=lprmax*2
          if(nit.gt.maxit) nit=maxit
          if(nit.lt.60) nit=60
          F(nit)=0.0_fp_kind
          F(nit-1)=1.0_fp_kind
          ! calculate descending series for dummy variable F      
          do j=nit-1,1,-1
                F(j-1)=float(2*j+1)*F(j)/x -F(j+1) !Eqn 10.1.19
                !! ensure overflows and underflows do not occur
                !if(F(j-1).gt.1.0e250_fp_kind) then
                !      do jj=1199,j-1,-1
                !            if(F(jj).gt.1.0e-250_fp_kind) then
                !                  F(jj)=F(jj)/1.0e50_fp_kind
                !            else
                !                  F(jj)=0.0_fp_kind
                !            endif
                !      enddo
                !endif
                !! ensure overflows and underflows do not occur
                if(F(j-1).gt.1.0e25_fp_kind) then
                      do jj=1199,j-1,-1
                            if(F(jj).gt.1.0e-15_fp_kind) then
                                  F(jj)=F(jj)/1.0e20_fp_kind
                            else
                                  F(jj)=0.0_fp_kind
                            endif
                      enddo
                endif
                
                
          enddo  !end loop over order

          bessmax=0.0_fp_kind
          ! determine maximum value for best normalisation            
          do j=0,4
                bessmax=max(bessmax,abs(F(j)))
                if (bessmax.eq.abs(F(j))) Imax=j
          enddo

          norm=sbtmp(Imax)/F(Imax)
          ! normalise spherical bessel functions
          do j=lprmax-2,5,-1
                if(abs(F(j)*norm).lt.1.0e-37_fp_kind) then
                      sbess(i,j)=0.0_fp_kind
                else
                      sbess(i,j)=F(j)*norm
                endif
          enddo
    enddo   !end loop over radial grid

    return
    end subroutine
    
    
    
    
    
      
    !--------------------------------------------------------------------------------
    !gspline and gsplint are repeated functions that are ALMOST redundant
    !The function calls are slightly different to the spline and splint subroutiunes
    !in numerical_recipies.f90
    !I intend to replace the appropriate calls so to remove redundancy
    !--------------------------------------------------------------------------------
    

    SUBROUTINE gspline(x,n,lpr,ilam,gint_local,dgint_local)
    implicit none
    INTEGER(4) :: n,NMAX
    REAL(fp_kind) :: yp1,ypn,x(n)
    PARAMETER (NMAX=500)
    integer(4) :: numwve,nsplpts
    parameter(numwve=100,nsplpts=300)
    INTEGER(4) :: i,k,lpr,ilam
    REAL(fp_kind) :: p,qn,sig,un,u(NMAX)
    real(fp_kind) :: gint_local(0:numwve,nsplpts,-1:1),dgint_local(0:numwve,nsplpts,-1:1)


    yp1=1.0e31_fp_kind                             
    ypn=1.0e31_fp_kind

    if (yp1.gt.0.99e30_fp_kind) then
        dgint_local(lpr,1,ilam)=0.0_fp_kind
        u(1)=0.0_fp_kind
    else
        dgint_local(lpr,1,ilam)=-0.5_fp_kind
        u(1)=(3.0_fp_kind/(x(2)-x(1)))*((gint_local(lpr,2,ilam)-gint_local(lpr,1,ilam)) &
                & /(x(2)-x(1))-yp1)
    endif
    do i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*dgint_local(lpr,i-1,ilam)+2.0_fp_kind
        dgint_local(lpr,i,ilam)=(sig-1.0_fp_kind)/p
        u(i)=(6.0_fp_kind*((gint_local(lpr,i+1,ilam)-gint_local(lpr,i,ilam)) &
                & /(x(i+1)-x(i))-(gint_local(lpr,i,ilam)-gint_local(lpr,i-1,ilam)) &
                & /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
    enddo
    if (ypn.gt..99e30_fp_kind) then
        qn=0.0_fp_kind
        un=0.0_fp_kind
    else
        qn=0.5_fp_kind
        un=(3.0_fp_kind/(x(n)-x(n-1)))*(ypn-(gint_local(lpr,n,ilam)-gint_local(lpr,n-1,ilam)) &
                & /(x(n)-x(n-1)))
    endif
    dgint_local(lpr,n,ilam)=(un-qn*u(n-1))/(qn*dgint_local(lpr,n-1,ilam)+1.)
    do k=n-1,1,-1
        dgint_local(lpr,k,ilam)=dgint_local(lpr,k,ilam)*dgint_local(lpr,k+1,ilam)+u(k)
    enddo
    continue
    return
    END subroutine

    SUBROUTINE gsplint(xa,n,x,y,lpr,ilam,gint_local,dgint_local)
    implicit none
    INTEGER(4) :: n
    REAL(fp_kind) :: x,y,xa(n)
    INTEGER(4) :: k,khi,klo,lpr,ilam
    REAL(fp_kind) :: a,b,h
    integer(4) :: numwve,nsplpts
    parameter(numwve=100,nsplpts=300)
    real(fp_kind) gint_local(0:numwve,nsplpts,-1:1),dgint_local(0:numwve,nsplpts,-1:1)                                    

    klo=1                  
    khi=n
1   if (khi-klo.gt.1) then
          k=(khi+klo)/2
          if(xa(k).gt.x)then
          khi=k
    else
          klo=k
    endif
    goto 1
    endif
    h=xa(khi)-xa(klo)
    if (h.eq.0.0_fp_kind) pause 'bad xa input in splint'
    a=(xa(khi)-x)/h
    b=(x-xa(klo))/h
    y=a*gint_local(lpr,klo,ilam)+b*gint_local(lpr,khi,ilam)+((a*a*a-a)*dgint_local(lpr,klo,ilam)+(b*b*b-b)*dgint_local(lpr,khi,ilam))*(h*h)/6.
    return
    END subroutine
!  (C) Copr. 1986-92 Numerical Recipes Software ]2+18Z9.


!--------------------------------------------------------------------------------


    SUBROUTINE FJN(L_val,XX,ETA1,F,FDOT,G,GDOT)
    !-------------------------------------------------
    ! RETURNS EITHER BESSEL/NEUMANN OR COULOMB FUCTIONS
    !-------------------------------------------------
    implicit none
    
    !IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    real(fp_kind) :: xx,xlmin,xlmax,paccq,eta1,g,gdot,f,fdot
    integer(4) :: npq,nfp,iexp,m1,mode1,ifail,kfn
    real(fp_kind),dimension(500) :: FC,GC,FCP,GCP
    
    
    integer(4) :: l_val
    COMMON /STEED/ PACCQ,NFP,NPQ,IEXP,M1

    XLMIN=float(L_val)
    XLMAX=float(L_val)
    MODE1=1
    KFN=0
    IFAIL=0
    CALL COULFG(XX,ETA1,XLMIN,XLMAX,FC,GC,FCP,GCP,MODE1,KFN,IFAIL)
    F=FC(M1)
    FDOT=FCP(M1)
    G=GC(M1)
    GDOT=GCP(M1)
    RETURN
    END subroutine

!--------------------------------------------------------------------------------

    SUBROUTINE COULFG(XX,ETA1,XLMIN,XLMAX, FC,GC,FCP,GCP,MODE1,KFN,IFAIL)
      
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C                                                                      C
!C  REVISED COULOMB WAVEFUNCTION PROGRAM USING STEED'S METHOD           C
!C                                                                      C
!C  A. R. BARNETT           MANCHESTER  MARCH   1981                    C
!C                                                                      C
!C  ORIGINAL PROGRAM 'RCWFN'      IN    CPC  8 (1974) 377-395           C
!C                 + 'RCWFF'      IN    CPC 11 (1976) 141-142           C
!C  FULL DESCRIPTION OF ALGORITHM IN    CPC 21 (1981) 297-314           C
!C  THIS VERSION WRITTEN UP       IN    CPC XX (1982) YYY-ZZZ           C
!C                                                                      C
!C  COULFG RETURNS F,G,F',G', FOR REAL XX.GT.0,REAL ETA1 (INCLUDING 0), C
!C   AND REAL LAMBDA(XLMIN) .GT. -1 FOR INTEGER-SPACED LAMBDA VALUES    C
!C   THUS GIVING POSITIVE-ENERGY SOLUTIONS TO THE COULOMB SCHRODINGER   C
!C   EQUATION,TO THE KLEIN-GORDON EQUATION AND TO SUITABLE FORMS OF     C
!C   THE DIRAC EQUATION ,ALSO SPHERICAL & CYLINDRICAL BESSEL EQUATIONS  C
!C                                                                      C
!C  FOR A RANGE OF LAMBDA VALUES (XLMAX - XLMIN) MUST BE AN INTEGER,    C
!C  STARTING ARRAY ELEMENT IS M1 = MAX (IDINT(XLMIN+ACCUR),0) + 1       C
!C      SEE TEXT FOR MODIFICATIONS FOR INTEGER L-VALUES                 C
!C                                                                      C
!C  IF 'MODE' = 1  GET F,G,F',G'   FOR INTEGER-SPACED LAMBDA VALUES     C
!C            = 2      F,G      UNUSED ARRAYS MUST BE DIMENSIONED IN    C
!C            = 3      F               CALL TO AT LEAST LENGTH (1)      C
!C  IF 'KFN'  = 0 REAL        COULOMB FUNCTIONS ARE RETURNED            C
!C            = 1 SPHERICAL   BESSEL      '      '     '                C
!C            = 2 CYLINDRICAL BESSEL      '      '     '                C
!C  THE USE OF 'MODE' AND 'KFN' IS INDEPENDENT                          C
!C                                                                      C
!C  PRECISION:  RESULTS TO WITHIN 2-3 DECIMALS OF 'MACHINE ACCURACY'    C
!C   IN OSCILLATING REGION X .GE. ETA1 + SQRT(ETA1**2 + XLM(XLM+1))     C
!C   COULFG IS CODED FOR REAL*8 ON IBM OR EQUIVALENT  ACCUR = 10**-16   C
!C   USE AUTODBL + EXTENDED PRECISION ON HX COMPILER  ACCUR = 10**-33   C
!C   FOR MANTISSAS OF 56 & 112 BITS. FOR SINGLE PRECISION CDC (48 BITS) C
!C   REASSIGN DSQRT=SQRT ETC.  SEE TEXT FOR COMPLEX ARITHMETIC VERSION  C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    implicit none
    !IMPLICIT REAL*8 (A-H,O-Z)
    real(fp_kind) :: xx,xlmin,xlmax,rt2dpi,paccq,xll,w,tm30,sign
    real(fp_kind) :: xl,wi,tk,ta,sl,rl,q,el,dr,dq,dp,di,c,br,bi,ar,ai
    real(fp_kind) :: xlm,x,gjwkb,fjwkb,xi,px,pk,fcl,eta,eta1
    real(fp_kind) :: zero,one,two,ten2,abort,tk1,pk1,p
    real(fp_kind) :: fpl,f,ek,e2mm1,df,dell,d,accur,acch,acc4,acc
    real(fp_kind) :: half,gpl,gcl1,gcl,gam,fcm,fcl1,beta,b,alpha,a
    integer(4) :: nfp,nfq,iexp,m1,mode,mode1,kfn,ifail,npq,lxtra,l1
    integer(4) :: maxl,lp,l
    
    
    
    
    real(fp_kind),dimension(500) :: FC,GC,FCP,GCP
    LOGICAL      ETANE0,XLTURN
    COMMON       /STEED/ PACCQ,NFP,NPQ,IEXP,M1
    !***  COMMON BLOCK IS FOR INFORMATION ONLY.  NOT REQUIRED IN CODE
    !***  COULFG HAS CALLS TO: DSQRT,DABS,DMOD,IDINT,DSIGN,DREAL,DMIN1
    DATA ZERO,ONE,TWO,TEN2,ABORT /0.0_fp_kind, 1.0_fp_kind, 2.0_fp_kind, 1.0_fp_kind, 2.0e4_fp_kind/
    DATA HALF,TM30 / 0.5_fp_kind, 1.0e-30_fp_kind /
      !DATA RT2DPI /0.79788456080286535587989211986876373_fp_kind/
    ! *** THIS CONSTANT IS  DSQRT(TWO/PI):  USE Q0 FOR IBM REAL*16: D0 FOR
    ! *** REAL*8 & CDC DOUBLE P:  E0 FOR CDC SINGLE P; AND TRUNCATE VALUE.
    rt2dpi = sqrt(2.0_fp_kind/pi)
    if(fp_kind.eq.Double) then
        ACCUR = 1.0d-16
    elseif(fp_kind.eq.Single) then
        ACCUR = 1.0e-16_fp_kind
    endif
    ! ***            CHANGE ACCUR TO SUIT MACHINE AND PRECISION REQUIRED
    MODE  = 1
    IF(MODE1 .EQ. 2 .OR. MODE1 .EQ. 3 ) MODE = MODE1
    IFAIL = 0
    IEXP  = 1
    NPQ   = 0
    ETA   = ETA1
    GJWKB = ZERO
    PACCQ = ONE
    IF(KFN .NE. 0) ETA = ZERO
    ETANE0  = ETA .NE. ZERO
    ACC   = ACCUR
    ACC4  = ACC*TEN2*TEN2
    ACCH  =  SQRT(ACC)
    ! ***    TEST RANGE OF XX, EXIT IF.LE.DSQRT(ACCUR) OR IF NEGATIVE
    !
    IF(XX .LE. ACCH)                          GO TO 100
    X     = XX
    XLM   = XLMIN
    IF(KFN .EQ. 2)  XLM = XLM - HALF
    IF(XLM .LE. -ONE .OR. XLMAX .LT. XLMIN)   GO TO 105
    E2MM1 = ETA*ETA + XLM*XLM + XLM
    XLTURN= X*(X - TWO*ETA) .LT. XLM*XLM + XLM
    DELL  = XLMAX - XLMIN + ACC
    IF( ABS(MOD(DELL,ONE)) .GT. ACC) continue !write(*,2040)XLMAX,XLMIN,DELL    !commented out dangerous error
    LXTRA =   INT(DELL)
    XLL   = XLM +  float(LXTRA)
    ! ***       LXTRA IS NUMBER OF ADDITIONAL LAMBDA VALUES TO BE COMPUTED
    ! ***       XLL  IS MAX LAMBDA VALUE, OR 0.5 SMALLER FOR J,Y BESSELS
    ! ***         DETERMINE STARTING ARRAY ELEMENT (M1) FROM XLMIN
    M1  = MAX (  INT(XLMIN + ACC),0) + 1
    L1  = M1 + LXTRA
    !
    ! ***    EVALUATE CF1  =  F   =  FPRIME(XL,ETA,X)/F(XL,ETA,X)
    !
    XI  = ONE/X
    FCL = ONE
    PK  = XLL + ONE
    PX  = PK  + ABORT
2   EK  = ETA / PK
    F   = (EK + PK*XI)*FCL + (FCL - ONE)*XI
    PK1 =  PK + ONE
    ! ***   TEST ENSURES B1 .NE. ZERO FOR NEGATIVE ETA; FIXUP IS EXACT.
    IF( ABS(ETA*X + PK*PK1) .GT. ACC)  GO TO 3
    FCL  = (ONE + EK*EK)/(ONE + (ETA/PK1)**2)
    PK   =  TWO + PK
    GO TO 2
3   D   =  ONE/((PK + PK1)*(XI + EK/PK1))
    DF  = -FCL*(ONE + EK*EK)*D
    IF(FCL .NE. ONE )  FCL = -ONE
    IF(D   .LT. ZERO)  FCL = -FCL
    F   =  F  + DF
    !
    ! ***   BEGIN CF1 LOOP ON PK = K = LAMBDA + 1
    !
    P     = ONE
4   PK    = PK1
    PK1 = PK1 + ONE
    EK  = ETA / PK
    TK  = (PK + PK1)*(XI + EK/PK1)
    D   =  TK - D*(ONE + EK*EK)
    IF( ABS(D) .GT. ACCH)             GO TO 5
    WRITE (6,1000) D,DF,ACCH,PK,EK,ETA,X
    P = P  +   ONE
    IF( P .GT. TWO )                  GO TO 110
5   D = ONE/D
    IF (D .LT. ZERO) FCL = -FCL
    DF  = DF*(D*TK - ONE)
    F   = F  + DF
    IF(PK .GT. PX)                    GO TO 110
    IF( ABS(DF) .GE.  ABS(F)*ACC)             GO TO 4
    NFP = PK - XLL - 1
    IF(LXTRA .EQ. 0)                          GO TO 7
    !     
    ! *** DOWNWARD RECURRENCE TO LAMBDA = XLM. ARRAY GC,IF PRESENT,STORES RL
    !
    FCL = FCL*TM30
    FPL = FCL*F
    IF(MODE .EQ. 1) FCP(L1) = FPL
    FC (L1) = FCL
    XL  = XLL
    RL  = ONE
    EL  = ZERO
    DO  LP = 1,LXTRA
        IF(ETANE0) EL = ETA/XL
        IF(ETANE0) RL =  SQRT(ONE + EL*EL)
        SL    =  EL  + XL*XI
        L     =  L1  - LP
        FCL1  = (FCL *SL + FPL)/RL
        FPL   =  FCL1*SL - FCL *RL
        FCL   =  FCL1
        FC(L) =  FCL
        IF(MODE .EQ. 1) FCP(L)  = FPL
        IF(MODE .NE. 3 .AND. ETANE0) GC(L+1) = RL
    ENDDO
6   XL = XL - ONE
    IF(FCL .EQ. ZERO) FCL = ACC
    F  = FPL/FCL
    ! ***    NOW WE HAVE REACHED LAMBDA = XLMIN = XLM
    ! ***    EVALUATE CF2 = P + I.Q  AGAIN USING STEED'S ALGORITHM
    ! ***    SEE TEXT FOR COMPACT COMPLEX CODE FOR SP CDC OR NON-ANSI IBM
7   IF( XLTURN ) CALL JWKB(X,ETA,MAX(XLM,ZERO),FJWKB,GJWKB,IEXP)
    IF( IEXP .GT. 1 .OR. GJWKB .GT. ONE/(ACCH*TEN2))  GO TO 9
    XLTURN = .FALSE.
    TA =  TWO*ABORT
    PK =  ZERO
    WI =  ETA + ETA
    P  =  ZERO
    Q  =  ONE - ETA*XI
    AR = -E2MM1
    AI =  ETA
    BR =  TWO*(X - ETA)
    BI =  TWO
    DR =  BR/(BR*BR + BI*BI)
    DI = -BI/(BR*BR + BI*BI)
    DP = -XI*(AR*DI + AI*DR)
    DQ =  XI*(AR*DR - AI*DI)
8   P  = P  + DP
    Q  = Q  + DQ
    PK = PK + TWO
    AR = AR + PK
    AI = AI + WI
    BI = BI + TWO
    D  = AR*DR - AI*DI + BR
    DI = AI*DR + AR*DI + BI
    C  = ONE/(D*D + DI*DI)
    DR =  C*D
    DI = -C*DI
    A  = BR*DR - BI*DI - ONE
    B  = BI*DR + BR*DI
    C  = DP*A  - DQ*B
    DQ = DP*B  + DQ*A
    DP = C
    IF(PK .GT. TA) GO TO 120
    IF( ABS(DP)+ ABS(DQ).GE.( ABS(P)+ ABS(Q))*ACC) GO TO 8
    NPQ   = PK/TWO
    PACCQ = HALF*ACC/MIN( ABS(Q),ONE)
    IF( ABS(P) .GT.  ABS(Q)) PACCQ = PACCQ* ABS(P)
    !
    ! *** SOLVE FOR FCM = F AT LAMBDA = XLM,THEN FIND NORM FACTOR W=W/FCM
    !
    GAM = (F - P)/Q
    IF(Q .LE. ACC4* ABS(P))             GO TO 130
    W   = ONE/ SQRT((F - P)*GAM + Q)
    GO TO 10
    ! *** ARRIVE HERE IF G(XLM) .GT. 10**6 OR IEXP .GT. 70 & XLTURN = .TRUE.
9   W   = FJWKB
    GAM = GJWKB*W
    P   = F
    Q   = ONE
    !
    ! *** NORMALISE FOR SPHERICAL OR CYLINDRICAL BESSEL FUNCTIONS
    !
10  ALPHA = ZERO
    IF(KFN  .EQ. 1) ALPHA = XI
    IF(KFN  .EQ. 2) ALPHA = XI*HALF
          BETA  = ONE
    IF(KFN  .EQ. 1) BETA  = XI
    IF(KFN  .EQ. 2) BETA  =  SQRT(XI)*RT2DPI
    FCM  =  SIGN(W,FCL)*BETA
    FC(M1)  = FCM
    IF(MODE .EQ. 3)           GO TO 11
    IF(.NOT. XLTURN)   GCL =  FCM*GAM
    IF(      XLTURN)   GCL =  GJWKB*BETA
    IF( KFN .NE. 0 )   GCL = -GCL
    GC(M1)  = GCL
    GPL =  GCL*(P - Q/GAM) - ALPHA*GCL
    IF(MODE .EQ. 2)           GO TO 11
    GCP(M1) = GPL
    FCP(M1) = FCM*(F - ALPHA)
11  IF(LXTRA .EQ. 0 ) RETURN
    ! *** UPWARD RECURRENCE FROM GC(M1),GCP(M1)  STORED VALUE IS RL
    ! *** RENORMALISE FC,FCP AT EACH LAMBDA AND CORRECT REGULAR DERIVATIVE
    ! ***    XL   = XLM HERE  AND RL = ONE , EL = ZERO FOR BESSELS
    W    = BETA*W/ ABS(FCL)
    MAXL = L1 - 1
    DO 12 L = M1,MAXL
          IF(MODE .EQ. 3)           GO TO 12
          XL = XL + ONE
          IF(ETANE0)   EL = ETA/XL
          IF(ETANE0)   RL = GC(L+1)
          SL = EL + XL*XI
          GCL1     = ((SL - ALPHA)*GCL - GPL)/RL
          GPL      =   RL*GCL -  (SL + ALPHA)*GCL1
          GCL      = GCL1
          GC(L+1)  = GCL1
          IF(MODE .EQ. 2)           GO TO 12
          GCP(L+1) = GPL
          FCP(L+1) = W*(FCP(L+1) - ALPHA*FC(L+1))
12  FC(L+1)     = W* FC(L+1)
    RETURN
1000 FORMAT(/' CF1 ACCURACY LOSS: D,DF,ACCH,K,ETA/K,ETA,X = ',1P7D9.2/)
!
! ***    ERROR MESSAGES
!
100 IFAIL = -1
      write(*,2000) XX,ACCH
2000  FORMAT(' FOR XX = ',1PD12.3,' TRY SMALL-X  SOLUTIONS',' OR X NEGATIVE'/ ,' SQUARE ROOT ACCURACY PARAMETER =  ',D12.3/)
      RETURN
105 IFAIL = -2
      WRITE (6,2005) XLMAX,XLMIN,XLM
2005  FORMAT(/' PROBLEM WITH INPUT ORDER VALUES:XLMAX,XLMIN,XLM = ',1P3D15.6/)
      RETURN
110 IFAIL =  1
      WRITE (6,2010) ABORT,F ,DF,PK,PX,ACC
2010  FORMAT(' CF1 HAS FAILED TO CONVERGE AFTER ',F10.0,' ITERATIONS',/' F,DF,PK,PX,ACCUR =  ',1P5D12.3//)
      RETURN
!120   iswitch=1
120    continue  
      return
!  120 IFAIL =  2
!	WRITE (6,2020) ABORT,P,Q,DP,DQ,ACC
2020  FORMAT(' CF2 HAS FAILED TO CONVERGE AFTER ',F7.0,' ITERATIONS',/' P,Q,DP,DQ,ACCUR =  ',1P4D17.7,D12.3//)
      RETURN
!130   iswitch=1
130   continue   
      return
!  130 IFAIL =  3
!	WRITE (8,2030) P,Q,ACC,DELL,LXTRA,M1
2030  FORMAT(' FINAL Q.LE.DABS(P)*ACC*10**4 , P,Q,ACC = ',1P3D12.3,4X,' DELL,LXTRA,M1 = ',D12.3,2I5 /)
      RETURN
2040  FORMAT(' XLMAX - XLMIN = DELL NOT AN INTEGER ',1P3D20.10/)
      END subroutine

    !--------------------------------------------------------------------------------
    SUBROUTINE JWKB(XX,ETA1,XL,FJWKB,GJWKB,IEXP) !may not be necessary AJD 2010
    implicit none
    
    !IMPLICIT REAL*8 (A-H,O-Z)
    integer(4) :: iexp
	REAL(fp_kind) :: XX,ETA1,XL,FJWKB,GJWKB
    REAL(fp_kind) :: DZERO
    REAL(fp_kind) :: zero,xll1,x,sl,rl2,phi10,phi,one,hll,hl,gh2,gh,eta
    REAL(fp_kind) :: ten,six,rl35,half,aloge
      
    ! *** COMPUTES JWKB APPROXIMATIONS TO COULOMB FUNCTIONS    FOR XL.GE. 0
    ! *** AS MODIFIED BY BIEDENHARN ET AL. PHYS REV 97 (1955) 542-554
    ! *** CALLS DMAX1,SQRT,LOG,EXP,ATAN2,REAL,INT         BARNETT FEB 1981
    DATA   ZERO,HALF,ONE,SIX,TEN/ 0.0E0_fp_kind, 0.5E0_fp_kind, 1.0E0_fp_kind, 6.0E0_fp_kind, 10.0E0_fp_kind /
    DATA  DZERO, RL35, ALOGE  /0.0E0_fp_kind, 35.0E0_fp_kind, 0.4342944819E0_fp_kind /
    X     = XX
    ETA   = ETA1
    GH2   = X*(ETA + ETA - X)
    XLL1  = MAX(XL*XL + XL,DZERO)
    IF(GH2 + XLL1 .LE. ZERO) RETURN
    HLL  = XLL1 + SIX/RL35
    HL   = SQRT(HLL)
    SL   = ETA/HL + HL/X
    RL2  = ONE + ETA*ETA/HLL
    GH   = SQRT(GH2 + HLL)/X
    PHI  = X*GH - HALF*( HL*LOG((GH + SL)**2/RL2) - LOG(GH) )
        IF(ETA .NE. ZERO) PHI = PHI - ETA*ATAN2(X*GH,X - ETA)
    PHI10 = -PHI*ALOGE
    IEXP  =  INT(PHI10)
    IF(IEXP .GT. 70) GJWKB = TEN**(PHI10 - float(IEXP))
    IF(IEXP .LE. 70) GJWKB = EXP(-PHI)
    IF(IEXP .LE. 70) IEXP  = 0
    FJWKB = HALF/(GH*GJWKB)
    RETURN
    END subroutine
    
    !--------------------------------------------------------------------------------
    ! This routine sets up the filter for significant transitions based on their 
    ! intensities. The transitions with lowest intensities are filtered out.
    ! The cutoff is set by the relative input threshold. This threshold sets the
    ! fraction of the total intensity that is omitted.
    ! Threshold is expected to be in the range [0,1] and usually less than 1.
    subroutine strength_filter(strength, valid, n, threshold)
    
    use m_numerical_tools
    
    implicit none
    
    real(fp_kind), parameter :: eps = 1.0e-6_fp_kind ! Small filter threshold
    
    integer(4), intent(in) :: n
    integer(4), intent(out) :: valid(n)
    real(fp_kind), intent(in) :: strength(n), threshold
    
    integer(4) :: i
    integer(4) :: index_sorted(n) ! Sorted indices (ascending by strength)
    real(fp_kind) :: sum_strength, cutoff_strength, accum
    
    if (n<1) return ! no items, nothing to filter
    
    ! Initialize valid array
    valid = 1 ! initialize all to be valid
    if (threshold < eps) return ! no filtering if threshold is too small
    
    ! Sort indices based on intensity (values remain unchanged)
    call heapsort_indices_only(strength, index_sorted, n)
    
    ! Compute total intensity and cutoff value
    sum_strength = SUM(strength)
    cutoff_strength = sum_strength * threshold
    accum = 0.0_fp_kind
    
    ! Accumulate from the least significant intensity upward
    do i = 1, n
        accum = accum + strength(index_sorted(i))
        if (accum < cutoff_strength) then
            valid(index_sorted(i)) = 0 ! mark this index as not valid
        else
            exit ! stop when threshold is reached
        end if
    end do
        
    return
    end subroutine strength_filter
    
    
    
    subroutine fill_tmatrix(tmatrix,m_l,lpr,m_lpr,qpos_local,qmin_local,gint_local,dgint_local)
    use global_variables
    use m_crystallography, only: trimr, trimi
        
    implicit none
    
    !input variables
    integer(4) :: lpr,m_l,m_lpr
    integer(4) :: lammin,lammax
    real(fp_kind) :: gint_local(0:numwve,nsplpts,-1:1),dgint_local(0:numwve,nsplpts,-1:1)
    real(fp_kind) :: qmin_local,const
    real(fp_kind) :: qpos_local(300)
   
    !calculation variables
    integer(4) :: shifty,shiftx,m1,m2,ny,nx 
    integer(4) :: m_lam,ilam,ilamp
    real(fp_kind) :: threej,tgint
    complex(fp_kind) ::  spherical_harm !,sph_harm
     
    !calculation vectors
    real(fp_kind) :: qy(3),qx(3),qz(3)
    real(fp_kind) :: qx_min,qy_min
    real(fp_kind) :: qx_mag,qy_mag,q_perp_mag,qz_mag,q_vec_mag
    real(fp_kind) :: q_perp(3), q_vec(3)
    real(fp_kind) :: qtmp1, ifx, ify
    real(fp_kind) :: lpr_temp, ilam_temp, l_temp, m_lam_temp, m_l_temp
    !calculation summing variables
    complex(fp_kind) :: csum,cval,fq
    !calculation angular variables
    real(fp_kind) :: theta,phi
    !output variables
    complex(fp_kind) :: tmatrix(nopiy,nopix)
    
    ifx = 1.0_fp_kind / real(ifactorx,fp_kind)
    ify = 1.0_fp_kind / real(ifactory,fp_kind)
    lammax = lorb + lpr
    lammin = abs(lorb - lpr)
    shifty = (nopiy-mod(nopiy,2))/2
    shiftx = (nopix-mod(nopix,2))/2
    
    qz(1:3) = qmin_local * orthog(1:3,3) / trimr(orthog(1:3,3),ss)
    qz_mag = qmin_local
    qx_min = trimi(ig1,ss) *ifx
    qy_min = trimi(ig2,ss) *ify
    
    !$OMP PARALLEL DO &
    !$OMP& PRIVATE(m1, m2, qx, qx_mag, qy, qy_mag, q_perp, q_perp_mag, q_vec, q_vec_mag, &
    !$OMP&         csum, m_lam, theta, ilam, cval, ilamp, qtmp1, threej, phi, spherical_harm, &
    !$OMP&         fq, tgint, lpr_temp, ilam_temp, l_temp, m_lam_temp, m_l_temp) &
    !$OMP& SHARED(tmatrix)
    do ny=1,nopiy
        m2 = mod( ny-1+shifty, nopiy) - shifty
        qy(1:3) = real(m2 * ig2(1:3), fp_kind)
        qy_mag = trimr(qy,ss)
        qy(1:3) = qy(1:3) * ify
        do nx=1,nopix
            m1 = mod( nx-1+shiftx, nopix) - shiftx
            qx(1:3) = real(m1 * ig1(1:3), fp_kind)
            qx_mag = trimr(qx,ss)
                
            !q_perp(1:3) = qy(1:3) + qx(1:3) * ifx      !qx has been divided by ifactor in the outer loop
            ! BDF FIXED 2016/02/22: changed to ifactorx
            q_perp(1:3) = qy(1:3) + qx(1:3) * ifx      !qx has been divided by ifactor in the outer loop
                
            q_perp_mag = trimr(q_perp,ss)
            q_vec(1:3) = q_perp(1:3) + qz(1:3)
            q_vec_mag = trimr(q_vec,ss)
            !BWL addition 
            if (q_vec_mag.ge.bwl_rad) then
                tmatrix(ny,nx)= 0.0_fp_kind !czero
                cycle
            endif
            !BWL addition - end
            csum = 0.0_fp_kind !czero
            m_lam = m_l - m_lpr
            theta = pi/2.0_fp_kind + atan2(qmin_local,q_perp_mag)
            do ilam = lammin,lammax
                if (mod(ilam+lorb+lpr,2).ne.0) cycle
                if (abs(m_lam).gt.ilam) cycle
                cval = (-ci) ** real(ilam,fp_kind) * sqrt(real(2*ilam+1,fp_kind))
                if (ilam.eq.lpr) then
                    ilamp = 0
                elseif (ilam.eq.(lpr-lorb)) then
                    ilamp = -1
                elseif (ilam.eq.(lpr+lorb)) then
                    ilamp = 1
                endif
                qtmp1 = q_vec_mag * tpb
                !numerics to make function call cleb work
                lpr_temp = real(lpr, fp_kind)
                ilam_temp = real(ilam, fp_kind)
                l_temp = real(lorb, fp_kind)
                m_lam_temp = real(m_lam, fp_kind)
                m_l_temp = real(m_l, fp_kind)
                call gsplint(qpos_local,nsplpts,qtmp1,tgint,lpr,ilamp,gint_local,dgint_local)
                threej = cleb(lpr_temp,ilam_temp,l_temp,0.0_fp_kind,0.0_fp_kind)* &
                            & (-1.0_fp_kind)**(lpr-ilam)/ sqrt(real(2*lorb+1,fp_kind))
                cval = cval * threej * tgint
                threej = cleb(lpr_temp,ilam_temp,l_temp,-m_lam_temp,-m_l_temp)* &
                            & (-1.0_fp_kind)**(-m_l+lpr-ilam) / sqrt(real(2*lorb+1,fp_kind))
                phi = atan2(sign(qy_mag,real(m2,fp_kind)),sign(qx_mag,real(m1,fp_kind)))
                spherical_harm=sph_harm2(ilam,m_lam,theta,phi)
                cval = cval * threej * conjg(spherical_harm)
                csum = csum + cval
            enddo
            fq = ((-1)**m_lpr)*((-1)**m_lam)*sqrt( 2.0_fp_kind*tp* &
                            & real(2*lpr+1,fp_kind)*real(2*lorb+1,fp_kind) ) * csum
            tmatrix(ny,nx) = (fq / (q_vec_mag*q_vec_mag))
        enddo
    enddo
    !$OMP END PARALLEL DO
    !The tmatrix in k-space at 0,0 need to shift for atom positions
    !Potential now in units of Angstrom . sqrt(eV)
    tmatrix = tmatrix * sqrt(real(nopiy*nopix,fp_kind))
    const = sqrt(2.0_fp_kind)* qx_min * qy_min * fsc * hbarc / pi
	tmatrix = const*tmatrix     ! (the extra factor of dsqrt(2) accounts for spin)
    !call ifft2(nopiy,nopix,tmatrix,nopiy,tmatrix,nopiy)    !Fourier transform to get it back to real space
    return
    end subroutine
    
    
    
    !----------------------------------------------------------------------------------------------
    !this subroutine sets up the transition matrices for each state.
    !This assumes that the atom is at the origin.
    !After making the transition matrix elements, they are then shifted using the 1D factorisation
    !The arrays Hn0_shiftx_coord and Hn0_shifty_coord contain the factorisation.
    
    function setup_ms_hn0_tmatrices(nopiy,nopix,nstates) result(tmatrix_states)
    
        use output
        use m_multislice
    
        implicit none

		integer*4,intent(in)::nopiy,nopix,nstates
		complex(fp_kind)::tmatrix_states(nopiy,nopix,nstates)
    
        !constant (hn0)
        integer(4) :: i,j
        integer(4) :: ml,lpr,mlpr
        integer(4) :: atom,natoms
        integer(4) :: ss_atomx,ss_atomy
        integer(4) :: counter
        real(fp_kind) :: x_coord,y_coord,x_cell,y_cell,ifx,ify
        real(fp_kind) :: state_intens(nstates) ! intensity of each state
        !real(fp_kind),dimension(nopiy,nopix) :: effective_potential

        complex(fp_kind)   :: c_1 = 9.7846113e-07_fp_kind ! c_1 = 1/(2mc^2) in eV^(-1)

        !Set the interaction constant
        alpha_n= (ci*relm)/(4.0_fp_kind*c_1*(hbarc**2)*pi*ak)   !this should be k_n
        
        !make the transition matrix elements and the 1D shift arrays
        !effective_potential = 0.0_fp_kind

        write(*,*) 'Calculating transition potentials...'
        
        state_intens = 0.0_fp_kind
        do i =1,nstates
            ml=state_vector(i,1)
            lpr=state_vector(i,2)
            mlpr=state_vector(i,3)
            write(6,100,ADVANCE='NO') achar(13),i,nstates,lorb,ml,lpr,mlpr
            flush(6)
100         format(a1,'   Transition ',i0,'/',i0,':  [',i0,',',i0,"]  ->  [",i0,',',i0,']   ')
101         format(a1,'                                                           ')            
            call calc_gsplint(ml,lpr,mlpr,qpos(:,i),qmin(i),gint(:,:,:,i),dgint(:,:,:,i))
            call fill_tmatrix(tmatrix_states(:,:,i),ml,lpr,mlpr,qpos(:,i),qmin(i), &
                    & gint(:,:,:,i),dgint(:,:,:,i)) !make the Hn0
            state_intens(i) = sum(abs(tmatrix_states(:,:,i))**2) !calculate the total intensity of the Hn0
            !effective_potential = effective_potential + abs(tmatrix_states(:,:,i))**2
        enddo
        write(6,101,ADVANCE='NO') achar(13)
        flush(6)
        
        if (0 < IAND(arg_debug_dump, 1)) then
            write(*,*)
            write(*,*) '  Creating dump of transition potentials (y,x,tran) [dump_tmat.bin]'
            open(unit=456,file='dump_tmat.bin',form='binary',status='replace',convert='big_endian')
            write(456) tmatrix_states
            close(456)
        endif
        
        ! Hn0 strength filter
        if (Hn0_filter_threshold > 0.0_fp_kind .and. Hn0_filter_threshold < 1.0_fp_kind) then
          write(*,*)
          write(*,*) '  Filtering insignificant transition ...'
          if (ALLOCATED(state_use)) deallocate(state_use)
          allocate(state_use(nstates))
          state_use = 1 ! initialize all states to be used
          call strength_filter(state_intens, state_use, nstates, Hn0_filter_threshold)
          if (SUM(state_use) < nstates) then
110         format('   ',i0,' of ',i0,' transitions used, filtered with rel. threshold ',E9.2)
            write(*,110) SUM(state_use), nstates, Hn0_filter_threshold
            if (SUM(state_use) < nstates) then
111           format('   The following ',i0,' transitions are ignored (total strength: ',E9.2,')')
              write(*,111) nstates-SUM(state_use), SUM(state_intens)
              do i=1, nstates
112             format('   #',i3,': l=',i3,', ml=',i3,"  ->  l'=",i3,", ml'=",i3,':   strength=',E9.2)
                if (state_use(i) == 0) write(*,112) i, lorb, state_vector(i,1), state_vector(i,2), &
                                                             state_vector(i,3), state_intens(i)
              enddo
            endif
          endif
        endif
        
        write(*,*)
        write(*,*) '  Preparing shift factors ...'
        
        !number of atoms in a slice
		if(allocated(natoms_slice_total)) deallocate(natoms_slice_total)
        allocate(natoms_slice_total(n_slices))
        !natoms_slice_total = nat_slice_unitcell(target_atom,:)*ycells*xcells
        natoms_slice_total = nat_slice_unitcell(target_atom,:)*ifactory*ifactorx
    
        !shift arrays
		if(allocated(Hn0_shiftx_coord)) deallocate(Hn0_shiftx_coord)
		if(allocated(Hn0_shifty_coord)) deallocate(Hn0_shifty_coord)
		if(allocated(ion_tau)) deallocate(ion_tau)
        allocate(Hn0_shiftx_coord(nopix,maxval(natoms_slice_total),n_slices))
        allocate(Hn0_shifty_coord(nopiy,maxval(natoms_slice_total),n_slices))
		allocate(ion_tau(2,sum(natoms_slice_total)))
        ifx = 1.0_fp_kind / real(ifactorx, fp_kind)
        ify = 1.0_fp_kind / real(ifactory, fp_kind)
!102     format('(dbg)  ',a,'= [',F6.3,',',F6.3,']')
!        write(*,102) 'tiling[x,y]', real(ifactorx), real(ifactory)
        do j=1,n_slices
          counter = 0
          do atom=1,nat_slice_unitcell(target_atom,j)     !loop over ionization species in unit cell
            if(pw_illum) then                           !if plane wave we only need to shift over the unit cell
              x_coord = tau_slice_unitcell(1,target_atom,atom,j)*ifx
              y_coord = tau_slice_unitcell(2,target_atom,atom,j)*ify
              counter = counter+1
              ion_tau(:,counter) = [x_coord,y_coord]
              call make_shift_oned(Hn0_shiftx_coord(:,counter,j),nopix,x_coord) !this is shifting fractionally on the supercell
              call make_shift_oned(Hn0_shifty_coord(:,counter,j),nopiy,y_coord) !this is shifting fractionally on the supercell
            else
              do ss_atomx=1,ifactorx                !tile the unit cell atom over the supercell
                x_cell = real(ss_atomx-1,fp_kind)
                do ss_atomy=1,ifactory              !tile the unit cell atom over the supercell
                  y_cell = real(ss_atomy-1,fp_kind)
                  x_coord = (tau_slice_unitcell(1,target_atom,atom,j)+x_cell)*ifx
                  y_coord = (tau_slice_unitcell(2,target_atom,atom,j)+y_cell)*ify
                  counter = counter+1
                  ion_tau(:,counter) = [x_coord,y_coord]
                  call make_shift_oned(Hn0_shiftx_coord(:,counter,j),nopix,x_coord) !this is shifting fractionally on the supercell
                  call make_shift_oned(Hn0_shifty_coord(:,counter,j),nopiy,y_coord) !this is shifting fractionally on the supercell
                enddo
              enddo
            endif
            ! below is code depricated with the implementation of the intensity-based inelastic wave function filter
!            elseif(all_atoms) then                      !if focussed probe we need to shift ionization potentials over supercell
!              do ss_atomx=1,ifactorx                  !tile the unit cell atom over the supercell
!                x_cell = real(ss_atomx-1,fp_kind)
!                do ss_atomy=1,ifactory              !tile the unit cell atom over the supercell
!                  y_cell = real(ss_atomy-1,fp_kind)
!!                  write(*,102) "[xc,yc]", x_cell, y_cell
!                  x_coord = (tau_slice_unitcell(1,target_atom,atom,j)+x_cell)*ifx
!                  y_coord = (tau_slice_unitcell(2,target_atom,atom,j)+y_cell)*ify
!                  counter = counter+1
!                  ion_tau(:,counter) = [x_coord,y_coord]
!                  call make_shift_oned(Hn0_shiftx_coord(:,counter,j),nopix,x_coord) !this is shifting fractionally on the supercell
!                  call make_shift_oned(Hn0_shifty_coord(:,counter,j),nopiy,y_coord) !this is shifting fractionally on the supercell
!!                  write(*,102) " UC-pos", tau_slice_unitcell(1,target_atom,atom,j), &
!!                        & tau_slice_unitcell(2,target_atom,atom,j)
!!                  write(*,102) " SC-pos", x_coord, y_coord
!
!                enddo
!              enddo
!            else      !user decided to only shift ionization potentials over a fraction small than the supercell
!              do ss_atomx=-xcells/2,xcells/2         !tile the unit cell atom over the supercell
!                if (ss_atomx.ge.0) then
!                  x_cell = real(ifactorx+ss_atomx,fp_kind)
!                else
!                  x_cell = real(ss_atomx,fp_kind)
!                endif
!                do ss_atomy=-ycells/2,ycells/2       !tile the unit cell atom over the supercell
!                  if (ss_atomy.ge.0) then
!                    y_cell = real(ifactory+ss_atomy,fp_kind)
!                  else
!                    y_cell = real(ss_atomy,fp_kind)
!                  endif
!                  y_cell = real(ifactorx+ss_atomy,fp_kind)
!                  x_coord = (tau_slice_unitcell(1,target_atom,atom,j)+x_cell)*ifx
!                  y_coord = (tau_slice_unitcell(2,target_atom,atom,j)+y_cell)*ify
!                  counter = counter+1
!                  ion_tau(:,counter) = [x_coord,y_coord]
!                  call make_shift_oned(Hn0_shiftx_coord(:,counter,j),nopix,x_coord) !this is shifting fractionally on the supercell
!                  call make_shift_oned(Hn0_shifty_coord(:,counter,j),nopiy,y_coord) !this is shifting fractionally on the supercell
!                enddo
!              enddo
!            endif
          enddo ! atom in slice
        enddo ! slices
        write(*,*)
		write(*,*) 'Finished with preparing transition potentials.'
        write(*,*)
    end function
    
    
    function fourier_shift_array(arrayin,yshift,xshift) result(arrayout)
    
        use CUFFT_wrapper
    
        implicit none
    
        real(fp_kind),intent(in)::arrayin(:,:),yshift,xshift
    
        integer*4::y,x
        complex(fp_kind),allocatable::array(:,:),shift_arrayy(:)
        complex(fp_kind),allocatable::shift_arrayx(:),shift_array(:,:)
        real(fp_kind),allocatable::arrayout(:,:)
    
        !Get size of input array
        y = size(arrayin,dim = 1)
        x = size(arrayin,dim = 2)
    
        !allocate working arrays
        allocate(array(y,x),shift_arrayy(y),shift_arrayx(x))
        allocate(shift_array(y,x),arrayout(y,x))
    
        !Get the FFT of the input array
        array = cmplx(arrayin,kind = fp_kind)
        call fft2(y,x,array,array)
    
        !Make the 1d shift arrays
        call make_shift_oned(shift_arrayy,y,yshift/y)
        call make_shift_oned(shift_arrayx,x,xshift/x)
    
        !Spread 1d shift arrays to make 2d shift array
        shift_array = spread(shift_arrayx,dim =1, ncopies = x)
        shift_array = shift_array*spread(shift_arrayy,dim = 2,ncopies = y)
    
        !Apply shift array in Fourier space
        array=array*shift_array
    
        !Inverse fourier transform
        call ifft2(y,x,array,array)
        arrayout = real(array,kind = fp_kind) 
    
    end function
	   
	subroutine choose_individual_targets(cells_out,targets,n_slices,n_cells,natoms_slice_total,ifactorx,ifactory,ion_tau,sumnatoms_slice_total)
		use output, only:nlines, array_from_txt_file
		implicit none
        
        integer*4::sumnatoms_slice_total
		logical::targets(sumnatoms_slice_total)
		integer*4::cells_out(n_cells),n_slices,n_cells,natoms_slice_total(n_slices,ifactorx,ifactory),ifactorx,ifactory
		real(fp_kind),allocatable::ion_tau(:,:)

		integer*4 :: k,numlines
		logical :: ifile_exists

		allocate(ion_tau(3,sum(natoms_slice_total)))

		write(*,*)
		write(*,*) "The program will search for a file named 'ionization_cells.txt'."
		write(*,*) "This will tell the program which cells you want to output, if this file"
		write(*,*) "does not exist then all ionizations will be outputted"
	
		inquire(file='ionization_cells.txt',exist = ifile_exists)
	
		if(ifile_exists) then
			numlines = min(nlines('ionization_cells.txt'),n_slices*n_cells)
			cells_out = 0
			cells_out(:numlines) = nint(array_from_txt_file('ionization_cells.txt',numlines))
		else
			cells_out = 1
		end if
		write(*,*)
		write(*,*) "The program will search for a file named 'targets.txt'."
		write(*,*) "This will tell the program which target atoms you want to output, if this file"
		write(*,*) "does not exist then all ionizations will be outputted"
		write(*,*)
		inquire(file='targets.txt',exist = ifile_exists)
		
		if(ifile_exists) then
			numlines = min(nlines('targets.txt'),sum(natoms_slice_total))
			targets = 0
			targets(:numlines) = nint(array_from_txt_file('targets.txt',numlines))==1
		else
			targets = .true.

		end if

		write(*,*)
		write(*,*) 'Diffraction patterns will be outputted for the target atoms'
		write(*,*) 'located at the following fractional coordinates).'
		do k=1,size(targets)
			if(targets(k)) write(*,800) k,char(10), mod(ion_tau(1,k)*ifactorx,1.0_fp_kind), mod(ion_tau(2,k)*ifactory,1.0_fp_kind), floor(ion_tau(1,k)*ifactorx,fp_kind)+1,ifactorx, floor(ion_tau(2,k)*ifactory,fp_kind)+1,ifactory
		enddo
		write(*,*)
800     format('Target no. = ',I3,':',a1,'  x = ',f5.2,'  y = ',f5.2,'  unit cell x: ', i2,'/',i2,'  unit cell y: ', i2,'/',i2)
    end subroutine

	 
    !Performs a circular shift average of the istem image     
    function tile_out_image(arrayin,tiley,tilex) result(arrayout)
    
        implicit none
        
        real(fp_kind),intent(in):: arrayin(:,:)
        integer*4,intent(in)::tilex,tiley

        integer*4:: x,y,i,j,deltax,deltay
        logical::even
        real(fp_kind)::dx,dy,shifty,shiftx
        real(fp_kind):: arrayout(size(arrayin,dim = 1),size(arrayin,dim = 2))
    
        y = size(arrayin,dim = 1)
        x = size(arrayin,dim = 2)
    
        !Check to see if pixel dimension is divisible by tiling
        even = (mod(x,tilex) == 0) .and. (mod(y,tiley) == 0)
    
        !Allocate and initialise the arrayout
        arrayout = 0.0_fp_kind
    
        if(even) then
            !Calculate integers that the array will be shifted by
            deltay = y/tiley
            deltax = x/tilex
                
            !Sum over all possible shifts
            do i=1,tiley    
            do j=1,tilex
                arrayout = arrayout + cshift(cshift(arrayin,deltay*i,1)&
												    &      ,deltax*j,2)
            enddo
            enddo
        
        else
            !Calculate non-integer shift values
            dy = real(y,fp_kind)/real(tiley,fp_kind)
            dx = real(x,fp_kind)/real(tilex,fp_kind)
        
            !Sum over all possible circular shifts
            do i=1,tiley    
            do j=1,tilex
                shifty = (i-1)*dy
                shiftx = (j-1)*dx
                arrayout = arrayout + fourier_shift_array(arrayin,shifty,shiftx)
            enddo
            enddo     
        endif
		arrayout= arrayout/ifactory/ifactorx
    end function
    
    
   
    end module
