!--------------------------------------------------------------------------------
!
!  Copyright (C) 2017-2025  L. J. Allen, H. G. Brown, A. J. D’Alfonso, &
!                           S.D. Findlay, B. D. Forbes, J. Barthel
!
!  Additional modifications:
!   2025-05-23 JB - custom ionization form-factor support
!   2025-06-07 JB - SE-imaging option
!   2025-06-12 JB - table of atomic and ionic radii added
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
 
module m_potential

    use m_precision, only: fp_kind
    implicit none
    
    integer(4) :: num_ionizations
    integer(4),allocatable::  atm_indices(:)
    character(20),allocatable::ion_description(:) ! increased to 20 chars for custom ionization lines (2025-05-21 JB)
      
    complex(fp_kind), allocatable :: ionization_mu(:,:,:)        !the ionization scattering factor array, calculated on the grid (supercell)
    complex(fp_kind), allocatable :: fz_adf(:,:,:,:)        !the adf scattering factor array, calculated on the grid (supercell)
    real(fp_kind),    allocatable :: adf_potential(:,:,:,:)
    real(fp_kind),    allocatable :: ionization_potential(:,:,:,:)
    real(fp_kind),    allocatable :: eels_correction_detector(:,:)
    real(fp_kind),    allocatable :: se_transf_aty_radius(:) ! list of effective radii for each atom type in the structure for SE transmission functions
    real(fp_kind),    allocatable :: se_transmission(:,:,:) ! SE transmission coefficients for each slice of the supercell
	!complex(fp_kind), allocatable :: inverse_sinc_new(:,:)
    
    integer(4) :: n_qep_grates,n_qep_passes,nran ! Start of random number sequence
    
    logical :: phase_ramp_shift
    logical(4) :: quick_shift
    
    real(fp_kind),parameter :: thr_cus_ekv = 1.0_fp_kind ! threshold for matching ekv (in keV) in custom ionization (2025-05-22 JB)
    real(fp_kind),parameter :: thr_se_transm = 0.01_fp_kind ! transmission threshold for stopping the SE calculation (2025-06-26 JB)
    real(fp_kind), parameter :: atomicRadius(99) = (/ &
        &  25.0, 120.0, 145.0, 105.0,  85.0,  70.0,  65.0,  60.0,  50.0, 160.0, & ! 1–10
        & 180.0, 150.0, 125.0, 110.0, 100.0, 100.0, 100.0,  71.0, 220.0, 180.0, & ! 11–20
        & 160.0, 140.0, 135.0, 140.0, 140.0, 140.0, 135.0, 135.0, 135.0, 135.0, & ! 21–30
        & 130.0, 125.0, 115.0, 115.0, 115.0,  88.0, 235.0, 200.0, 180.0, 155.0, & ! 31–40
        & 145.0, 145.0, 135.0, 130.0, 135.0, 140.0, 160.0, 155.0, 155.0, 145.0, & ! 41–50
        & 145.0, 140.0, 140.0, 108.0, 260.0, 215.0, 195.0, 185.0, 185.0, 185.0, & ! 51–60
        & 185.0, 185.0, 185.0, 180.0, 175.0, 175.0, 175.0, 175.0, 175.0, 175.0, & ! 61–70
        & 175.0, 155.0, 145.0, 135.0, 135.0, 130.0, 135.0, 135.0, 135.0, 150.0, & ! 71–80
        & 190.0, 180.0, 160.0, 190.0, 127.0, 120.0, 348.0, 215.0, 195.0, 180.0, & ! 81–90
        & 180.0, 175.0, 175.0, 175.0, 175.0, 176.0, 170.0, 186.0, 186.0         & ! 91–99
        /) ! Atomic radii in pm, 1-99 from https://en.wikipedia.org/wiki/Atomic_radii_of_the_elements_(data_page)#Atomic_radius
    real(fp_kind), parameter :: bondRadius(99) = (/ &
        &  32.0,  64.0, 133.0, 102.0,  85.0,  75.0,  71.0,  63.0,  64.0,  67.0, & ! 1–10
        & 155.0, 139.0, 126.0, 116.0, 111.0, 103.0,  99.0,  96.0, 196.0, 171.0, & ! 11–20
        & 148.0, 136.0, 134.0, 122.0, 119.0, 116.0, 111.0, 110.0, 112.0, 118.0, & ! 21–30
        & 124.0, 121.0, 121.0, 116.0, 114.0, 117.0, 210.0, 185.0, 163.0, 154.0, & ! 31–40
        & 147.0, 138.0, 128.0, 125.0, 125.0, 120.0, 128.0, 136.0, 142.0, 140.0, & ! 41–50
        & 140.0, 136.0, 133.0, 131.0, 232.0, 196.0, 180.0, 163.0, 176.0, 174.0, & ! 51–60
        & 173.0, 172.0, 168.0, 169.0, 168.0, 167.0, 166.0, 165.0, 164.0, 170.0, & ! 61–70
        & 162.0, 152.0, 146.0, 137.0, 131.0, 129.0, 122.0, 123.0, 124.0, 133.0, & ! 71–80
        & 144.0, 144.0, 151.0, 145.0, 147.0, 142.0, 223.0, 201.0, 186.0, 175.0, & ! 81–90
        & 169.0, 170.0, 171.0, 172.0, 166.0, 166.0, 168.0, 168.0, 165.0         & ! 91–99
        /) ! Atomic radii in pm, 1-99 from "Molecular Single-Bond Covalent Radii for Elements 1–118"
        ! Pyykkö and Atsumi, Chem. Eur. J. 2009, 15, 186-192, https://doi.org/10.1002/chem.200800987
    
    interface
        subroutine make_site_factor_generic(site_factor, tau)
            use m_precision, only: fp_kind
            complex(fp_kind),intent(out) :: site_factor(:, :)
            real(fp_kind),intent(in) :: tau(:,:)
        end subroutine make_site_factor_generic
    end interface
    
    contains
    
    subroutine prompt_high_accuracy
        
        use m_user_input, only: get_input
        use global_variables, only: high_accuracy
        use m_string
        
        implicit none
        
        integer :: i
        
        call command_line_title_box('Potential calculation method')
        write(*,*) 'Two choices are available for the calculation of potentials.'
        write(*,*) 'The reciprocal space method is accurate but may be slower.'
        write(*,*) 'The hybrid method due to Van den Broek et al. is faster but'
        write(*,*) 'is an approximate approach. Note that if "on-the-fly"'
        write(*,*) 'scattering potentials are used the calculation defaults to the'
        write(*,*) 'hybrid approach. '
        write(*,*) '(Van den Broek et al., Ultramicroscopy 158 (2015) pp. 89-97)'
        write(*,*)
        write(*,*) 'Note: if there is insufficient GPU memory, this choice will be overridden.'
        write(*,*)
    10  write(*,*) 'Please choose a method:'
        write(*,*) '<1> Reciprocal space (accuracy)'
        write(*,*) '<2> Hybrid (speed)'
        call get_input('Scattering factor accuracy', i)   
        write(*,*) 
    
        if (i.eq.1) then
            high_accuracy = .true.
            
        elseif (i.eq.2) then
            high_accuracy = .false.
            
        else
            goto 10
            
        endif
        
    end subroutine prompt_high_accuracy
    

    
    subroutine precalculate_scattering_factors
        
        use m_crystallography
	    use m_precision, only: fp_kind
        use global_variables
        use m_electron, only: elsa_ext,peng_ionic_ff,element
        use m_absorption!, only: complex_absorption, setup_absorptive_array, max_int, delta_kstep, tdsbr, fz_abs,calculate_absorption_mu,include_absorption
		use m_numerical_tools, only: cubspl,ppvalu
        use output
        use m_string
  
	    implicit none
    
        integer(4) :: i, j, k,Z
        real(fp_kind) :: xkstep, temp,el_scat,ax,ay,g2,s2,sky,skx
        real(fp_kind),allocatable :: xdata(:),tdsbrc(:,:,:) 
        real(fp_kind) :: factor, eps, g_vec_array(3,nopiy,nopix)
        
        write(*,*) 'Precalculating scattering factors ...'
  
        if(allocated(inverse_sinc)) deallocate(inverse_sinc)
        if(allocated(fz)) deallocate(fz)
        if(allocated(fz_DWF)) deallocate(fz_DWF)
    
        allocate(fz(nopiy,nopix,nt))
        allocate(inverse_sinc(nopiy,nopix))
        allocate(fz_DWF(nopiy,nopix,nt))
    
        ax = (a0(1)*float(ifactorx))/(float(nopix)*2.0_fp_kind)
        ay = (a0(2)*float(ifactory))/(float(nopiy)*2.0_fp_kind)
        factor = 1.0_fp_kind
        eps = tiny(0.0_fp_kind)
        call make_g_vec_array(g_vec_array,ifactory,ifactorx)
        
	    do j = 1, nopix;do i = 1, nopiy;
            skx = trimr([g_vec_array(1,i,j),0.0_fp_kind,0.0_fp_kind],ss)
            sky = trimr([0.0_fp_kind,g_vec_array(2,i,j),0.0_fp_kind],ss)
            g2 =  trimr(g_vec_array(:,i,j),ss)**2
            s2 = g2 / 4.0_fp_kind
                
            do k = 1, nt
                ! Multiply by fractional occupancy
                if (.not. ionic) el_scat = elsa_ext(nt,k,atomf,s2) * atf(2,k)   
			    if(ionic) el_scat = Peng_ionic_FF(s2,nint(atf(1,k)),dZ(k)) * atf(2,k)        
                    
                
	            ! Fill the potential matrix. Note: these are U(g)/2K
                fz(i,j,k) = cmplx( el_scat, 0.0_fp_kind ,fp_kind)
                fz_DWF(i,j,k) = cmplx( exp( -tp**2.0_fp_kind*g2*atf(3,k) / 2.0_fp_kind ), 0.0_fp_kind,fp_kind ) 
            enddo
            
            !Sinc 
            inverse_sinc(i,j) = cmplx((tp*skx*ax+eps)/(sin(tp*skx*ax)+eps)*(tp*sky*ay+eps)/(sin(tp*sky*ay)+eps),0.0_fp_kind,fp_kind)
        enddo; enddo
        ! Currently have U(g)/2K, so multiply by 2K
        fz = 2*ak*fz
        
        ! Normalise the sinc function
        inverse_sinc = inverse_sinc*float(nopiy)*float(nopix)
        
        write(*,*)
    
    end subroutine precalculate_scattering_factors
    

    
    subroutine make_site_factor_matmul(site_factor, tau)

        use m_precision, only: fp_kind
        use global_variables, only: nopiy, nopix, tp
        use m_crystallography
    
	    implicit none
    
        !output
        complex(fp_kind),intent(out) :: site_factor(:, :)
        
        !input
        real(fp_kind),intent(in) :: tau(:,:)
    
        integer :: i, j
    
        integer :: g_vec_array(3,nopiy,nopix)
        
        call make_g_vec_array(g_vec_array)
        
        !$OMP PARALLEL PRIVATE(i, j)
        !$OMP DO
	    do i = 1, nopiy
            do j = 1, nopix
                site_factor(i,j) = sum(exp(cmplx(0.0_fp_kind, -tp*matmul(g_vec_array(:,i,j), tau), fp_kind)))
            enddo
        enddo
	    !$OMP END DO
        !$OMP END PARALLEL
        
    end subroutine make_site_factor_matmul
    
    
#ifdef GPU    
    subroutine make_site_factor_cuda(site_factor, tau)

        use m_precision, only: fp_kind
        use global_variables, only: nopiy, nopix, tp
        use cuda_potential, only: cuda_site_factor
        use cuda_array_library, only: blocks, threads
		use m_crystallography,only:make_g_vec_array
        
	    implicit none
    
        !output
        complex(fp_kind),intent(out) :: site_factor(:, :)
        
        !input
        real(fp_kind),intent(in) :: tau(:,:)
    
        integer :: g_vec_array(3,nopiy,nopix)
    
        integer,device :: g_vec_array_d(3,nopiy,nopix)
        real(fp_kind),device :: tau_d(3, size(tau,2))
        complex(fp_kind),device :: site_factor_d(nopiy,nopix)
                          
        call make_g_vec_array(g_vec_array)
		g_vec_array_d = g_vec_array
        tau_d = tau
        
        call cuda_site_factor<<<blocks,threads>>>(site_factor_d, tau_d, g_vec_array_d, nopiy, nopix)        

        site_factor = site_factor_d
    end subroutine make_site_factor_cuda
#endif    
    
    subroutine make_site_factor_hybrid(site_factor, tau)

        use m_precision, only: fp_kind
        use global_variables, only: nopiy, nopix,inverse_sinc
        use CUFFT_wrapper, only: fft2
        
        implicit none
    
        complex(fp_kind),intent(out) :: site_factor(:,:)
        
        real(fp_kind),intent(in) :: tau(:,:)
        
        integer :: j
        integer :: xpixel, ypixel
        real(fp_kind) :: xpos, ypos, fracx, fracy
    
        site_factor = 0.0_fp_kind
        
        !$OMP PARALLEL PRIVATE(xpos, ypos, j, xpixel,ypixel,fracx,fracy)
        !$OMP DO
        do j = 1, size(tau, 2)
            xpos = tau(1,j)*nopix
            ypos = tau(2,j)*nopiy
            
            ! Ensure that the pixel positions are in range
            
            if (ceiling(xpos).gt.nopix) then
                xpos = xpos - float(nopix)
            elseif (floor(xpos).lt.1) then
                xpos = xpos + float(nopix)
            endif
            
            if (ceiling(ypos).gt.nopiy) then
                ypos = ypos - float(nopiy)
            elseif (floor(ypos).lt.1) then
                ypos = ypos + float(nopiy)
            endif
            
            !fraction of the pixel top right
            xpixel = ceiling(xpos)
            ypixel = ceiling(ypos)
            fracx = mod(xpos, 1.0_fp_kind)
            fracy = mod(ypos, 1.0_fp_kind)
            
            call pixel_check(xpixel, ypixel)
            site_factor(ypixel,xpixel) = site_factor(ypixel,xpixel) + fracx*fracy
            
            !fraction of the pixel top left
            xpixel = floor(xpos)
            ypixel = ceiling(ypos) 
            fracx = 1.0_fp_kind - mod(xpos, 1.0_fp_kind)
            fracy = mod(ypos, 1.0_fp_kind)
            
            call pixel_check(xpixel, ypixel)
            site_factor(ypixel,xpixel) = site_factor(ypixel,xpixel) + fracx*fracy
            
            !fraction of the pixel bottom right
            xpixel = ceiling(xpos)
            ypixel = floor(ypos)
            fracx = mod(xpos, 1.0_fp_kind)
            fracy = 1.0_fp_kind - mod(ypos,1.0_fp_kind)
            
            call pixel_check(xpixel, ypixel)
            site_factor(ypixel,xpixel) = site_factor(ypixel,xpixel) + fracx*fracy
            
            !fraction of the pixel bottom left
            xpixel = floor(xpos)
            ypixel = floor(ypos)
            fracx = 1.0_fp_kind - mod(xpos, 1.0_fp_kind)
            fracy = 1.0_fp_kind - mod(ypos, 1.0_fp_kind)
            
            call pixel_check(xpixel, ypixel)
            site_factor(ypixel,xpixel) = site_factor(ypixel,xpixel) + fracx*fracy
        enddo
        !$OMP end do
        !$OMP end parallel
        !fix pixel offset
        site_factor = cshift(site_factor,SHIFT = -1,DIM=1)
        site_factor = cshift(site_factor,SHIFT = -1,DIM=2)

        call fft2(nopiy, nopix, site_factor, site_factor)
        site_factor = site_factor * inverse_sinc/sqrt(float(nopiy)*float(nopix))!_new * sqrt(float(nopiy)*float(nopix))
        
        contains
        
        subroutine pixel_check(x, y)
            ! Wrap pixel coordinates around so that they remain in range.
    
            implicit none
    
            integer(4) :: x,y
    
            if(x.eq.0) x = nopix
            if(x.eq.nopix+1) x = 1
            if(y.eq.0) y = nopiy
            if(y.eq.nopiy+1) y = 1
    
        end subroutine pixel_check
    
    end subroutine make_site_factor_hybrid

    subroutine setup_inelastic_ionization_types()
        use global_variables, only: EELS, EDX, SEI
        use m_user_input
        use m_string
        implicit none
        integer*4::i_eels
      
        call command_line_title_box('Ionization')
		i_eels = 0
        do while(i_eels<1.or.i_eels>3) ! 2025-06-07 JB: added SE-Imaging option
		    write(*,*) char(10),' <1> EELS',char(10),' <2> EDX',char(10),' <3> SE-imaging',char(10)
            call get_input('Ionization choice', i_eels)
        enddo
        EELS = i_eels.eq.1
        EDX = i_eels.eq.2
        SEI = i_eels.eq.3
		call local_potential()

    end subroutine
    
    function get_custom_ionization_total() result(num_custom)
    !Returns the total number of custom ionization data sets
    ! in the custom ionization file
    ! 2025-05-21 JB
        implicit none
        integer*4 :: num_custom, reason , lcnt
        character*120 :: filename,line
        num_custom=0
        filename = 'custom_ionization.dat'
		open(unit=35,file=filename,status='old',err=970)
        lcnt = 0
        do
            read(35, '(a)', iostat=reason) line
            if (0==reason) then ! successful read
                lcnt = lcnt + 1
            else ! stop reading
                exit
            end if
        end do
        close(35)
        num_custom = RSHIFT(lcnt,1) ! divide by to (assuming two lines per data set)
970     return ! could not open file, thus no custom ionization data
    end function
    
    
    function get_custom_ionization_num() result(n)
    !Returns the number of matching custom ionization data sets and fills information
    ! with the data from the custom ionization file (custom_ionization.dat)
    ! This is not loading the data, only the number of matching data sets is determined
    ! Use load_custom_ionization_data() to load the data
    !
    ! Warning: The loading in this function assumes that the compiler treats tabs as spaces
    !           (which is not the case for gfortran, but is for ifort, ifx, pgifortran).
    !           The solution is to replace the tabs with spaces in the file or in ldata.
    !           Since the same problem exists with the default data files, we do not handle
    !           this here. Note, char(9) is a tab character. ;)
    !
    ! 2025-05-23 JB
        use global_variables, only: ekv, nt, atf
        implicit none
        integer*4, parameter :: funit = 35 ! unit number for file access
        integer*4 :: n ! returned number of available (matching) custom ionization data sets
        ! temps and parameters from file
        character(10)::lsym,lshell
        integer*4::reason,i,z,lz
        character*120::filename,lhead
        real(fp_kind)::lekv
        ! init
        n = 0
        filename = 'custom_ionization.dat'
        ! open the file
        open(unit=funit,file=trim(filename),status='old',err=970)
        ! start parsing the file
        do ! reading lines loop
            read(funit, '(a)', iostat=reason) lhead ! read next line (limited to 120 chars)
            if (0==reason) then ! successful read
                ! find matching header data (this must be equal to whats in load_custom_ionization(...))
                if (INDEX(lhead, "EEL") > 0) then ! EELS header line
                    ! check matches to atom types and ekv, using first four columns, which should always be there
                    read(lhead, *) lsym, lz, lshell, lekv
                    if (ABS(lekv-ekv)<thr_cus_ekv) then ! beam energy matches
                        do i=1, nt ! loop over structure atom types
                            if (lz==NINT(atf(1,i))) then ! atom type matches
                                n = n + 1 ! we have to make a choice available for each atom type, so accumulate each match
                            end if
                        enddo 
                    end if
                else
                    cycle ! next line
                end if
            else ! stop reading (probably end of file)
                exit
            endif
        end do ! reading lines loop
        ! close the file
        close(unit=funit)
970     return ! could not open file, thus no custom ionization data
    end function
    
    
    subroutine load_custom_ionization(num,aty,orb,de,linenum)
    !Loads num data sets of matching custom ionization parameters into the interfaced arrays.
    ! First, call num = get_custom_ionization_num() to get the number n of matching data sets.
    !
    ! This routine loads header data from file 'custom_ionization.dat'.
    ! To load the actual ionization form factors, call the subroutine
    ! load_custom_ionization_data(...) and provide the line number in the file,
    ! which is returned by this routine in array linenum.
    !
    ! 2025-05-23 JB
        use global_variables, only: ekv, nt, atf
        implicit none
        integer*4, parameter :: funit = 36 ! unit number for file access
        integer*4, intent(in) :: num ! number of matching custom ionization data sets
        integer*4, intent(inout) :: aty(num), linenum(num) ! atom type indices for each data set
        character(2), intent(inout) :: orb(num) ! ionization shell symbol for each data set
        real*4, intent(inout) :: de(2,num) ! energy-loss window and ff-parameters
        ! temps and parameters from file
        character(10)::lsym,lshell
        integer*4::reason,i,ii,z,lz,iline
        character*120::filename,lhead,ldum
        character*1014::ldum2
        real(fp_kind)::lekv,lde1,lde2
        ! init
        ii = 0
        iline = 0
        filename = 'custom_ionization.dat'
        ! open the file
        open(unit=funit,file=trim(filename),status='old',err=970)
        ! start parsing the file
        do ! reading lines loop
            read(funit, '(a)', iostat=reason) lhead ! read next line assuming header
            iline = iline + 1 ! increment line number
            if (reason/=0) then ! eof or error, stop reading
                close(unit=funit)
                return
            end if
            read(funit, '(a)', iostat=reason) ldum ! read next line assuming data
            iline = iline + 1 ! increment line number
            if (0==reason) then ! successful data read
                !
                ! find matching header data (this must be equal to what's in get_custom_ionization_num())
                if (INDEX(lhead, "EEL") > 0) then ! EELS header line
                    ! check matches to atom types and ekv, using first four columns, which should always be there
                    read(unit=lhead, fmt=*) lsym, lz, lshell, lekv, ldum, lde1, lde2
                    if (ABS(lekv-ekv)<thr_cus_ekv) then ! beam energy matches
                        do i=1, nt ! loop over structure atom types
                            if (lz==NINT(atf(1,i))) then ! atom type also matches
                                ii = ii + 1 ! increment store index
                                ! store the header data
                                aty(ii) = i ! atom type index
                                orb(ii) = lshell ! ionization shell symbol
                                de(1,ii) = lde1 ! energy-loss window start
                                de(2,ii) = lde2 ! energy-loss window end
                                linenum(ii) = iline ! read the parameters
                                ! This way we only load the data once, but make it availabale to each
                                ! matching atom type of the structure. This keeps the custom table
                                ! consistent with how the default table is used.
                            end if
                        enddo 
                    end if
                else
                    cycle ! next line
                end if
            else ! stop reading (probably end of file)
                exit
            endif
        end do ! reading lines loop
        ! close the file
        close(unit=funit)
970     return ! could not open file, thus no custom ionization data
    end subroutine
    
    
    subroutine load_custom_ionization_data(linenum, params)
    !Loads the custom ionization form factor parameters from the file 'custom_ionization.dat'
    ! This routine loads the parameters from the line number linenum in the file.
    !
    ! Call this subroutine after calling load_custom_ionization(...) to get the line number.
    !
    ! Warning: The loading in this function assumes that the compiler treats tabs as spaces
    !           (which is not the case for gfortran, but is for ifort, ifx, pgifortran).
    !           The solution is to replace the tabs with spaces in the file or in ldata.
    !           Since the same problem exists with the default data files, we do not handle
    !           this here. Note, char(9) is a tab character. ;)
    !
    ! 2025-05-23 JB
        implicit none
        integer*4, parameter :: funit = 37 ! unit number for file access
        integer*4, intent(in) :: linenum ! line number in the file to read from
        real(fp_kind), intent(inout) :: params(29) ! ionization form factor parameters taken from the file at linenum
        integer*4::reason,i
        character*120::filename ! file name buffer
        character*1024::ldata ! data line (1024 chars, should be enough for 29 numbers and separators)
        ! init
        filename = 'custom_ionization.dat'
        ! open the file
        open(unit=funit,file=trim(filename),status='old',err=970)
        ! start parsing the file
        do i=1,linenum! reading lines loop
            read(funit, '(a)', iostat=reason) ldata ! read linenum lines assuming data is in line linenum
            if (reason/=0) then ! eof or error, stop reading
                close(unit=funit)
                return
            end if
        end do
        ! close the file
        close(unit=funit)
        read(ldata,*) params ! read the parameters from the line
970     return ! could not open file, thus no custom ionization data
    end subroutine
    
	
    function get_ionization_shell_line(shell,atno) result(lineno)
        !Get the line for the given ionization shell and atom number,returns 
        !-1 if that shell is not contained in the parameterization
        use m_string
        character(2),intent(in)::shell
		integer*4,intent(in)::atno
        
        character*31::filename,line
        character(len=:),allocatable ::string 
        
        integer*4::reason,lineno,l
    
        !open the pertinent data files
        ! 2025-06-12 JB: modified path seprator to '/', should work with Windows and Linux
        filename = 'ionization_data/EELS_EDX_'//shell//'.dat'
		open(unit=35,file=filename,status='old',err=970)
        
        l = len('Z = '//to_string(int(atno))//' ')
        allocate(character(l)::string)
        string = 'Z = '//to_string(int(atno))//' '
        lineno = 1
        DO
           READ(35,960,IOSTAT=Reason) line
960        format(a31)
           IF (Reason > 0)  THEN
                write(*,*) 'Problem reading ',filename           
           !If end of file is reach return -1 (shell is not in parametrization)
           ELSE IF (Reason < 0) THEN
              lineno = -1
              close(35)
              return
           ELSE
              !If the substring Z = atno, then return this line 
              if(index(line,string)>0) then
                  close(35)
                  return
              endif
              lineno = lineno+1
           END IF
        END DO
        close(35)
        return
970 write(*,*) 'Problem reading ',filename           
  end function
        
  function get_ionization_parameters(shell,atno,DE,EDX) result (EELS_EDX_params)
  !Read ionisation form factor parameters from ionization_data files
  !shell should be a string describing the orbital ie, '1s', '2s', '2p' etc
  !atno is the atomic number
  !DE is the energy window for EELS and is ignored if the EDX parameterization is requested
  !EDX is a boolean variable, pass .TRUE. for EDX parameterization and .FALSE. for EELS
  !
  ! 2025-May-14, modified to run linear eels range interpolation on command line option linpoleels
    use m_numerical_tools
    use m_string
    use global_variables,only: ekv,ak1,nt,atf,substance_atom_types,linpoleels
        
    character(2),intent(in)::shell
    integer*4,intent(in)::atno
	real(fp_kind),intent(in)::DE
	logical,intent(in):: EDX
		
	real(fp_kind)::params(29,8,5),EELS_PARAM_SET2(5,29),xdata(8),bscoef_(4,8)
    real(fp_kind) :: dedata(5),bscoef2_(4,5),p(29),EELS_EDX_params(29)
        
    integer*4::iatom,atno_check,i,ii,m,ishell,iz,lineno,n,mm
	character(10) junk
	character(1) cjunk1,cjunk2,cjunk3

	m = 5
	if(EDX) m=1
		
	!write(*,*) 'reading inelastic scattering parameterization:'
	!write(*,*) '- file: '//'ionization_data/EELS_EDX_'//shell//'.dat'
        
    n = str2int(shell(1:1))
    lineno = get_ionization_shell_line(shell,atno)
	!write(*,*) '- shell: '//shell//' - id = '//to_string(n)
	!write(*,*) '- from line: '//to_string(lineno)
		
    !open the pertinent data files and read to relevant line 
    ! 2025-06-12 JB: modified path seprator to '/', should work with Windows and Linux
	open(unit=16,file='ionization_data/EELS_EDX_'//shell//'.dat',status='old',err=970)
	do iz = 1,lineno
	  read(16,*) junk
    enddo
    p = 0
    !Later parametrizations only contain 28 datapoints
    mm=29;if (n>2) mm=28
    !write(*,*) '- number of items per line: '//to_string(mm)
    !Read parameters
    do i=1,8 !Loop over accelerating voltages
	  read(16,*) junk ! E=xx kev header
	  do ii=1,6 !Loop over energy loss above threshhold (EELS) and EDX 
	    read(16,*) p(1:mm)
		!write(*,*) '--- #'//to_string(ii)//': '//to_string(p(1))
		if((.not.EDX).and.(ii<6)) params(:,i,ii) = p(:) !EELS is the first 5 lines
		if(EDX.and.(ii==6)) params(:,i,1) = p(:) !EDX is the last line
      enddo
    enddo
    !Can close parameters file now
    close(16)
    
    
    !Interpolate to accelerating voltage used
    !data in files is in steps of 50 keV with 8 points
    !this is stored in xdata
    xdata =(/(i*50, i=1,8,1)/)
    do ii=1,m
      do i=1,mm
		bscoef_(1,:) = params(i,1:8,ii)
		call cubspl(xdata, bscoef_(:,:), 8, 0, 0)
		EELS_param_set2(ii,i) = ppvalu(xdata,bscoef_(:,:),7,4,ekv,0)
      enddo
    enddo
	   
	!If EDX then no energy window interpolation is needed
	if (EDX) then
	  EELS_EDX_params = EELS_param_set2(1,:)
      return
    endif
	  
    ! contained within EELS_param_set2(i,ii) is the 29 data points (first index) interpolated
    ! to the correct incident energy there are 5 rows Interpolate to energy window desired
    dedata = real([1,10,25,50,100],kind=fp_kind)

    !f(s)/DE is mostly flat and interpolates more simply
    EELS_EDX_params=0
    do i=1,mm
      if (linpoleels) then ! linear interpolation mode, 2025-May-14, JB
        bscoef2_(1,1:m) = EELS_param_set2(1:m,i) ! using m instead of fix number 5 (2025-05-21 JB)
        EELS_EDX_params(i) = linpol(DE, dedata, bscoef2_(1,:), m)
      else
		do ii=1,m
		  bscoef2_(1,ii) = EELS_param_set2(ii,i) / dedata(ii) ! setup for cupspl with coeffs/de
		enddo
		call cubspl(dedata, bscoef2_(:,:), m, 0, 0) ! cubic spline coefficients for m points
        EELS_EDX_params(i) = DE*ppvalu(dedata,bscoef2_(:,:),m-1,4,DE,0)
      endif
    enddo
		
	return
970 write(*,*) ' Cannot access data file EELS_EDX_'//shell//'.dat'
	stop
  end function

  
  !********************************************************************************
  !     subroutine EELS_local_potential()
  !     reads in the scattering factors and performs the interpolation
  !     necessary to accommodate arbitrary geometries and energy windows
  !     
  !     2020-11-26 / JB / support for more orbitals added (4f, 5s, 5p)
  !                       now 12 files named 'EELS_EDX_{orb}.dat
  !                       are expected in the ionization_data folder
  !     2025-05-14 / JB / added linear interpolation option for EELS range
  !     2025-05-23 / JB / added custom ionization parameters to EELS
  !                       this requires the file custom_ionization.dat to be
  !                       present in the execution directory
  !     2025-06-07 / JB / added SE imaging, flags EDX and SEI are now in
  !                       global_variables module, and set before this subroutine
  !     
  !********************************************************************************
  subroutine local_potential()

    use m_string
    use m_numerical_tools, only: cubspl,ppvalu
    use m_multislice
    use global_variables
    use m_user_input

    implicit none

    integer(4) i,ii,iii,j,kval,m,nchoices,ZZ,nshells,norbitals,k
    integer(4),allocatable::available_shells(:),available_atoms(:)

    real(fp_kind),allocatable:: DE(:)
    real(fp_kind) eels_inner,eels_outer
    character(2) shell_name_EELS(12),orb
    character(3) shells(12) !(9)
    !character(13):: contributions(4)
    character(32)::tmpi
    logical,allocatable::choices(:)
    logical::k_shell,l_shell,EDXpresent(nt,4)
    
    ! addional variables for custom ionization parameters (2025-05-21 JB)
    integer*4 :: ncustom ! number of available custom parameter sets
    logical, allocatable :: choicus(:) ! menu selection for custom parameters
    integer*4, allocatable :: cus_aty(:) ! structure atom type indices for custom parameters
    character(2), allocatable :: cus_orb(:) ! orbital codes for custom parameters
    real(fp_kind), allocatable :: cus_de(:,:) ! energy windows for custom parameters
    integer*4, allocatable :: cus_line(:) ! custom ionization parameter line number in the parameterization file

    shell_name_EELS = [ '1s', '2s', '2p', '3s', '3p', '3d', '4s', '4p', '4d', '4f', '5s', '5p']
    shells =          [  'K', 'L1','L23', 'M1','M23','M45', 'N1','N23','N45','N67', 'O1','O23']
    !contributions = ['1s','2s and 2p','3s, 3p and 3d','4s, 4p and 4d']
    norbitals = size(shell_name_EELS)
    nshells = size(shells)
    
	!If EELS, setup detector
	if(EELS) then
	    write(*,91)
91      format(1x,'The EELS calculations assume the local approximation, which may be', /, &
			  &1x,'inappropriate when the EELS detector does not have a very large', /, &
              &1x,'acceptance angle. To account for the finite detector size, a correction', /, &
              &1x,'is applied. For more details see Y. Zhu et al. APL 103 (2013) 141908.', /)
		
        eels_inner = 0.0_fp_kind
		  
        write(*,*) 'EELS detector outer angle (mrad):'
        call get_input('Outer EELS angle', eels_outer)
		  
        write(*,*)

        eels_inner = ak1*tan(eels_inner/1000.0_fp_kind)
        eels_outer = ak1*tan(abs(eels_outer)/1000.0_fp_kind)
        if(allocated(eels_correction_detector)) deallocate(eels_correction_detector)
        allocate(eels_correction_detector(nopiy,nopix))
        eels_correction_detector = make_detector(nopiy,nopix,ifactory,ifactorx,ss,eels_inner,eels_outer)
    endif
    
    !If SEI, setup detector, 2025-06-07 JB
	if(SEI) then
	    write(*,92)
92      format(/,1x,'The SE image simulation assumes an isotropic SE emission. A limited', /, &
                &1x,'acceptance angle for the SE detector is taken into account as scaling factor.',/)
		
        eels_inner = 0.0_fp_kind
		  
        write(*,*) 'Enter SE detector acceptance semi angle (mrad):'
        call get_input('SE acceptance angle', eels_outer)
		  
        write(*,*)

        se_det_scale = abs(eels_outer) / 1.256637E+04_fp_kind ! convert to radians and devide by 4pi
        eels_outer = ak1*tan(abs(eels_outer)/1000.0_fp_kind) ! convert to 1/A (only for consistency with EELS)
        
        write(*,93)
93      format(  1x,'The SE detector usually accepts only slow electrons. Their signal', /, &
		        &1x,'quickly attenuates depending on the depth of the ionization event', /, &
		        &1x,'in the material. Effective transmission functions are calculated', /, &
				&1x,'based on atomic radii and the inelastic mean free path for slow', /, &
                &1x,'secondary electrons. Figures 1(a) and (d) in  Seah & Dench,', /, &
                &1x,'Surf. and Interf. Anal. 1 (1979) 2-11 could help to estimate a value.',/)
        
        se_imfp = 0.001_fp_kind
        do while(se_imfp < 0.1_fp_kind) ! 2025-06-07 lower limit 0.1 A
            write(*,*) 'Enter the effective inelastic mean free path (>0.1) in A for SE:'
            call get_input('SE inelastic mean free path', se_imfp)
            se_imfp = ABS(se_imfp) ! ensure positive value
        enddo
        
        write(*,*)
        
        ! setup table of effective atomic radii used for SE transmission
        allocate(se_transf_aty_radius(nt))
        do i = 1, nt
            ZZ=nint(ATF(1,i)); se_transf_aty_radius(i) = bondRadius(ZZ) ! use bond radius as default
        enddo
        ! implement a menu that allows the user to change the effective atomic radius
        write(*,94)
94      format(  1x,'The SE absorption in the structure is modelled by effective radii', /, &
		        &1x,'defining a Gaussian like peak at the equilibrium position for all', /, &
		        &1x,'atoms in the input structure. The radii are preset from a table of', /, &
				&1x,'bond radii, see Pyykkoe and Atsumi, Chem. Eur. J. 2009, 15, 186-192.', /, &
                &1x,'These radii can be modified now:')
95      format('Index  Atom| Z  | Radius (pm)')
96      format('------------------------------')
97      format(1x,'<',i2,'>',2x,a2,2x,'|',1x,i2,1x,'|',1x,f6.1)
        kval = -1
        do while (kval.ne.0)
            write(*,*) char(10),'Radii used for SE transmission functions:',char(10)
            write(*,95)
            write(*,96)
            do i=1, nt
                write(*,97) i,trim(adjustl(substance_atom_types(i))), &
                    & nint(atf(1,i)), se_transf_aty_radius(i)
            enddo
            write(*,120)
            write(*,96)
            call get_input('Select to change radius <0> continue', kval)
            if ((kval>0).AND.(kval<=nt)) then
                ! get new radius for the selected atom type
                write(*,*) 'Enter new effective atomic radius in pm:'
                call get_input('Set atom radius', se_transf_aty_radius(kval))
            endif
        enddo
        write(*,*)
        
    endif
    
    ! Setup for default ionization parameters
    !
    ! Count available orbitals by checking what is available for the given atoms in
    ! the parametrization files
    ii = 0
    do i = 1, nt; ZZ=nint(ATF(1,i)); do j=1,norbitals
        if(get_ionization_shell_line(shell_name_EELS(j),ZZ)>-1) ii=ii+1 ! count default ionization data sets
    enddo; enddo
          
    ! allocations to hold user selections
    nchoices = ii
    allocate(choices(ii),DE(ii))
    DE = 0
    choices = .false.
    num_ionizations = 0
    
    do while (num_ionizations.eq.0) ! loop until user selects ionization shells
	
    kval = -1
    write(*,*) char(10),' muSTEM calculates signals for ionization of electrons to the continuum.'
    write(*,*) 'At this point bound->bound (white line) transitions are not taken'
    write(*,*) 'into account.',char(10)
    write(*,*) 'K, L, M and N shell ionizations are available though users should be aware '
    write(*,*) 'that quantitative agreement between simulation and theory has only been '
    write(*,*) 'demonstrated for K and L shells (see Y. Zhu and C. Dwyer, Microsc. Microanal.'
    write(*,*) '20 (2014), 1070-1077)',char(10)
    
    do while (kval.ne.0) ! Default table EELS selection menu loop
        
100 format(/,' Ionization choices',/,/,'Index  Atom| Z  |',a,'| Included(y/n)')
101 format(  '---------------------------------------------------------------')
    if(EDX) write(*,100) ' Orbital | Shell '
110 format(1x,'<',i2,'>',2x,a2,2x,'|',1x,i2,1x,'|',2x,a2,5x,'|',1x,a3,3x,'|',1x,a1,6x)
    
    if(.not.EDX) write(*,100) ' Orbital | Shell | Window (eV)'
    write(*,101)
111 format(1x,'<',i2,'>',2x,a2,2x,'|',1x,i2,1x,'|',2x,a2,5x,'|',1x,a3,3x,'|',1x,f5.1,6x,'|',1x,a1,6x)
120 format(' < 0> continue')
    
    !Display choices for default ionization tables
    ii=1
    do i = 1, nt;ZZ=nint(ATF(1,i)); do j=1,norbitals
      if(get_ionization_shell_line(shell_name_EELS(j),ZZ)>-1) then
        if(EDX) write(*,110) ii,trim(adjustl(substance_atom_types(i))), &
                 & int(ZZ),shell_name_EELS(j),shells(j),logical_to_yn(choices(ii))
        if(.not.EDX) write(*,111) ii,trim(adjustl(substance_atom_types(i))), &
                 & int(ZZ),shell_name_EELS(j),shells(j),DE(ii),logical_to_yn(choices(ii))
        ii=ii+1
      endif
    enddo;enddo
    
    write(*,120)
    write(*,101)
    call get_input('Shell choice <0> continue', kval)
	!Update choice
	if ((kval.gt.0).and.(kval.le.nchoices)) then
        choices(kval) = .not.choices(kval)
        
        !If EELS or SEI get energy window
        if (.not.EDX) then
          if (EELS) write(*,*) 'Define the electron energy-loss window for the selected shell.'
          if (SEI) write(*,*) 'Define the energy range accepted by the SE detector.'
          DE(kval) =-1
          do while ((DE(kval).lt.1).or.(DE(kval).gt.100)) ! 2025-05-05 lower limit corrected to 1 eV (JB/LJA)
            write(*,*) 'Enter energy range above ionization threshold in eV (between 1 and 100):',char(10)
            call get_input('Energy window', DE(kval))
          enddo 
        end if
    end if
    num_ionizations = count(choices)
    enddo ! eels selection menu
    
    
    ! Custom EELS ionization data initializtation (currently not allowed for EDX)
    ncustom = get_custom_ionization_num() ! get matching number of custom ionization data sets
    if ((.not.EDX) .and. (ncustom>0)) then
        if (ALLOCATED(cus_aty)) deallocate(cus_aty, cus_orb, cus_de, cus_line, choicus)
        allocate(cus_aty(ncustom), cus_orb(ncustom), cus_de(2,ncustom), &
                & cus_line(ncustom), choicus(ncustom))
        choicus = .false.
        ! load custom ionization header data
        call load_custom_ionization(ncustom, cus_aty, cus_orb, cus_de, cus_line)
        write(unit=tmpi,fmt='(i)') ncustom
        write(*,*) char(10),'Custom ionization selection ('//TRIM(adjustl(tmpi))// &
                & ' data sets available)'
        
        
        ! allow the user to select custom ionization data sets
        kval = -1        
200 format('Index  Atom| Z  | Orbital | Shell | Window (eV)     | Incl.(y/n)')
201 format('----------------------------------------------------------------')
211 format(1x,'<',i2,'>',2x,a2,2x,'|',1x,i2,1x,'|',2x,a2,5x,'|',1x,a3,3x,'|', &
                & 1x,f6.1,' - ',f6.1,1x,'|',1x,a1,6x)
212 format(' <99> toggle all')
                
        do while (kval.ne.0) !custom EELS selection menu loop
            write(*,*) char(10),'Ionization choices',char(10)
            write(*,200)
            write(*,201)
            do i=1, ncustom
                do j=1,norbitals
                    if (0<INDEX(cus_orb(i),shell_name_EELS(j))) exit
                end do
                write(*,211) i,trim(adjustl(substance_atom_types(cus_aty(i)))), &
                    & nint(atf(1,cus_aty(i))), cus_orb(i), shells(j), cus_de(1,i), &
                    & cus_de(2,i), logical_to_yn(choicus(i))
            end do
            write(*,212)
            write(*,120)
            write(*,201)
            call get_input('Shell choice <0> continue', kval)
            if((kval.gt.0).and.(kval.le.ncustom)) choicus(kval) = .not.choicus(kval)
            if(kval.eq.99) choicus = .not.choicus ! toggle all choices
        end do
        
        num_ionizations = count(choices) + count(choicus) ! total number of ionizations (default + custom)
        
    else ! either EDX or ncustom==0
        ncustom = 0 ! no custom ionization data for EDX
    end if
    
    if (num_ionizations.eq.0) then
        write(*,*) char(10),'You need to choose at least 1 ionization.',char(10)
    end if
    
    enddo ! loop user EELS/EDX/SE ionization choices
    
    write(*,*) char(10),'Number of ionization scattering factors:',num_ionizations,char(10)
    allocate(ionization_mu(nopiy,nopix,num_ionizations),atm_indices(num_ionizations), &
                & Ion_description(num_ionizations))
    ionization_mu = 0
    ii=1
    iii=1
    !Now read in default EELS or EDX parameters
    do i = 1, nt; ZZ=nint(ATF(1,i))
      do j=1,norbitals
        if(get_ionization_shell_line(shell_name_EELS(j),ZZ)>-1) then; if(choices(ii)) then
          ionization_mu(:,:,iii) = make_fz_EELS_EDX(shell_name_EELS(j),zz,DE(ii),0)* atf(2,i)*fz_DWF(:,:,i)
          atm_indices(iii) = i ! store atom type index for this ionization, to be used for output
          ion_description(iii) = shell_name_EELS(j) ! store shell name for this ionization, to be used for output
          iii=iii+1;endif; ii=ii+1; endif
      enddo
    enddo
    !... and now for the custom ionization data sets
300 format(a2,'_',I0.4,'-',I0.4,'eV')
    if ((.NOT.EDX).AND.(ncustom>0)) then ! only if custom EELS data is available and not EDX
    do i=1, ncustom
      if (choicus(i)) then ! only if the user selected this data set
        zz = nint(atf(1,cus_aty(i))) ! atom type index -> atomic number
        ionization_mu(:,:,iii) = make_fz_EELS_EDX(cus_orb(i),zz,cus_de(2,i),cus_line(i))* &
                & atf(2,cus_aty(i))*fz_DWF(:,:,cus_aty(i))
        atm_indices(iii) = cus_aty(i) ! store atom type index for this ionization, to be used for output
        ! store shell name and energy window for this ionization, to be used for output
        write(unit=ion_description(iii),fmt=300) TRIM(adjustl(cus_orb(i))), nint(cus_de(1,i)), nint(cus_de(2,i))
        iii=iii+1
      end if
    end do
    endif
 
  end subroutine
      
  !Subroutine to make the Fz_mu needs to have prefactors accounted for (volume fo the unit cell etc.)
  !needs to be multiplied by the DWF for the pertinent atom type
  !Modified interface by adding lc to allow for custom ionization data (2025-05-21 JB)
  !Modified interface, EDX flag removed as this is now in global_variables module
  function make_fz_EELS_EDX(orbital,zz,DE,lc) result(fz_mu)
	use m_precision
    use global_variables
	use m_numerical_tools, only: cubspl,ppvalu
    use m_crystallography,only:trimr,make_g_vec_array
    use m_electron, only:element
	implicit none
    
    character(2),intent(in)::orbital
    integer*4,intent(in)::zz,lc
    real(fp_kind),intent(in)::DE
    real(fp_kind):: g_vec_array(3,nopiy,nopix)
    
    complex(fp_kind):: fz_mu(nopiy,nopix)
    
    !dummy variables
    integer(4) i,j
    real(fp_kind) sval

    real(fp_kind) svals(29),EELS_EDX_bscoef(4,29)
    !DATA POINTS USED FOR THE INTERPOLATION S-VALUES (q/2)
    data svals / 0.0_fp_kind,0.025_fp_kind,0.05_fp_kind,0.1_fp_kind,0.2_fp_kind,0.3_fp_kind,0.4_fp_kind,0.5_fp_kind,0.625_fp_kind,&
               & 0.75_fp_kind,0.875_fp_kind,1.0_fp_kind,1.5_fp_kind,2.0_fp_kind,2.5_fp_kind,3.0_fp_kind,3.5_fp_kind,4.0_fp_kind,  &
               & 5.0_fp_kind,6.0_fp_kind,7.0_fp_kind,8.0_fp_kind,9.0_fp_kind,10.0_fp_kind,12.0_fp_kind,14.0_fp_kind,16.0_fp_kind, &
               & 18.0_fp_kind,20.0_fp_kind /
    
    if (lc > 0) then
        write(*,*) 'Filling the custom ionization scattering factor grid for ' &
                    & //element(zz)//' '//orbital//', please wait...',char(10)
    else
        write(*,*) 'Filling the ionization scattering factor grid for ' &
                    & //element(zz)//' '//orbital//', please wait...',char(10)
    end if
    	
    ! 2025-05-22 JB: custom ionization data
    if ((lc>0) .AND. (.NOT.EDX)) then
        !read in the custom ionization form factors
        call load_custom_ionization_data(lc, EELS_EDX_bscoef(1,:))
    else
        !read the default EELS or EDX parameterization data file
        EELS_EDX_bscoef(1,:) = get_ionization_parameters(orbital,zz,DE,EDX)
    end if
    
    !pppack interpolation over s=q/2 to q-grid cubspl -> ppvalu
	call cubspl(svals,EELS_EDX_bscoef(:,:), 29, 0, 0 )
    
    fz_mu = 0.0_fp_kind
    call make_g_vec_array(g_vec_array,ifactory,ifactorx)
    !!$OMP PARALLEL PRIVATE(i, j, m2, m1, sky, skx, tempval, sval), SHARED(fz_mu) 
    !!$OMP DO
	do i=1, nopiy;do j=1, nopix
      sval = trimr(g_vec_array(:,i,j),ss) / 2.0_fp_kind
      if (sval.le.20.0_fp_kind) fz_mu(i,j) = cmplx(ppvalu(svals,EELS_EDX_bscoef(:,:),28,4,sval,0),0.0_fp_kind ,fp_kind) / (tp*ak1)
    enddo;enddo
    !!$OMP END DO
	!!$OMP END PARALLEL

    return
  end function

  !--------------------------------------------------------------------------------------
  !   make_mu_matrix() makes the mu matrices for each HOLZ slice
  !   subroutine to take the unit cell input, 
  !   and slice based on holz
  subroutine make_local_inelastic_potentials(ionization)
      
    use m_multislice
    use global_variables!, only:adf,nopiy,nopix,high_accuracy,nt,ss,ndet,ig1,ig2,ifactory
    use m_absorption
    !use m_string
    !use output
      
    implicit none
      
    logical,intent(in)::ionization
    integer(4)   i,j,k,nat_
    real(fp_kind) :: potential_matrix_complex(nopiy,nopix),vol
    real(8)::thmin,thmax,phmin,phmax
    complex(fp_kind)::fz_adf(nopiy,nopix,nt,ndet)

    write(*,134)        
134 format(/,' Calculating effective inelastic potentials.',/)

    if(allocated(adf_potential)) deallocate(adf_potential)
    if(allocated(ionization_potential)) deallocate(ionization_potential)

    allocate(adf_potential(nopiy,nopix,n_slices,ndet)) !the adf potential
    adf_potential=0
    if(ionization) allocate(ionization_potential(nopiy,nopix,num_ionizations,n_slices))  !the ionization potential
      
    !if(adf) call make_fz_adf()
	  if(adf.and.complex_absorption) then  
      do k=1,ndet/nseg
        thmin =  atan(inner((k-1)/nseg+1)/ak)
        thmax =  atan(outer((k-1)/nseg+1)/ak)
        !Note that the absorptive calculations do not take into account the directionality of inelastic scattering, the absorptive scattering
        !factors are assumed isotropic and this is only an approximation for inelastic scattering to segmented detectors
        fz_adf(:,:,:,(k-1)*nseg+1:k*nseg) = spread(absorptive_scattering_factors( &
                            & ig1,ig2,ifactory,ifactorx,nopiy,nopix,nt,a0,ss,atf,nat,&
                            & ak, relm, orthog,thmin,thmax),dim=4,ncopies=nseg)/nseg
      enddo
    endif

    do j = 1, n_slices
	    vol = ss_slice(7,j)
	    !calculate the ionization potential
	    if(ionization) then
		    do i=1,num_ionizations
                nat_ = nat_slice(atm_indices(i),j)
                ionization_potential(:,:,i,j) = real(potential_from_scattering_factors( &
                            & ionization_mu(:,:,i),tau_slice(:,atm_indices(i),:nat_,j), &
                            & nat_,nopiy,nopix,high_accuracy)/vol)
  			enddo
        endif  
	    !calculate the ADF potential
        if(adf.and.complex_absorption) then  
            do i=1,nt
                nat_ = nat_slice(i,j)
                do k=1,ndet
                    adf_potential(:,:,j,k) = adf_potential(:,:,j,k) + real( &
                            & potential_from_scattering_factors(fz_adf(:,:,i,k), &
                            & tau_slice(:,i,:nat_,j),nat_,nopiy,nopix,high_accuracy)/vol*ss(7)*4*pi)
                enddo
            enddo
        endif
    enddo	!ends loop over the number of potential subslices
      
  end subroutine
      
  !--------------------------------------------------------------------------------------
  function potential_from_scattering_factors(scattering_factor,atom_posn,nat_layer,nopiy,nopix,high_accuracy) result(slice_potential)
    use m_precision
    use cufft_wrapper
    use output
    
    implicit none
    
    integer(4),intent(in) :: nat_layer,nopiy,nopix
    real(fp_kind),intent(in) :: atom_posn(3,nat_layer)
    complex(fp_kind),intent(in)::scattering_factor(nopiy,nopix)
	logical,intent(in),optional::high_accuracy
    
    complex(fp_kind),dimension(nopiy,nopix) :: potential, site_term,slice_potential
	logical::high_accuracy_

    procedure(make_site_factor_generic),pointer :: make_site_factor
    
	high_accuracy_ = .false.;if(present(high_accuracy)) high_accuracy_= high_accuracy
#ifdef GPU
    make_site_factor => make_site_factor_cuda
#else
    make_site_factor => make_site_factor_matmul
#endif        
    
    slice_potential = 0.0_fp_kind
    
    if (nat_layer.ne.0) then
        if (high_accuracy_) then
            call make_site_factor(site_term, atom_posn)
        else
            call make_site_factor_hybrid(site_term, atom_posn)        
        endif
		!call ifft2(nopiy,nopix,site_term,nopiy,site_term,nopiy)
		!call binary_out(nopiy,nopix,site_term,'site_term')
		!stop
        slice_potential = site_term*scattering_factor
        ! Get realspace potential
        call ifft2(nopiy,nopix,slice_potential,slice_potential)
        slice_potential = slice_potential*sqrt(float(nopiy*nopix))
    endif
    
  end function
  
  
  function make_absorptive_grates(nopiy,nopix,n_slices) result(projected_potential)
    
    use m_precision, only: fp_kind
	use cufft_wrapper, only: fft2, ifft2
    use global_variables, only: ig1,ig2,ifactory,ifactorx,nt, relm, tp, ak, atf, high_accuracy, ci, pi, bwl_mat,fz,fz_DWF,ss,a0,nat,orthog
    use m_absorption!, only: transf_absorptive,fz_abs
    use m_multislice, only: nat_slice, ss_slice, tau_slice ! JB 2022-08-03 error fixed for ifort compile
    use m_string, only: to_string
    use output
        
    implicit none
        
    integer*4,intent(in)::nopiy,nopix,n_slices
    complex(fp_kind)::projected_potential(nopiy,nopix,n_slices)
        
    integer(4) :: j, m, n,nat_layer
    real(fp_kind) :: ccd_slice,V_corr
    complex(fp_kind),dimension(nopiy,nopix) :: scattering_pot,temp,effective_scat_fact
    complex(fp_kind)::fz_abs(nopiy,nopix,nt)
    
    real(fp_kind) :: t1, delta,amplitude(nopiy,nopix),phase(nopiy,nopix)
    
    procedure(make_site_factor_generic),pointer :: make_site_factor
    projected_potential= 0 
    t1 = secnds(0.0_fp_kind)
    fz_abs=0
    if(include_absorption) fz_abs = absorptive_scattering_factors(ig1,ig2,ifactory,ifactorx,nopiy,nopix,nt,a0,ss,atf,nat, ak, relm, orthog, 0.0_8, 4.0d0*atan(1.0d0))*2*ak
            
    do j = 1, n_slices
        write(*,'(1x, a, a, a, a, a)') 'Calculating transmission functions for slice ', to_string(j), '/', to_string(n_slices), '...'
        
198	    write(*,199) to_string(sum(nat_slice(:,j)))
199     format(1x, 'Number of atoms in this slice: ', a, /) 

		ccd_slice = relm / (tp * ak * ss_slice(7,j))
            V_corr = ss(7)/ss_slice(7,j)
        do m=1,nt
            nat_layer = nat_slice(m,j)
            effective_scat_fact = CCD_slice*fz(:,:,m)*fz_DWF(:,:,m)+cmplx(0,1)*fz_abs(:,:,m)*V_corr
            projected_potential(:,:,j) = projected_potential(:,:,j)+potential_from_scattering_factors(effective_scat_fact,tau_slice(:,m,:nat_layer,j),nat_layer,nopiy,nopix,high_accuracy)
        enddo
                
	enddo ! End loop over slices
	
	delta = secnds(t1)
        
	if(timing) then
        write(*,*) 'The calculation of transmission functions for the absorptive model took ', delta, 'seconds.'
        write(*,*)
		open(unit=9834, file=trim(adjustl(output_prefix))//'_timing.txt', access='append')
		write(9834, '(a, g, a, /)') 'The calculation of transmission functions for the absorptive model took ', delta, 'seconds.'
		close(9834)
    endif    
        
  end function make_absorptive_grates
  
  
  subroutine make_se_transmission_grates()
    
    use m_precision, only: fp_kind
	use cufft_wrapper, only: fft2
    use global_variables
    use m_multislice, only: n_slices, nat_slice, tau_slice, prop_distance, save_grates
    use m_string, only: to_string
    use output
        
    implicit none
            
    integer*4 :: i, j, m, n, nat_layer, ix, iy, jx, jy, nx2, ny2
    real(fp_kind) :: t1, delta, rad, val, x, y, x0, y0, rinv2
    real(fp_kind), allocatable :: layer_density(:,:)
    complex(fp_kind), allocatable :: atom_ff(:,:,:)
    character(32) :: snum
    
    t1 = secnds(0.0_fp_kind)
    if (ALLOCATED(se_transmission)) deallocate(se_transmission) ! deallocate previous transmission function if it exists
    allocate(se_transmission(nopiy,nopix,n_slices)) ! allocate new transmission function array
    allocate(layer_density(nopiy,nopix))
    allocate(atom_ff(nopiy,nopix,nt))
    
    se_transmission = 1.0_fp_kind ! initialize transmission function to 1.0
    layer_density = 0.0_fp_kind ! initialize layer density to zero
    nx2 = (nopix + modulo(nopix,2)) / 2 ! x nyquist pixel number 512 -> 256, 513 -> 257, 514 -> 257
    ny2 = (nopiy + modulo(nopiy,2)) / 2 ! y nyquist pixel number
    
    !write(*,*) '  total grid size in A: ', a0(1)*ifactorx, ' x ', a0(2)*ifactory
    ! prepare the atom form factor (real-space to fourier space)
    do i = 1, nt ! loop over atom types
        rad = se_transf_aty_radius(i) * 0.01 ! atomic radius of this atom type in A
        !write(*,*) '  preparing SE form factor for '//trim(adjustl(substance_atom_types(i)))// &
        !        &   ' with radius ', rad, ' A'
        rinv2 = -0.69_fp_kind / rad**2 ! exponential form parameter
        ! setup the atom form factor for this atom type in real space
        !$OMP PARALLEL DO PRIVATE(iy, ix, jy, jx, y, x) &
        !$OMP& SHARED(atom_ff, a0, ifactory, ifactorx, rinv2, nx2, ny2, nopix, nopiy, i) &
        !$OMP& COLLAPSE(2)
        do iy=1, nopiy; do ix=1,nopix
            jy = modulo(iy + ny2 - 1, nopiy) - ny2
            y = a0(2) * ifactory * jy / nopiy ! y coordinate in A
            jx = modulo(ix + nx2 - 1, nopix) - nx2
            x = a0(1) * ifactorx * jx / nopix ! x coordinate in A
            atom_ff(iy, ix, i) = EXP(rinv2 * (x**2 + y**2)) ! Gaussian form factor
        enddo; enddo
        !$OMP END PARALLEL DO
        call fft2(nopiy, nopix, atom_ff(:,:,i), atom_ff(:,:,i))
    enddo
    
    do j = 1, n_slices ! loop over slices of the supercell
        write(*,'(1x, a, a, a, a, a)') 'Calculating SE transmission function for slice ', &
                & to_string(j), '/', to_string(n_slices), '...'
	    write(*,301) to_string(sum(nat_slice(:,j)))
301     format(1x, 'Number of atoms in this slice: ', a, /) 
        layer_density = 0.0_fp_kind ! initialize layer density to zero
        do m=1,nt ! loop over atom types
            nat_layer = nat_slice(m,j) ! number of atoms of this type in the slice
            if (nat_layer == 0) cycle ! skip if no atoms of this type in the slice
            !write(*,*) '  adding', nat_layer, 'atoms of type ', trim(adjustl(substance_atom_types(m)))
            ! add the layer density contribution from this atom type
            layer_density(:,:) = layer_density(:,:) + real(potential_from_scattering_factors( &
                                & atom_ff(:,:,m),tau_slice(:,m,:nat_layer,j),nat_layer, &
                                & nopiy,nopix,high_accuracy))
        enddo
        layer_density = layer_density / sqrt(REAL(nopiy*nopix,fp_kind)) ! normalize layer density
        !write(*,*) '  max. layer density for slice ', j, ':', maxval(layer_density(:,:))
        ! debug storing of the layer density
        !write(unit=snum,fmt='(i0.4)') j
        !open(unit=3984, file=trim(adjustl(output_prefix))//"_SE-layer_"//trim(adjustl(snum))// &
        !                        & ".bin", form='binary', convert='big_endian')
        !write(3984) layer_density
        !close(3984)
        ! calculate SE transmission function from density, slice thickness and inelastic mean free path
        se_transmission(:,:,j) = exp(-1.0_fp_kind * layer_density * prop_distance(j) / se_imfp)
	enddo ! End loop over slices
	
	delta = secnds(t1)
        
	if(timing) then
        write(*,*) 'The calculation of SE transmission functions took ', delta, 'seconds.'
        write(*,*)
		open(unit=9834, file=trim(adjustl(output_prefix))//'_timing.txt', access='append')
		write(9834, '(a, g, a, /)') 'The calculation of SE transmission functions took ', delta, 'seconds.'
		close(9834)
    endif
    
    deallocate(layer_density, atom_ff) ! deallocate temporary arrays
    
    if (0 < IAND(arg_debug_dump,2)) then ! store the SE transmission functions
      write(*,*) 'Saving SE transmission functions to file...',char(10)
      open(unit=3984, file=trim(adjustl(output_prefix))//"_SE-transmission.bin", form='binary', convert='big_endian')
      write(3984) se_transmission
      close(3984)
    endif
        
  end subroutine make_se_transmission_grates
    
           
    integer function seed_rng() result(idum)
            
	    use m_numerical_tools, only: ran1
        use m_precision, only: fp_kind
            
        implicit none

        integer :: i
        real(fp_kind) :: random
            
	    idum = -1
            
	    do i = 1, nran
		    random = ran1(idum)
	    enddo
            
    end function seed_rng

      
    subroutine make_propagator(nopiy,nopix,prop,dz,ak1,ss,ig1,ig2,claue,ifactorx,ifactory,exponentiate)

        use m_precision, only: fp_kind
        use m_crystallography, only: trimr,make_g_vec_array
            
        implicit none

        integer(4) :: nopiy,nopix
        complex(fp_kind) :: prop(nopiy,nopix)        
        real(fp_kind) :: ak1, ss(7), claue(3), dz, g_vec_array(3,nopiy,nopix)
        integer(4) :: ifactorx, ifactory, ig1(3), ig2(3)
		logical,intent(in),optional::exponentiate
		
		logical::exp_
        real(fp_kind),parameter :: pi = atan(1.0d0)*4.0d0
        integer(4) :: ny, nx
        
		exp_ = .true.
		if(present(exponentiate)) exp_ = exponentiate
        
        call make_g_vec_array(g_vec_array,ifactory,ifactorx)

        do ny = 1, nopiy;do nx = 1, nopix
            ! calculate the propagator phase for each q and write it to the real part of prop            
            prop(ny,nx) = cmplx(-pi*dz*trimr(g_vec_array(:,ny,nx)-claue,ss)**2/ak1, 0.0d0, fp_kind )
        enddo;enddo

        ! exponentiate the propagator if requested, but calculate exp( i*phase )
		if(exp_) prop = exp( cmplx(0.0d0, 1.0d0, fp_kind) * prop)
        
        ! if exponentiate is not requested, prop is just the phase factor
        !        sitting in the real part of prop, imaginary part is zero

    end subroutine
	      
    subroutine make_propagator_components(nopiy,nopix,propy,propx,dz,ak1,ss,ig1,ig2,ifactorx,ifactory)

        use m_precision, only: fp_kind
        use m_crystallography, only: trimr,make_g_vec_array
            
        implicit none

        integer(4) :: nopiy,nopix
        complex(fp_kind) :: propy(nopiy),propx(nopix)
        real(fp_kind) :: ak1, ss(7), dz, g_vec_array(3,nopiy,nopix)
        integer(4) :: ifactorx, ifactory, ig1(3), ig2(3)

        real(fp_kind),parameter :: pi = atan(1.0d0)*4.0d0
        integer(4) :: ny, nx
        
        
        call make_g_vec_array(g_vec_array,ifactory,ifactorx)

        do ny = 1, nopiy;
            propy(ny) = exp(cmplx(0.0_fp_kind,-pi*dz*trimr(g_vec_array(:,ny,1),ss)**2/ak1))
        enddo

		do nx = 1, nopix
			propx(nx) = exp(cmplx(0.0_fp_kind,-pi*dz*trimr(g_vec_array(:,1,nx),ss)**2/ak1))
		enddo

    end subroutine
       
    function make_qep_grates(idum) result(projected_potential)
    
        use m_precision, only: fp_kind
	      use cufft_wrapper, only: fft2, ifft2
        use global_variables, only: nopiy, nopix, nt, relm, tp, ak, ak1, atf, high_accuracy, ci, pi, bwl_mat,fz
        use m_multislice
        use m_string, only: to_string
        use output, only: output_prefix,timing,binary_in
        use m_numerical_tools, only: displace

        implicit none
        
        integer(4),intent(inout) :: idum
    
        complex(fp_kind) :: projected_potential(nopiy,nopix,n_qep_grates,n_slices),temp(nopiy,nopix),scattering_pot(nopiy,nopix,nt)
        integer(4), allocatable :: handled(:,:)
		integer(4):: save_list(2,nt),match_count, i, j, m, n,ii,jj,jjj,kk,iii
        real(fp_kind) :: tau_holder(3),tau_holder2(3),ccd_slice,ums,amplitude(nopiy,nopix),phase(nopiy,nopix)
        real(fp_kind) :: mod_tau(3,nt,maxnat_slice,n_slices,n_qep_grates),t1, delta, msd
        logical::fracocc
    
        procedure(make_site_factor_generic),pointer :: make_site_factor
        
	    
        ! Search for fractional occupancy
        fracocc = any(atf(2,:).lt.0.99d0)
        
        t1 = secnds(0.0_fp_kind)

        do j = 1, n_slices
	        write(*,'(1x, a, a, a, a, a)') 'Calculating transmission functions for slice ', to_string(j), '/', to_string(n_slices), '...'
        
    	    write(*,199) to_string(sum(nat_slice(:,j)))
    199     format(1x, 'Number of atoms in this slice: ', a) 

            ccd_slice = relm / (tp * ak * ss_slice(7,j))

	        do i = 1, n_qep_grates
#ifdef GPU
    200	        format(a1, 1x, i3, '/', i3)
	            write(6,200, advance='no') achar(13), i, n_qep_grates
#else
    201         format(1h+, 2x, i3, '/', i3)
                write(6,201) i, n_qep_grates
#endif
                flush(6)
                
                ! Randomly displace the atoms
				if (.not.fracocc) then ! no fractional occ. displace all atoms
 	                do m = 1, nt
	                    do n = 1, nat_slice(m,j)
                            msd = atf(3,m) ! (default) init from structure model
                            !
			                call displace(tau_slice(1:3,m,n,j),mod_tau(1:3,m,n,j,i),sqrt(msd),a0_slice,idum)
	                    end do
                    end do
				else ! displace but handle shared sites
					allocate( handled(nt,maxnat_slice) )
					handled = 0
					do ii=1, nt
           
					    do jj = 1, nat_slice(ii,j)
						    if (handled(ii,jj).eq.1) cycle ! skip, shared site has been handled
						    tau_holder(1:3) = tau_slice(1:3,ii,jj,j)

						    save_list = 0
						    match_count = 0
                            msd = atf(3,ii) ! (default) init from structure model
                            !
						    ums = msd ! atf(3,ii)
						    do iii=ii+1,nt
                            do jjj=1,nat_slice(iii,j)
						 	if (same_site(tau_holder,tau_slice(1:3,iii,jjj,j))) then
						 	    match_count = match_count+1
							    save_list(1,match_count)=iii
							    save_list(2,match_count)=jjj
                                msd = atf(3,iii) ! (default) init from structure model
                                !
							    ums = ums + msd ! atf(3,iii)
							    cycle
                            end if
                            end do
                            end do

						    ums = ums / dfloat(match_count+1)
					        call displace(tau_holder(1:3),tau_holder2(1:3),sqrt(ums),a0_slice,idum)
						    mod_tau(1:3,ii,jj,j,i) = tau_holder2(1:3)
						    handled(ii,jj) = 1
						    do kk=1,match_count
							    mod_tau(1:3,save_list(1,kk),save_list(2,kk),j,i) = tau_holder2(1:3)
							    handled(save_list(1,kk),save_list(2,kk)) = 1
						    enddo
						   
					    enddo
					enddo

					deallocate( handled )
				endif
				
				projected_potential(:,:,i,j) = 0
				do m = 1, nt
					projected_potential(:,:,i,j) = projected_potential(:,:,i,j) + &
                        & real(potential_from_scattering_factors(CCD_slice*fz(:,:,m), &
                        & mod_tau(:,m,1:nat_slice(m,j),j,i),nat_slice(m,j),nopiy,nopix,high_accuracy))
                enddo
	        enddo ! End loop over grates
        
            write(*,*)
            write(*,*)
        
	    enddo ! End loop over slices
	
	    delta = secnds(t1)
        
        write(*,*) 'The calculation of transmission functions for the QEP model took ', delta, 'seconds.'
        write(*,*)

  	    if(timing) then
			open(unit=9834, file=trim(adjustl(output_prefix))//'_timing.txt', access='append')
			write(9834, '(a, g, a, /)') 'The calculation of transmission functions for the QEP model took ', delta, 'seconds.'
			close(9834)
        endif    
    
    end function make_qep_grates
    
    

	logical(4) function same_site(site1,site2)
      
        implicit none
      
        real(fp_kind) site1(3),site2(3)
        real(fp_kind) tol
      
        tol = 1.0d-6
        same_site = all(abs(site1-site2).lt.tol)
      
        return
    end function
    
    
end module
