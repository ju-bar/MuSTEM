!--------------------------------------------------------------------------------
!
!  Copyright (C) 2023  J. Barthel, L. J. Allen
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
!  Phonon density of states Monte-Carlo to determine a random mean squared
!  displacement for atom types in a QEP calculation.
!
!--------------------------------------------------------------------------------
!
!  How to use this module:
!
!  * call pd_scan_xtl(iunit, nt, substance_atom_types)
!    This scans the open file "iunit" for PDOS tables, initializes storage for
!    "nt" atom types and stores the tables under type indices determined from
!    matching type names in "substance_atom_types"
!  * call pd_report(substance_atom_types, atf) to print a report on loaded PDOS
!    tables
!  * call get_rand_msd(seed, iaty, msd) to update "msd" for atom type "iaty"
!    by a random "msd" sampled from the PDOS.
!
!--------------------------------------------------------------------------------
!
!  PDOS data structure (read from XTL files)
!
!  PDOS                (begin of PDOS data section)
!  Si                  (atom type name as in XTL structure definition)
!  28.0855             (atom mass in Daltons = 1 u)
!  300                 (temperature in K)
!  40                  (length of PDOS table)
!  0.002000   0.010406 (energy [eV], PDOS sample)
!  0.004000   0.193359 ( ... )
!  0.006000   0.772637
!  0.008000   1.884095
!  0.010000   3.575660    
!  ...
!
!--------------------------------------------------------------------------------
  
  
module pdos

	use m_precision, only: fp_kind
  use m_numerical_tools, only: r8_distr_rej, trapz

	! declaration of module members
	implicit none
	! - module constants
  integer(4),public,parameter :: pd_pdos_maxlen = 32768 ! max. length of the input PDOS
  real(8),public,parameter :: pd_ep_thr = 0.00001 ! phonon energy lower limit [eV]
  ! - module parameters (variables)
	integer(4),public :: pd_naty ! max. number of atom types to handle
  integer(4),public,allocatable :: pd_pdos_len(:) ! stores the length of each set of PDOS data (one set per atom type)
  real(fp_kind),public,allocatable :: pd_aty_mass(:) ! stores atomic mass in Dalton
  real(fp_kind),public,allocatable :: pd_aty_temp(:) ! stores temperature in K
  !
  ! 3 dimensions for storing the input PDOS
  ! dimension 1 : length 2 (energy, PDOS(energy))
  ! dimension 2 : length pd_pdos_maxlen (length of the PDOS table)
  ! dimension 3 : length pd_naty (number of atom types)
	real(8),public,allocatable :: pd_pdos_samples(:,:,:) ! samples of the phonon density of states
  
  contains
  
  
  !
  ! returns the msd for a given phonon energy, temperature and mass
  ! "ep" phonon energy [eV]
  ! "t" temperature [K]
  ! "m" mass [u]
  !
  function pd_msd(ep, t, m)
  
    implicit none
    
    real(fp_kind),parameter :: f1 = 0.0020900796_fp_kind ! prefactor hbar^2 / ( 2 e u ) [A^2 eV Dalton]
    real(fp_kind),parameter :: ft = 8.617333262145179e-05_fp_kind ! prefactor for thermal energy [eV K^(-1)]
    
    real(fp_kind), intent(in) :: ep, t, m
    real(fp_kind) :: pd_msd
    real(fp_kind) :: et
    
    et = t * ft ! thermal energy in eV
    !write(*,*) 'ep:', ep, '    et:',et
    pd_msd = f1 / (ep * m * tanh(0.5_fp_kind * ep / et)) ! MSD at temperature
    !write(*,*) 'pd_msd:',pd_msd
    return
    
  end function
  
  
  !
  ! initializes the module arrays for num_aty atom types
  !
  subroutine pd_init(num_aty)
    implicit none
    integer( kind=4 ), intent(in) :: num_aty
    integer( kind=4 ) :: ialloc
    if (allocated(pd_pdos_len)) deallocate(pd_pdos_len, stat=ialloc)
    if (allocated(pd_aty_mass)) deallocate(pd_aty_mass, stat=ialloc)
    if (allocated(pd_aty_temp)) deallocate(pd_aty_temp, stat=ialloc)
    if (allocated(pd_pdos_samples)) deallocate(pd_pdos_samples, stat=ialloc)
    pd_naty = num_aty
    if (pd_naty > 0) then
      allocate(pd_pdos_len(pd_naty), stat=ialloc)
      pd_pdos_len = 0
      allocate(pd_aty_mass(pd_naty), stat=ialloc)
      allocate(pd_aty_temp(pd_naty), stat=ialloc)
      allocate(pd_pdos_samples(2, pd_pdos_maxlen, pd_naty), stat=ialloc)
    end if
  end subroutine pd_init
  
  
  !
  ! scans the file unit "iunit" for pdos table input and
  ! stores it in module tables for "nt" atom types
  ! identified by the strings "substance_atom_types"
  !
  subroutine pd_scan_xtl(iunit, nt, substance_atom_types)
    implicit none
    integer( kind=4 ), intent(in) :: iunit, nt
    character*10, intent(in) :: substance_atom_types(nt)
    integer( kind=4 ) :: ioerr, iaty, i
    character*1024 :: sdummy, stype
    
    call pd_init(nt) ! initialize PDOS module for this XTL file
    
100 continue
    sdummy = ""
    do while (trim(adjustl(sdummy)) /= "PDOS" .and. trim(adjustl(sdummy)) /= "pdos")
      read(iunit,*,iostat=ioerr) sdummy
      if (ioerr /= 0) goto 900 ! error or eof, get out
    end do
    
    write(6,*)
    write(6,*) '  Found PDOS table in the XTL file.'
    
    ! determine atom type index in substance_atom_types
    read(iunit,*,iostat=ioerr) stype
    if (ioerr /= 0) goto 900 ! error or eof, get out
    iaty = 0
    do i=1, nt
      if (0 < index(substance_atom_types(i), trim(adjustl(stype)))) then
        iaty = i
        exit
      end if
    end do
    if (iaty == 0) then
      write(unit=6,fmt='(A,I5)') '     Error: atom type "'//trim(adjustl(stype))//'" not in structure.'
      write(6,*) '    Table ignored.'
      goto 100 ! atom type not identified, try next PDOS table
    end if
    
    write(6,*) '  * atom type: '//trim(substance_atom_types(iaty))
    
    ! index iaty is the index of the atom type
    
    ! read the mass of the atom type
    read(iunit,*,iostat=ioerr) pd_aty_mass(iaty)
    if (ioerr /= 0) goto 100 ! atom mass could not be read, try next PDOS table
    write(unit=sdummy,fmt='(F12.5)') pd_aty_mass(iaty)
    write(6,*) '  * atomic mass: '//trim(adjustl(sdummy))//' u'
    
    ! read the temperature
    read(iunit,*,iostat=ioerr) pd_aty_temp(iaty)
    if (ioerr /= 0) goto 100 ! temperature could not be read, try next PDOS table
    write(unit=sdummy,fmt='(F6.1)') pd_aty_temp(iaty)
    write(6,*) '  * temperature: '//trim(adjustl(sdummy))//' K'
    
    ! read length of the PDOS table
    read(iunit,*,iostat=ioerr) pd_pdos_len(iaty)
    if (ioerr /= 0) goto 100 ! length of PDOS could not be read, try next PDOS table
    write(unit=sdummy,fmt=*) pd_pdos_len(iaty)
    write(6,*) '  * table of length: '//trim(adjustl(sdummy))
    
    if (pd_pdos_len(iaty) > pd_pdos_maxlen) then
      write(unit=6,fmt='(A,I5)') '   Error: unsupported length, max. allowed: ', pd_pdos_maxlen
      if (ioerr /= 0) goto 100 ! length of PDOS not supported, try next PDOS table
    end if
    
    if (pd_pdos_len(iaty)==0) then
      write(6,*) '    Zero table length means working with the original MSD.'
      goto 100 ! table of zero length, try reading another table
    end if
    
    ! read the PDOS table
    do i=1, pd_pdos_len(iaty)
      read(iunit,*,iostat=ioerr) pd_pdos_samples(1, i, iaty), pd_pdos_samples(2, i, iaty)
      if (ioerr /= 0) goto 100 ! PDOS sample could not be read, try next PDOS table
    end do
    
    write(6,*) '    Finished reading the PDOS.'
    
    goto 100 ! try reading another table
  
    
900 continue
    return
    
  end subroutine pd_scan_xtl
  
  
  !
  ! prints a report on which msd model will be used for
  ! which atom type
  !
  subroutine pd_report(nt, substance_atom_types, atf)
    
    implicit none
    
    integer( kind=4 ), intent(in) :: nt
    real(fp_kind), intent(in) :: atf(3,nt)
    character*10, intent(in) :: substance_atom_types(nt)
    
    integer(4) :: i, j, n, ialloc
    real(fp_kind) :: msd, ep, msdj, pdosj, wpdos
    character*10 :: smodel, stype
    real(fp_kind), allocatable :: p_msd(:,:)
    
    write(6,*)
    write(6,*) '  Thermal vibration parameters'
    if (nt /= pd_naty) then
      write(6,*) '  Error: inconsistency in number of atom types.'
      write(6,fmt='(A,I)') '    XTL structure data: ', nt
      write(6,fmt='(A,I)') '    number of PDOS tables: ', pd_naty
      stop 5
    end if
    
    write(6,101)
    write(6,102)
    
    if (nt > 0) then
      do i=1, nt ! loop over atom types
        stype = substance_atom_types(i)
        smodel = 'Einstein'
        msd = atf(3,i)
        n = pd_pdos_len(i) ! length of the PDOS table for this atom type
        if (n == 0) then ! MSD will be taken from the structure file
          write(6,103) i, substance_atom_types(i), smodel, atf(3,i)  
        end if
        if (n == 1) then ! Einstein model but modified msd
          ep = real( pd_pdos_samples(1, 1, i), kind=fp_kind ) ! single energy value -> Einstein model
          msd = pd_msd(ep, pd_aty_temp(i), pd_aty_mass(i))
          write(6,104) i, substance_atom_types(i), smodel, atf(3,i), msd  
        end if
        if (n > 1) then ! PDOS model
          smodel = 'PDOS'
          ! calculate avg. MSD from PDOS
          if (allocated(p_msd)) deallocate(p_msd)
          allocate(p_msd(2,n), stat=ialloc)
          ! calculate the weight of the PDOS first
          do j=1, n ! loop over PDOS points
            ! store phonon energy in eV
            p_msd(1, j) = real( pd_pdos_samples(1, j, i), kind=fp_kind )
            ! store PDOS
            p_msd(2, j) = real( pd_pdos_samples(2, j, i), kind=fp_kind )
          end do
          wpdos = trapz(n, p_msd) ! total weight of the PDOS for normalization
          do j=1, n ! loop over PDOS points
            ! store phonon energy in eV
            p_msd(1, j) = real( pd_pdos_samples(1, j, i), kind=fp_kind )
            ! get (harmonic oscillator) MSD at temperature for each phonon energy
            msdj = pd_msd(p_msd(1, j), pd_aty_temp(i), pd_aty_mass(i))
            ! get normalized value of the PDOS
            pdosj = real( pd_pdos_samples(2, j, i), kind=fp_kind ) / wpdos
            ! store MSD at energy weighted with PDOS
            p_msd(2, j) = msdj * pdosj
          end do
          msd = trapz(n, p_msd) ! get avg. msd by integrating over phonon energies
          write(6,104) i, substance_atom_types(i), smodel, atf(3,i), msd, pd_aty_temp(i)
        end if
      end do
      
    end if

101 format('   Type Name       Model       Structure      PDOS           Temperature')
102 format('                               <u>**2 [A**2]  <u>**2 [A**2]  T [K]')
103 format('  ',i2,4x,a10,1x,a10,1x,g12.5,'    -              -')
104 format('  ',i2,4x,a10,1x,a10,1x,g12.5,3x,g12.5,3x,f6.1)

    return  
  end subroutine pd_report
  
  !
  ! updates the "msd" value by a random msd based on the
  ! pdos stored for atom type index "iaty"
  !
  ! "seed" is a running rng seed number
  !
  subroutine get_rand_msd(seed, iaty, msd)
  
    implicit none
  
    integer(4), parameter :: ntry_max = 100
    integer(4), intent(inout) :: seed
    integer(4), intent(in) :: iaty
    real(fp_kind), intent(inout) :: msd
    
    integer(4) :: n, ntry
    real(8) :: ran_ep ! random phonon energy
    real(fp_kind) :: ep
    
    !write(*,*) 'msd in :' , msd
    ntry = 0
    if (iaty > 0 .and. iaty <= pd_naty) then ! valid atom type index
      n = pd_pdos_len(iaty) ! table length
      if (n > 0) then ! valid length of pdos
        ! get random phonon energy (assuming PDOS input in eV)
        ran_ep = 0.0D+0
        do while (ran_ep <= pd_ep_thr) ! lower limit of phonon energy
          if (ntry > ntry_max) return ! infinity loop catch
          ntry = ntry + 1
          ran_ep = r8_distr_rej( seed, n, pd_pdos_samples(1:2, 1:n, iaty) ) ! random phonon energy from PDOS [eV]
          !write(*,*) 'ran_ep :' , ran_ep
        end do
        ep = real( ran_ep, kind=fp_kind )
        !write(*,*) 'ep :' , ep
        ! calculate MSD from ran_ep in [A^2] at temperature pd_aty_temp(iaty)
        ! and for mass pd_aty_mass(iaty)
        msd = pd_msd(ep, pd_aty_temp(iaty), pd_aty_mass(iaty))
      end if
    end if
    !write(*,*) 'msd out:' , msd
    return
  end subroutine get_rand_msd
  
  
  
end module pdos
  