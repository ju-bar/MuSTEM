!--------------------------------------------------------------------------------
!
!  Copyright (C) 2025  L. J. Allen, H. G. Brown, A. J. D’Alfonso, S.D. Findlay, B. D. Forbes, J. Barthel
!
!  modified:
!    2019-Dec-13, J. Barthel, added plasmon monte-carlo triggers and calls
!    2025-Jun-02, J. Barthel, replaced get_sum using manual reduction kernels
!    
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

module cuda_ms
     
    use global_variables
	use plasmon
    use cuda_array_library
    use cudafor
    
    implicit none
    
    integer, parameter :: max_blocks = 1048576
    real(fp_kind), device :: block_sums_d(max_blocks)
    real(fp_kind) :: block_sums(max_blocks)

    interface get_sum
        module procedure get_sum_complex
        module procedure get_sum_real
    end interface

    interface cuda_stem_detector
        module procedure cuda_stem_detector_wavefunction
        module procedure cuda_stem_detector_cbed
    end interface
    
    interface cuda_multislice_iteration
        module procedure cuda_multislice_iteration_fractorized_prop
        module procedure cuda_multislice_iteration_single_prop
        module procedure cuda_multislice_iteration_new
        module procedure cuda_multislice_iteration_cshifted_transf
        end interface
    
    contains

	attributes(host) subroutine cuda_phase_shift_from_1d_factor_arrays(arrayin,arrayout,shift_arrayy_d,shift_arrayx_d,nopiy,nopix,plan)
		use CUFFT_wrapper
		use cuda_array_library, only: blocks, threads
		complex(fp_kind),intent(in),device::arrayin(nopiy,nopix),shift_arrayy_d(nopiy),shift_arrayx_d(nopix)
		integer(4),intent(in)::nopiy,nopix
		integer,value::plan

		complex(fp_kind),intent(out),device::arrayout(nopiy,nopix)

		complex(fp_kind),device::shift_array_d(nopiy,nopix)

		integer*4::shifty,shiftx		

			call cuda_make_shift_array<<<blocks,threads>>>(shift_array_d,shift_arrayy_d,shift_arrayx_d,nopiy,nopix)     !make the qspace shift array
            call cuda_multiplication<<<blocks,threads>>>(arrayin,shift_array_d, arrayout,1.0_fp_kind,nopiy,nopix) !multiply by the qspace shift array
            call cufftExec(plan,arrayout,arrayout,CUFFT_INVERSE)!inverse fourier transform
	end subroutine

	attributes(host) subroutine cuda_image(psi_d,ctf_d,image_d,normalisation, nopiy, nopix,plan,rspacein)
		use CUFFT_wrapper
        use cuda_array_library, only: blocks, threads
		implicit none

		complex(fp_kind),intent(in),dimension(nopiy,nopix),device::psi_d,ctf_d
		real(fp_kind),intent(out),dimension(nopiy,nopix),device::image_d
		integer(4),intent(in)::nopiy,nopix
		real(fp_kind),intent(in),value::normalisation
		integer,value::plan
		logical,intent(in)::rspacein

		complex(fp_kind),dimension(nopiy,nopix),device::psi_temp_d

		if(rspacein)  then 
			call cufftExec(plan, psi_d, psi_temp_d, CUFFT_FORWARD)
		else
			psi_temp_d = psi_d
		endif
		call cuda_multiplication<<<blocks,threads>>>(psi_temp_d, ctf_d, psi_temp_d, sqrt(normalisation), nopiy, nopix)
        call cufftExec(plan, psi_temp_d, psi_temp_d, CUFFT_INVERSE)
		call cuda_mod<<<blocks,threads>>>(psi_temp_d, image_d, normalisation, nopiy, nopix)
	end subroutine

	attributes(host) subroutine cuda_multislice_iteration_cshifted_transf(psi_d, transf_d, prop_d, normalisation, nopiy, nopix,nucy,nucx,plan)
		use CUFFT_wrapper
		use cuda_array_library, only: blocks, threads
		! (taken out because unused, JB-191213) use output
		implicit none

		complex(fp_kind),device,dimension(nopiy,nopix)::psi_d, transf_d, prop_d, psi_out_d
		integer(4):: nopiy,nopix,nucy,nucx
		integer(4):: iexc,isx,isy ! local plasmon scatter results added here, JB-191213
		real(fp_kind),intent(in),value::normalisation
		! (taken out because unused, JB-191213) complex(fp_kind)::out(nopiy,nopix)
		integer,value::plan
			
		call cuda_multiplication<<<blocks,threads>>>(psi_d,transf_d, psi_out_d,1.0_fp_kind,nopiy,nopix,nucy,nucx)
		call cufftExec(plan,psi_out_d,psi_d,CUFFT_FORWARD)
		! plasmon trigger, JB-191213
		if (plasmonmc) then
			call pl_scatt_slc(pl_tslc,pl_grid_sqx,pl_grid_sqx,iexc,isx,isy,pl_idum)
			if (iexc>0 .and. (isx/=0 .or. isy/=0)) then ! scattering happened and shift is needed
				psi_out_d = psi_d ! extra copy required to avoid access conflicts on shifting and to remain consistent with the sequence of using psi_d and psi_out_d
				call cuda_cshift<<<blocks,threads>>>(psi_out_d,psi_d,nopiy,nopix,isy,isx) ! shift from psi_out_d to psi_d by pixel distances (isy,isx) on device
			endif
		endif
		call cuda_multiplication<<<blocks,threads>>>(psi_d,prop_d, psi_out_d,normalisation,nopiy,nopix)
		call cufftExec(plan,psi_out_d,psi_d,CUFFT_INVERSE)

	end subroutine

	attributes(host) subroutine cuda_multislice_iteration_single_prop(psi_d, transf_d, prop_d, normalisation, nopiy, nopix,plan)
		use CUFFT_wrapper
		use cuda_array_library, only: blocks, threads
		implicit none

		complex(fp_kind),device,dimension(nopiy,nopix)::psi_d, transf_d, prop_d, psi_out_d
		integer(4):: nopiy,nopix
		real(fp_kind),intent(in),value::normalisation
		integer,value::plan
		integer(4):: iexc,isx,isy ! local plasmon scatter results added here, JB-191213

		call cuda_multiplication<<<blocks,threads>>>(psi_d,transf_d, psi_out_d,1.0_fp_kind,nopiy,nopix)
		call cufftExec(plan,psi_out_d,psi_d,CUFFT_FORWARD)
		! plasmon trigger, JB-191213
		if (plasmonmc) then
			call pl_scatt_slc(pl_tslc,pl_grid_sqx,pl_grid_sqx,iexc,isx,isy,pl_idum)
			if (iexc>0 .and. (isx/=0 .or. isy/=0)) then ! scattering happened and shift is needed
				psi_out_d = psi_d ! extra copy required to avoid access conflicts on shifting and to remain consistent with the sequence of using psi_d and psi_out_d
				call cuda_cshift<<<blocks,threads>>>(psi_out_d,psi_d,nopiy,nopix,isy,isx) ! shift from psi_out_d to psi_d by pixel distances (isy,isx) on device
			endif
		endif
		call cuda_multiplication<<<blocks,threads>>>(psi_d,prop_d, psi_out_d,normalisation,nopiy,nopix)
		call cufftExec(plan,psi_out_d,psi_d,CUFFT_INVERSE)

	end subroutine

	attributes(host) subroutine cuda_multislice_iteration_fractorized_prop(psi_d, transf_d, propy_d, propx_d,normalisation, nopiy, nopix,plan)
		use CUFFT_wrapper
		use cuda_array_library, only: blocks, threads
		implicit none

		complex(fp_kind),device,dimension(nopiy,nopix)::psi_d, transf_d, psi_out_d
		complex(fp_kind),device::propy_d(nopiy),propx_d(nopix)
		integer(4):: nopiy,nopix
		real(fp_kind),intent(in),value::normalisation
		integer,value::plan
		integer(4):: iexc,isx,isy ! local plasmon scatter results added here, JB-191213

		call cuda_multiplication<<<blocks,threads>>>(psi_d,transf_d, psi_d,1.0_fp_kind,nopiy,nopix)
		call cufftExec(plan,psi_d,psi_d,CUFFT_FORWARD)
		! plasmon trigger, JB-191213
		if (plasmonmc) then
			call pl_scatt_slc(pl_tslc,pl_grid_sqx,pl_grid_sqx,iexc,isx,isy,pl_idum)
			if (iexc>0 .and. (isx/=0 .or. isy/=0)) then ! scattering happened and shift is needed
				psi_out_d = psi_d ! extra copy required to avoid access conflicts on shifting and to remain consistent with the sequence of using psi_d and psi_out_d
				call cuda_cshift<<<blocks,threads>>>(psi_out_d,psi_d,nopiy,nopix,isy,isx) ! shift from psi_out_d to psi_d by pixel distances (isy,isx) on device
			endif
		endif
		call cuda_multiplication<<<blocks,threads>>>(psi_d,propy_d,propx_d, psi_d,normalisation,nopiy,nopix)
		call cufftExec(plan,psi_d,psi_d,CUFFT_INVERSE)

	end subroutine

	attributes(host) subroutine cuda_multislice_iteration_new(psi_d, transf_d, prop_d,prop_distance,even_slicing, normalisation, nopiy, nopix,plan)
		use CUFFT_wrapper
		use cuda_array_library, only: blocks, threads
		implicit none

		complex(fp_kind),device,dimension(nopiy,nopix)::psi_d, transf_d, prop_d, psi_out_d
		complex(fp_kind),device,allocatable:: propp_d(:,:)
		integer(4):: nopiy,nopix
		real(fp_kind),intent(in),value::normalisation,prop_distance
		logical,intent(in)::even_slicing
		integer,value::plan
		integer(4):: iexc,isx,isy ! local plasmon scatter results added here, JB-191213
			
		call cuda_multiplication<<<blocks,threads>>>(psi_d,transf_d, psi_out_d,1.0_fp_kind,nopiy,nopix)
		call cufftExec(plan,psi_out_d,psi_d,CUFFT_FORWARD)
		! plasmon trigger, JB-191213
		if (plasmonmc) then
			call pl_scatt_slc(pl_tslc,pl_grid_sqx,pl_grid_sqx,iexc,isx,isy,pl_idum)
			if (iexc>0 .and. (isx/=0 .or. isy/=0)) then ! scattering happened and shift is needed
				psi_out_d = psi_d ! extra copy required to avoid access conflicts on shifting and to remain consistent with the sequence of using psi_d and psi_out_d
				call cuda_cshift<<<blocks,threads>>>(psi_out_d,psi_d,nopiy,nopix,isy,isx) ! shift from psi_out_d to psi_d by pixel distances (isy,isx) on device
			endif
		endif
		if(even_slicing) then
			call cuda_multiplication<<<blocks,threads>>>(psi_d,prop_d, psi_out_d,normalisation,nopiy,nopix)
		else
			allocate(propp_d(nopiy,nopix))
			call cuda_complex_exponentiation<<<blocks,threads>>>(propp_d, prop_d,prop_distance, nopiy, nopix)
			call cuda_multiplication<<<blocks,threads>>>(psi_d,propp_d, psi_out_d,normalisation,nopiy,nopix)
		endif
		call cufftExec(plan,psi_out_d,psi_d,CUFFT_INVERSE)

	end subroutine
	 
	
    attributes(global) subroutine cuda_phase_shift(shift_array, ifactory, ifactorx, ig1, ig2, n, m, coord)
		
        implicit none
    
        integer(4), value :: n,m,ix,iy
        integer(4), value :: ifactory,ifactorx
        integer(4) :: m1,m2
        integer(4),device,dimension(3) :: ig1,ig2,kx,ky
        real(fp_kind),device,dimension(3) :: kr
        real(fp_kind),device,dimension(3) :: coord
        real(fp_kind), value :: tp

        complex(fp_kind),dimension(n,m) :: shift_array

        tp = atan(1.0_fp_kind)*8.0_fp_kind
        
        !calculate the number of threads and blocks
        ix = (blockIdx%x-1)*blockDim%x + threadIdx%x
        iy = (blockIdx%y-1)*blockDim%y + threadIdx%y

        if(ix <= n) then
            if (iy <= m) then
                m1 = mod(ix+(n-1)/2-1,n)-(n-1)/2
                ky = m1*ig1
                m2 = mod(iy+(m-1)/2-1,m)-(m-1)/2           
                kx = m2*ig2
                kr = kx/float(ifactorx)+ky/float(ifactory)
                shift_array(ix,iy) = exp(cmplx(0.0_fp_kind,-tp*dot_product( kr, coord) ))
            endif
        endif   

    end subroutine cuda_phase_shift

	attributes(host) subroutine cuda_accumulate_intensity(psi,intensity, nopiy, nopix)
		use cuda_array_library, only: blocks, threads
		complex(fp_kind),device,dimension(nopiy,nopix)::psi
		real(fp_kind),device,dimension(nopiy,nopix)::intensity,temp
		integer(4):: nopiy, nopix
		call cuda_mod<<<blocks,threads>>>(psi, temp, 1.0_fp_kind, nopiy, nopix)
        call cuda_addition<<<blocks,threads>>>(intensity, temp, intensity, 1.0_fp_kind, nopiy, nopix)                  
	end subroutine
   !
	!attributes(host) subroutine cuda_multislice_iteration(psi_d,transf_d,prop_d,plan)
	!	use CUFFT_wrapper
   !    use cuda_array_library, only: blocks, threads
   !
	!	implicit none
   !
	!	integer,value::plan
	!	complex(fp_kind),device,dimension(nopiy,nopix)::psi_d,transf_d,prop_d
	!	complex(fp_kind),device::psi_out_d(nopiy,nopix)
   !
	!	! Transmission
	!	call cuda_multiplication<<<blocks,threads>>>(psi_d, transf_d, psi_out_d, 1.0_fp_kind, nopiy, nopix)
   !            
	!	! Propagate
	!	call cufftExec(plan, psi_out_d, psi_d, CUFFT_FORWARD)
	!	call cuda_multiplication<<<blocks,threads>>>(psi_d, prop_d, psi_out_d, normalisation, nopiy, nopix)
	!	call cufftExec(plan, psi_out_d, psi_d, CUFFT_INVERSE)
	!end subroutine
    
    attributes(host) function cuda_stem_detector_wavefunction(psi_d, mask_d)
    
        use cuda_array_library, only: blocks, threads
    
        implicit none
    
        real(fp_kind) :: cuda_stem_detector_wavefunction
        complex(fp_kind), device, dimension(nopiy,nopix),intent(in) :: psi_d
        real(fp_kind), device, dimension(nopiy,nopix) :: mask_d, cbed_d

        call cuda_mod<<<blocks,threads>>>(psi_d,cbed_d,normalisation,nopiy,nopix)
        call cuda_multiplication<<<blocks,threads>>>(cbed_d,mask_d, cbed_d ,1.0_fp_kind,nopiy,nopix)

        cuda_stem_detector_wavefunction = get_sum(cbed_d)
        
    end function cuda_stem_detector_wavefunction
    
	!attributes(global) subroutine cuda_stem_detector_wavefunction_on_the_fly(psi_d,result_d, inner,outer,dimy,dimx,n,m) 
  !
  !    implicit none
  ! 
  !    integer*4:: l
	!	real(fp_kind),device :: get_sum_complex_d,xx_d,yy_d,dimx_d,dimy_d,r_d,inner_d,outer_d
	!	real(fp_kind),value:: result_d,inner,outer,dimy,dimx,xx,yy,r
	!	complex(fp_kind),device::a
  !    integer :: istat
	!	integer,value:: n,m
	!	complex(fp_kind), device, dimension(n,m) :: psi_d
  !    
	!	
	!	integer(4), value :: ix,iy
	!	real(fp_kind),value :: iix,iiy
  !
  !    ix = (blockIdx%y-1)*blockDim%y + threadIdx%y
  !    iy = (blockIdx%x-1)*blockDim%x + threadIdx%x
  !
  !    if (ix <= m) then
  !         if (iy <= n) then
  !            !iix = abs(modulo(ix+m/2,m)-m/2)/(m/2.0_fp_kind)
	!			!iiy = abs(modulo(iy+n/2,n)-n/2)/(n/2.0_fp_kind)
	!			!iix = sqrt(iix**2+iiy**2)
	!			!if(.not.(iix<2.0_fp_kind/3)) then
	!			!	arrayin(iy,ix) = 0
	!			!else
	!			!	arrayin(iy,ix) = arrayin(iy,ix)*scale
	!			!endif
	!			
	!			xx = (modulo(ix+m/2-1,m)-m/2+1)/dimx
	!			yy = (modulo(iy+n/2-1,n)-n/2+1)/dimy
	!			r = sqrt(xx**2+yy**2)
	!			!psi_d(iy,ix) = r
  !            if(r<outer.and.r>inner) then
	!				!a = psi_d(iy,ix)
	!				r = abs(a)**2
	!				istat = atomicAdd(result_d,r)
	!			endif
  !        endif
  !    endif   
  !end subroutine cuda_stem_detector_wavefunction_on_the_fly

    !attributes(host) function cuda_stem_detector_wavefunction_on_the_fly(psi_d, inner,outer,dimy,dimx) 
    !
	!	use cudafor
    !    use cuda_array_library, only: blocks, threads
    !
    !    implicit none
	!	
	!	integer*4:: l,m
    !    real(fp_kind) :: cuda_stem_detector_wavefunction_on_the_fly
	!	real(fp_kind),device :: get_sum_complex_d,xx_d,yy_d,dimx_d,dimy_d,r_d,inner_d,outer_d
	!	real(fp_kind):: inner,outer,dimy,dimx,xx,yy,r
    !    complex(fp_kind), device, dimension(nopiy,nopix),intent(in) :: psi_d
    !    integer :: istat
	!
    !   get_sum_complex_d = 0.0_fp_kind
	!	write(*,*) 'hi'
	!	!write(*,*) inner;inner_d=inner
	!	!write(*,*) outer;outer_d=outer
	!	!write(*,*) dimx;dimx_d=dimx
	!	!write(*,*) dimy;dimy_d=dimy
	!	
	!	write(*,*) 'hi'
    !    !$cuf kernel do (2) <<<*,*>>>
    !    do m = 1, nopix
    !        do l = 1, nopiy
	!			!xx_d = 
	!			!xx = (modulo(m+nopix/2-1,m)-m/2+1)/dimx
	!			!yy = (modulo(l+nopiy/2-1,l)-l/2+1)/dimy
	!			!r = xx**2+yy**2
    !            !if(r<outer.and.r_d>inner) get_sum_complex_d = get_sum_complex_d + abs(psi_d(l,m))**2
	!			istat = atomicAdd(get_sum_complex_d, abs(psi_d(l,m))**2)
	!			!get_sum_complex_d = get_sum_complex_d + abs(psi_d(l,m))**2
    !        enddo
    !    enddo
    !    
    !    cuda_stem_detector_wavefunction_on_the_fly = get_sum_complex_d
    !    
    !end function cuda_stem_detector_wavefunction_on_the_fly


    
    attributes(host) function cuda_stem_detector_cbed(cbed_d, mask_d)
    
        use cuda_array_library, only: blocks, threads
        
        implicit none
    
        real(fp_kind) :: cuda_stem_detector_cbed
        real(fp_kind), device, dimension(nopiy,nopix) :: mask_d, image_d, cbed_d

        call cuda_multiplication<<<blocks,threads>>>(cbed_d,mask_d, image_d ,1.0_fp_kind,nopiy,nopix)
    
        cuda_stem_detector_cbed = get_sum(image_d)
    
    end function cuda_stem_detector_cbed
    
    
    ! -------------------------------
    ! Step 1: Per-block reduction
    ! -------------------------------
    attributes(global) subroutine reduce_real_kernel(A, n, block_sums)
      implicit none
      integer, value :: n
      real(fp_kind), device :: A(n)
      real(fp_kind), device :: block_sums(*)
      real(fp_kind), shared :: cache(256)
      integer :: tid, i, idx, gid
      real(fp_kind) :: temp

      tid = threadIdx%x
      gid = blockIdx%x
      idx = (gid - 1) * blockDim%x + tid ! 1-based indexing

      if (idx <= n) then
         temp = A(idx)
      else
         temp = 0.0_fp_kind
      end if
      cache(tid) = temp

      call syncthreads()

      i = blockDim%x / 2
      do while (i > 0)
         if (tid <= i) then
            cache(tid) = cache(tid) + cache(tid + i)
         end if
         call syncthreads()
         i = i / 2
      end do

      if (tid == 1) then
         block_sums(gid) = cache(1)
      end if
    end subroutine reduce_real_kernel
    
    attributes(global) subroutine reduce_complex_kernel(A, n, block_sums)
      implicit none
      integer, value :: n
      complex(fp_kind), device :: A(n)
      real(fp_kind), device :: block_sums(*)
      real(fp_kind), shared :: cache(256)
      integer :: tid, i, idx, gid
      complex(fp_kind) :: temp
      real(fp_kind) :: mod_sqr

      tid = threadIdx%x
      gid = blockIdx%x
      idx = (gid - 1) * blockDim%x + tid ! 1-based indexing

      if (idx <= n) then
         temp = A(idx)
         mod_sqr = real(temp)**2 + aimag(temp)**2
      else
         mod_sqr = 0.0_fp_kind
      end if
      cache(tid) = mod_sqr

      call syncthreads()

      i = blockDim%x / 2
      do while (i > 0)
         if (tid <= i) then
            cache(tid) = cache(tid) + cache(tid + i)
         end if
         call syncthreads()
         i = i / 2
      end do

      if (tid == 1) then
         block_sums(gid) = cache(1)
      end if
    end subroutine reduce_complex_kernel
    

    !rewritten to work with manual reduction kernel
    ! attributes(host) function get_sum_complex(psi_d)
    
        ! implicit none
    
        ! integer(4) l,m
        ! real(fp_kind) :: get_sum_complex
        ! real(fp_kind),device :: get_sum_complex_d
        ! complex(fp_kind), device, dimension(nopiy,nopix) :: psi_d
    
        ! get_sum_complex_d = 0.0_fp_kind
    
        ! !$cuf kernel do (2) <<<*,*>>>
        ! do m = 1, nopix
            ! do l = 1, nopiy
                ! get_sum_complex_d = get_sum_complex_d + abs(psi_d(l,m))**2
            ! enddo
        ! enddo
        
        ! get_sum_complex = get_sum_complex_d
    
    ! end function get_sum_complex
    
    !new sum for complex taking mod square on the run, 2025-06-02 JB
    attributes(host) function get_sum_complex(psi_d)
    
        implicit none
        
        integer, parameter :: threads = 256
        complex(fp_kind), device, dimension(nopiy,nopix) :: psi_d
        real(fp_kind) :: get_sum_complex
        integer :: n, blocks

        !init numbers
        n = nopiy * nopix
        blocks = (n + threads - 1) / threads
        if (blocks > max_blocks) stop 'Too many blocks'
        
        !Launch reduction kernel, this will init and synchronize itself
        call reduce_complex_kernel<<<blocks, threads>>>(psi_d, n, block_sums_d)
    
        !Copy result from device
        block_sums(1:blocks) = block_sums_d(1:blocks)
        
        !Final sum on host
        get_sum_complex = sum(block_sums(1:blocks))
        
    end function get_sum_complex
    
    
    ! rewritten to work with manual reduction kernel
    ! attributes(host) function get_sum_real(psi_d)
    
        ! implicit none
        
        ! integer(4) l, m
        ! real(fp_kind) :: get_sum_real
        ! real(fp_kind), device, dimension(nopiy,nopix) :: psi_d
        ! real(fp_kind),device :: get_sum_real_d
    
        ! get_sum_real_d = 0.0_fp_kind
    
        ! !$cuf kernel do (2) <<<*,*>>>
        ! do m = 1, nopix
          ! do l = 1, nopiy
              ! get_sum_real_d = get_sum_real_d + psi_d(l,m)
          ! enddo
        ! enddo
        
        ! get_sum_real = get_sum_real_d
    
    ! end function get_sum_real
    
    ! new sum for real, 2025-06-02 JB
    attributes(host) function get_sum_real(psi_d)
    
        implicit none
        
        integer, parameter :: threads = 256
        real(fp_kind), device, dimension(nopiy,nopix) :: psi_d
        real(fp_kind) :: get_sum_real
        integer :: n, blocks

        ! init numbers
        n = nopiy * nopix
        blocks = (n + threads - 1) / threads
        if (blocks > max_blocks) stop 'Too many blocks'

        ! Launch reduction kernel, this will init and synchronize itself
        call reduce_real_kernel<<<blocks, threads>>>(psi_d, n, block_sums_d)

        ! Copy result from device
        block_sums(1:blocks) = block_sums_d(1:blocks)
        
        ! Final sum on host
        get_sum_real = sum(block_sums(1:blocks))

    end function get_sum_real

    

end module
