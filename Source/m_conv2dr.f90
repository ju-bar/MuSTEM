!=======================================================================
! m_conv2dr.f90
! Reusable 2D real-to-real convolution module with hidden FFTs
! - GPU backend: cuFFT + CUDA Fortran kernels (device-only)
! - CPU backend: MKL (FFTW3 single-precision interface)
!
! Select backend at compile time:
!   -DGPU      -> build with cuFFT (NVHPC/PGI CUDA Fortran)
!   (default)  -> CPU with MKL/FFTW interface
!
! Public API:
!   call conv2dr_init(max_kernels) - initialize context with maximum number of kernels
!   call conv2dr_set_kernel(slot, kernel, ny, nx) - setup the kernel in the specified slot
!   call conv2dr_apply(slot, input, output, ny, nx) - apply the kernel to the input data
!   call conv2dr_unint() - cleanup context and deallocate resources
!
!=======================================================================
! Notes:
! - The data alignment is assumed to be data(ny, nx) for input, output and kernels.
! - Kernel coefficients are stored on a reduced grid (ny/2 + 1, nx)
! - The module stores kernel Fourier coefficients on host and on device (if GPU).
!   The kernel coefficients are not normalized with respect to array size.
!   Since the data coefficients will be normalized, the resulting convolution
!   preserves the original data intensity (not its L2 norm) if the kernel is normalized to 1.
!=======================================================================

module m_conv2dr

  use m_precision, only: fp_kind
#ifdef GPU
  use cudafor
  use cufft
  use cuda_array_library
#endif
  use CUFFT_wrapper

  implicit none

  !-------------------------------------------------------
  ! Derived types
  !-------------------------------------------------------
  type :: kernel_ft_t
    integer :: nx, ny                   ! kernel dimensions
#ifdef GPU
    complex(fp_kind), device, allocatable, dimension(:,:) :: d_ft ! device-side ft coefficients
    integer(c_int) :: cufft_plan_f      ! forward FFT plan
    integer(c_int) :: cufft_plan_b      ! backward FFT plan
#endif
    complex(fp_kind), allocatable, dimension(:,:) :: ft ! ft coefficients

  end type kernel_ft_t

  type :: conv_context_t
    integer :: max_kernels              ! maximum kernels allocated
    type(kernel_ft_t), allocatable :: kernels(:)
  end type conv_context_t

  !-------------------------------------------------------
  ! Module variable
  !-------------------------------------------------------
  type(conv_context_t) :: ctx
  
  !-------------------------------------------------------
  ! Module routines
  !-------------------------------------------------------
  public :: conv2dr_init, conv2dr_uninit
  public :: conv2dr_set_kernel
  private :: conv2dr_free_kernel
  
  !-------------------------------------------------------
  ! Module interfaces
  !-------------------------------------------------------
  interface conv2dr_apply
#ifdef GPU
    module procedure conv2dr_apply_d
#endif
    module procedure conv2dr_apply_h
  end interface conv2dr_apply
  
  contains

    
  !-------------------------------------------------------
  ! Initialize convolution context with maximum number of kernels
  !-------------------------------------------------------
  subroutine conv2dr_init(max_kernels)
    implicit none
    integer, intent(in) :: max_kernels
    integer :: i
    
    call conv2dr_uninit()  ! cleanup if already initialized

    ctx%max_kernels = max_kernels
    allocate(ctx%kernels(max_kernels))

    ! Initialize each kernel slot empty
    do i = 1, max_kernels
      ctx%kernels(i)%nx   = 0
      ctx%kernels(i)%ny   = 0
#ifdef GPU
      ctx%kernels(i)%cufft_plan_f = 0
      ctx%kernels(i)%cufft_plan_b = 0
#endif
    end do
  end subroutine conv2dr_init
  
  !-------------------------------------------------------
  ! Free the kernel in the specified slot
  !-------------------------------------------------------
  subroutine conv2dr_free_kernel(slot)
    implicit none
    integer, intent(in) :: slot
    
    if (ctx%max_kernels <= 0) return ! Nothing to do if not initialized
    if (.not. allocated(ctx%kernels)) return ! No kernels allocated
    if (slot < 1 .or. slot > ctx%max_kernels) return ! invalid slot

#ifdef GPU
    ! Destroy cuFFT plans
    if (ctx%kernels(slot)%cufft_plan_f /= 0) then
      call cufftDestroy(ctx%kernels(slot)%cufft_plan_f)
      ctx%kernels(slot)%cufft_plan_f = 0
    end if
    if (ctx%kernels(slot)%cufft_plan_b /= 0) then
      call cufftDestroy(ctx%kernels(slot)%cufft_plan_b)
      ctx%kernels(slot)%cufft_plan_b = 0
    end if
    ! Deallocate device-side FFT coefficients
    if (allocated(ctx%kernels(slot)%d_ft)) deallocate(ctx%kernels(slot)%d_ft)
#endif
    ! Deallocate host-side FFT coefficients
    if (allocated(ctx%kernels(slot)%ft)) deallocate(ctx%kernels(slot)%ft)

    ctx%kernels(slot)%nx   = 0
    ctx%kernels(slot)%ny   = 0

    return
  end subroutine conv2dr_free_kernel
  
  
  !-------------------------------------------------------
  ! This subroutine cleans up the context and deallocates resources.
  !-------------------------------------------------------
  subroutine conv2dr_uninit()
    implicit none
    integer :: i

    if (ctx%max_kernels <= 0) return ! Nothing to do if not initialized
    if (.NOT. allocated(ctx%kernels)) return ! Nothing to do if no kernels are available

    do i = 1, ctx%max_kernels
      call conv2dr_free_kernel(i) ! Free each kernel slot
    end do

    if (allocated(ctx%kernels)) deallocate(ctx%kernels)
    ctx%max_kernels = 0
    
    return
  end subroutine conv2dr_uninit
  
  
  !-------------------------------------------------------
  ! Setup the kernel in the specified slot
  !-------------------------------------------------------
  subroutine conv2dr_set_kernel(slot, kernel, ny, nx)
    implicit none
    integer, intent(in) :: slot ! slot index (1-based), less than or equal to ctx%max_kernels
    real(fp_kind), intent(in) :: kernel(ny,nx) ! kernel coefficients in real space
    integer, intent(in) :: nx, ny ! kernel dimensions
    real(fp_kind) :: norm_factor ! normalization factor for the kernel
    ! Locals
#ifdef GPU
    ! device arrays for FFT
    real(fp_kind), device, allocatable :: ker_d(:,:)
#endif
    
    if (ctx%max_kernels <= 0) return ! Nothing to do if not initialized
    if (.not. allocated(ctx%kernels)) return ! No kernels allocated
    if (slot < 1 .or. slot > ctx%max_kernels) return ! invalid slot
    
    call conv2dr_free_kernel(slot) ! Free existing kernel in the slot
    
#ifdef GPU
    ! Allocate storage
    allocate(ker_d(ny,nx)) ! Allocate device array for kernel
    ker_d = kernel ! Copy to device array
    ! Allocate device-side FFT coefficients
    allocate(ctx%kernels(slot)%d_ft(ny/2+1,nx))
    if (fp_kind==4) then
      call cufftPlan(ctx%kernels(slot)%cufft_plan_f, nx, ny, CUFFT_R2C)
      call cufftPlan(ctx%kernels(slot)%cufft_plan_b, nx, ny, CUFFT_C2R)
    else
      call cufftPlan(ctx%kernels(slot)%cufft_plan_f, nx, ny, CUFFT_D2Z)
      call cufftPlan(ctx%kernels(slot)%cufft_plan_b, nx, ny, CUFFT_Z2D)
    endif
    call cufftExec(ctx%kernels(slot)%cufft_plan_f, ker_d, ctx%kernels(slot)%d_ft)
    deallocate(ker_d) ! Free device array after execution
#endif
    norm_factor = sqrt(real(nx*ny, kind=fp_kind)) ! set result back to unnormalized Fourier coefficients
    allocate(ctx%kernels(slot)%ft(ny/2+1,nx))
    call fft2(ny, nx, kernel*norm_factor, ctx%kernels(slot)%ft)

    ctx%kernels(slot)%nx   = nx
    ctx%kernels(slot)%ny   = ny
    return
  end subroutine conv2dr_set_kernel
  
  
  
  !-------------------------------------------------------
  ! Apply convolution with a kernel on host
  !-------------------------------------------------------
  subroutine conv2dr_apply_h(slot, data_in, data_out, ny, nx)
    implicit none
    integer, intent(in) :: slot
    integer, intent(in) :: ny, nx
    real(fp_kind), intent(in)  :: data_in(ny,nx)
    real(fp_kind), intent(out) :: data_out(ny,nx)
    integer :: ny2
    complex(fp_kind), allocatable :: data_fft(:,:)

    ! Check context
    if (ctx%max_kernels <= 0) return
    if (.not. allocated(ctx%kernels)) return
    if (slot < 1 .or. slot > ctx%max_kernels) return
    if (ctx%kernels(slot)%nx /= nx .or. ctx%kernels(slot)%ny /= ny) return
    if (.NOT. allocated(ctx%kernels(slot)%ft)) return ! No kernel set
    
    ny2 = ny / 2 + 1 ! Half the size for FFT output

    allocate(data_fft(ny2,nx))
    call fft2(ny, nx, data_in, data_fft) ! applies also normalization 1/sqrt(nx*ny)
    data_fft = data_fft * ctx%kernels(slot)%ft
    call ifft2(ny, nx, data_fft, data_out) ! applies also normalization 1/sqrt(nx*ny)
    
    deallocate(data_fft)

    return
  end subroutine conv2dr_apply_h
  
  
  
#ifdef GPU
  
  !-------------------------------------------------------
  ! Apply convolution with a kernel on device
  !-------------------------------------------------------
  subroutine conv2dr_apply_d(slot, data_in, data_out, ny, nx)
    implicit none
    integer, intent(in) :: slot
    integer, intent(in) :: ny, nx
    real(fp_kind), device, intent(in)  :: data_in(ny,nx)
    real(fp_kind), device, intent(out) :: data_out(ny,nx)
    complex(fp_kind), device, allocatable :: data_fft(:,:)
    integer :: ny2
    real(fp_kind) :: norm_factor
    type(dim3) :: local_blocks

    ! Check context
    if (ctx%max_kernels <= 0) return
    if (.not. allocated(ctx%kernels)) return
    if (slot < 1 .or. slot > ctx%max_kernels) return
    if (ctx%kernels(slot)%nx /= nx .or. ctx%kernels(slot)%ny /= ny) return
    if (.NOT. allocated(ctx%kernels(slot)%d_ft)) return ! No kernel set
    
    ny2 = ny / 2 + 1 ! Half the size for FFT output
    
    local_blocks = dim3(ceiling(real(ny2)/threads%x), ceiling(real(nx)/threads%y), 1)
    norm_factor = 1.0_fp_kind / real(nx*ny, kind=fp_kind) ! FFT normalization factor
    allocate(data_fft(ny2,nx))
    call cufftExec(ctx%kernels(slot)%cufft_plan_f, data_in, data_fft)
    call cuda_multiplication<<<local_blocks,threads>>>(data_fft,ctx%kernels(slot)%d_ft, &
                            & data_fft,norm_factor,ny2,nx)
    call cufftExec(ctx%kernels(slot)%cufft_plan_b, data_fft, data_out)
    deallocate(data_fft)

    return
  end subroutine conv2dr_apply_d
  
#endif

end module m_conv2dr
