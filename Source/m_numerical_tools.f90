!--------------------------------------------------------------------------------
!
!  Copyright (C) 2017  L. J. Allen, H. G. Brown, A. J. D�Alfonso, S.D. Findlay, B. D. Forbes
!
!  modified:
!    2019-Dec-13, J. Barthel : code formatting, unirand, poirand1
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

!This module contains routines for statistical sampling of a Gaussian function (normal
!distribution), adapted from the library prob.f90 available from the website 
!<http://people.sc.fsu.edu/~jburkardt/f_src/prob/prob.html> and for cubic interpolation, 
!adapted from the library pppack.f90 abailable from the website
!<http://people.sc.fsu.edu/~jburkardt/f_src/pppack/pppack.html>. Only the required routines 
!were taken from these libraries and no modifications were made to those routines.
!The original MuSTEM codes used equivalent routines described in the book "Numerical Recipes 
!in FORTRAN: The Art of Scientific Computing" by Press et al. (1986), however the copyright 
!conditions of that work prevented those routines being included in this open source 
!implementation of MuSTEM.

module m_numerical_tools
	use m_precision

	interface ran1
		module procedure r8_uniform_01	
	end interface

contains        


	subroutine displace(tauin, tauout, urms, a0, idum)

		use m_precision

		implicit none

		real(fp_kind) tauin(3), tauout(3), a0(3)
		real(fp_kind) urms
		integer(4) idum
		real(fp_kind) xran, yran, zran
		!real(fp_kind) gasdev

		xran = gasdev(idum)
		yran = gasdev(idum)
		!zran = gasdev(idum)

		xran = xran * urms / a0(1)
		yran = yran * urms / a0(2)
		!zran = zran * urms / a0(3)

		tauout(1) = tauin(1) + xran
		tauout(2) = tauin(2) + yran
		tauout(3) = tauin(3)
		!	tauout(3) = tauin(3) + zran

	end subroutine displace
	
	
	function unirand(seed)
		use m_precision, only: fp_kind
		implicit none
		real(fp_kind) :: unirand
		integer(kind=4) :: seed
		real(kind=8) :: x
		x = r8_uniform_01(seed)
		unirand = real(x, kind=fp_kind)
	end function
	
	
	function poirand1(m, seed)
!*****************************************************************************80
!
!   POIRAND1 samples of the Poisson PDF with given mean value
!
!   Discussion:
!		Donold Knuth's algorithm is fast for small mean values < 30
!
!   Input:
!       real(fp_kind) m = PDF mean value
!       integer(4) seed = rng seed
!
!   Output:
!       integer(4) Poisson random sample
!
		use m_precision, only: fp_kind
		implicit none
		integer(kind=4) :: poirand1
		real(fp_kind) :: m
		integer(kind=4) :: seed
		integer(kind=4) :: k
		real(kind=8) :: md, p, l
		md = dabs(dble(m))
		l = dexp(-md)
		k = 0
		p = 1.D+0
		do while (p > l)
			k = k + 1
			p = p * r8_uniform_01(seed)
		end do
		poirand1 = k-1
	end function
	
	
		
	subroutine normal_sample ( a, b, seed, x )

!*****************************************************************************80
!
!! NORMAL_SAMPLE samples the Normal PDF.
!
!  Discussion:
!
!    The Box-Muller method is used.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0 < B.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number
!    generator.
!
!    Output, real ( kind = 8 ) X, a sample of the PDF.
!
		implicit none

		real ( kind = 8 ) a
		real ( kind = 8 ) b
		integer ( kind = 4 ) seed
		real ( kind = 8 ) x

		call normal_01_sample ( seed, x )

		x = a + b * x

		return
	end subroutine normal_sample

	!This function is a wrapper to the routine normal_01_sample that mimics
	!the interface to the gasdev function from Numerical Recipes (see comment
	!at the beginning of this file).
	function gasdev(seed)
		real(fp_kind)::gasdev
		real*8::x
		integer*4,intent(inout)::seed

		call normal_01_sample(seed,x)
		gasdev = real(x, fp_kind)
	end function

	subroutine normal_01_sample ( seed, x )

!*****************************************************************************80
!
!! NORMAL_01_SAMPLE samples the standard normal probability distribution.
!
!  Discussion:
!
!    The standard normal probability distribution function (PDF) has
!    mean 0 and standard deviation 1.
!
!    The Box-Muller method is used, which is efficient, but
!    generates two values at a time.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 August 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number
!    generator.
!
!    Output, real ( kind = 8 ) X, a sample of the standard normal PDF.
!
		implicit none

		real ( kind = 8 ) r1
		real ( kind = 8 ) r2
		real ( kind = 8 ), parameter :: r8_pi = 3.141592653589793D+00
		!real ( kind = 8 ) r8_uniform_01
		integer ( kind = 4 ) seed
		real ( kind = 8 ) x

		r1 = r8_uniform_01 ( seed )
		r2 = r8_uniform_01 ( seed )
		x = sqrt ( -2.0D+00 * log ( r1 ) ) * cos ( 2.0D+00 * r8_pi * r2 )

		return
	end subroutine normal_01_sample
	
	

	function r8_uniform_01 ( seed )

!*****************************************************************************80
!
!! R8_UNIFORM_01 returns a unit pseudorandom R8.
!
!  Discussion:
!
!    An R8 is a real ( kind = 8 ) value.
!
!    For now, the input quantity SEED is an integer ( kind = 4 ) variable.
!
!    This routine implements the recursion
!
!      seed = 16807 * seed mod ( 2^31 - 1 )
!      r8_uniform_01 = seed / ( 2^31 - 1 )
!
!    The integer arithmetic never requires more than 32 bits,
!    including a sign bit.
!
!    If the initial seed is 12345, then the first three computations are
!
!      Input     Output      R8_UNIFORM_01
!      SEED      SEED
!
!         12345   207482415  0.096616
!     207482415  1790989824  0.833995
!    1790989824  2035175616  0.947702
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley Interscience, page 95, 1998.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which should
!    NOT be 0. On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R8_UNIFORM_01, a new pseudorandom variate,
!    strictly between 0 and 1.
!
		implicit none

		!integer ( kind = 4 ) i4_huge
		integer ( kind = 4 ) k
		real ( kind = 8 ) r8_uniform_01
		integer ( kind = 4 ) seed

		if ( seed == 0 ) then
			write ( *, '(a)' ) ' '
			write ( *, '(a)' ) 'R8_UNIFORM_01 - Fatal error!'
			write ( *, '(a)' ) '  Input value of SEED = 0.'
			stop 1
		end if

		k = seed / 127773

		seed = 16807 * ( seed - k * 127773 ) - k * 2836

		if ( seed < 0 ) then
		seed = seed + i4_huge ( )
		end if
		!
		!  Although SEED can be represented exactly as a 32 bit integer,
		!  it generally cannot be represented exactly as a 32 bit real number!
		!
		r8_uniform_01 = real ( seed, kind = 8 ) * 4.656612875D-10

		return
  end function
  
  
!*****************************************************************************80
!
!! R8_PDF_REJ returns a pseudorandom R8 following a probability distribution
!             function (PDF) given by samples x, y in arrays of length n.
!
!  Discussion:
!
!    An R8 is a real ( kind = 8 ) value.
!
!    For now, the input quantity SEED is an integer ( kind = 4 ) variable.
!
!    This routine implements the rejection sampling method.
!
!    It makes use of function R8_UNIFORM_01 to generate uniform random numbers
!    for the range of samples x and y.
!
!    A random is x is accepted as return value if the random y is smaller than
!    or equal to the value of y at x.
!
!    Linear interpolation is used between samples of x.
!
!    x is expected to be in ascending order with corresponding y = PDF(x).
!
!    It is not required to input a normalized PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 August 2023
!
!  Author:
!
!    Juri Barthel
!
!  Reference:
!
!    rejection sampling:  https://en.wikipedia.org/wiki/Rejection_sampling
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which should
!    NOT be 0. On output, SEED has been updated.
!
!    Input, integer ( kind = 4 ) N, number of samples in following input.
!    N must be greater than 0. N = 1 assumes a delta ditribution.
!
!    Input, real ( kind = 8 , dimension = (2,N) ) X,Y samples of Y=PDF(X) corresponding
!    to the input X.
!
!    Output, real ( kind = 8 ) R8_PDF_REJ, a new pseudorandom variate,
!    strictly in the range of minval(pdf(1,:)) and maxval(pdf(1,:)).
!
!
  function r8_distr_rej ( seed, n, pdf )
  
		implicit none
    
    integer ( kind = 4 ), parameter :: num_rej_max = 1048576

		integer ( kind = 4 ), intent(in) :: n
    real ( kind = 8 ), intent(in), dimension(2,n) :: pdf
		integer ( kind = 4 ), intent(inout) :: seed
    real ( kind = 8 ) :: r8_distr_rej
    real ( kind = 8 ) :: x0, x1, y0, y1 ! x and y ranges
    real ( kind = 8 ) :: urn_x, urn_y ! uniform random numbers
    real ( kind = 8 ) :: pdf_urn_x ! interpolated PDF value at urn_x
    integer ( kind = 4 ) :: ip ! interpolation base
    real ( kind = 8 ) :: frc ! interpolation step
    integer ( kind = 4 ) :: irej ! rejection count
    
    if ( n == 0 ) then
			write ( *, '(a)' ) ' '
			write ( *, '(a)' ) 'R8_PDF_REJ - Fatal error!'
			write ( *, '(a)' ) '  Input value of N = 0.'
			stop 1
    end if
    
    if ( n == 1 ) then ! delta distribution
      
      r8_distr_rej = pdf(1,1)
      
    end if
    
    if ( n > 1 ) then
      
      ! determine range of x samples
      x0 = minval(pdf(1,:))
      x1 = maxval(pdf(1,:))
      ! determine range of y samples
      y0 = minval(pdf(2,:))
      y1 = maxval(pdf(2,:))
      !
      if (y0 < 0.0D+0 .or. y1 < 0.0D+0 ) then
        write ( *, '(a)' ) ' '
			  write ( *, '(a)' ) 'R8_PDF_REJ - Fatal error!'
			  write ( *, '(a)' ) '  Negative probability input encountered.'
			  stop 2
      end if
      
      if (x1 > x0) then ! good range in x
        
        if (y1 > y0) then ! good range in y
          
          ! implement rejection sampling
          do irej=1, num_rej_max
            
            ! get a random x, uniform sample of the x range
            urn_x = x0 + (x1 - x0) * r8_uniform_01( SEED )
            ! get a random y, uniform sample of the y range
            urn_y = y0 * r8_uniform_01( SEED )
            
            ! get the value of the PDF a urn_x by linear interpolation
            ! on the samples (x, y)
            pdf_urn_x = 0.0D+0 ! init some value of the PDF, use zero to catch problems
            do ip=1, n
              
              ! find first x sample that is bigger then urn_x
              if ( pdf(1,ip) > urn_x ) then
                
                ! handle strange cases first
                if ( ip == 1 ) then
                  
                  pdf_urn_x = pdf(2,1)
                  
                  exit ! exit loop over ip
                  
                else if ( ip > 1 .and. ip <= n ) then ! in interpolation range
                  
                  ! linear interpolation
                  pdf_urn_x = pdf(2,ip-1) + ( pdf(2,ip) - pdf(2,ip-1) )  &
                     / ( pdf(1,ip) - pdf(1,ip-1) ) * ( urn_x - pdf(1,ip-1) )
                  
                  exit ! exit loop over ip
                  
                end if
                
              end if
              
            end do
            
            ! loopexit, check whether this was because none of the if conditions aplied
            if (ip > n) then
              
               pdf_urn_x = pdf(2,n) ! use end value
               
            end if
            
            ! acceptance ?
            if ( urn_y <= pdf_urn_x ) then ! yes, accept this random x
              
              r8_distr_rej = urn_x
              
              return
              
            end if
            
          end do
          
          ! if the code gets to this point the number of rejections
          ! is too large. Something is fishy here. Report a problem.
          write ( *, '(a)' ) ' '
			    write ( *, '(a)' ) 'R8_PDF_REJ - Fatal error!'
			    write ( *, '(a)' ) '  Too many rejections.'
			    stop 3
          
        else ! invalid range in y, assuming uniform distribution
          
          r8_distr_rej = x0 + (x1 - x0) * r8_uniform_01( SEED )
          
          return
          
        end if
        
      else ! invalid range in x, assuming delta distribution
        
        r8_distr_rej = pdf(1,1)
        
        return
        
      end if
    end if
  
  end function
	
	

	function i4_huge ( )

!*****************************************************************************80
!
!! I4_HUGE returns a "huge" I4.
!
!  Discussion:
!
!    On an IEEE 32 bit machine, I4_HUGE should be 2^31 - 1, and its
!    bit pattern should be
!
!     01111111111111111111111111111111
!
!    In this case, its numerical value is 2147483647.
!
!    Using the Dec/Compaq/HP Alpha FORTRAN compiler FORT, I could
!    use I4_HUGE() and HUGE interchangeably.
!
!    However, when using the G95, the values returned by HUGE were
!    not equal to 2147483647, apparently, and were causing severe
!    and obscure errors in my random number generator, which needs to
!    add I4_HUGE to the seed whenever the seed is negative.  So I
!    am backing away from invoking HUGE, whereas I4_HUGE is under
!    my control.
!
!    Explanation: because under G95 the default integer type is 64 bits!
!    So HUGE ( 1 ) = a very very huge integer indeed, whereas
!    I4_HUGE ( ) = the same old 32 bit big value.
!
!    An I4 is an integer ( kind = 4 ) value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 January 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) I4_HUGE, a "huge" I4.
!
		implicit none

		integer ( kind = 4 ) i4
		integer ( kind = 4 ) i4_huge

		i4_huge = 2147483647

		return
	end function



	subroutine cubspl ( tau, c, n, ibcbeg, ibcend )

!*****************************************************************************80
!
!! CUBSPL defines an interpolatory cubic spline.
!
!  Discussion:
!
!    A tridiagonal linear system for the unknown slopes S(I) of
!    F at TAU(I), I=1,..., N, is generated and then solved by Gauss
!    elimination, with S(I) ending up in C(2,I), for all I.
!
!  Modified:
!
!    14 February 2007
!
!  Author:
!
!    Original FORTRAN77 version by Carl de Boor.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Carl de Boor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) TAU(N), the abscissas or X values of
!    the data points.  The entries of TAU are assumed to be
!    strictly increasing.
!
!    Input, integer ( kind = 4 ) N, the number of data points.  N is
!    assumed to be at least 2.
!
!    Input/output, real ( kind = 8 ) C(4,N).
!    On input, if IBCBEG or IBCBEG is 1 or 2, then C(2,1)
!    or C(2,N) should have been set to the desired derivative
!    values, as described further under IBCBEG and IBCEND.
!    On output, C contains the polynomial coefficients of
!    the cubic interpolating spline with interior knots
!    TAU(2) through TAU(N-1).
!    In the interval interval (TAU(I), TAU(I+1)), the spline
!    F is given by
!      F(X) = 
!        C(1,I) + 
!        C(2,I) * H +
!        C(3,I) * H^2 / 2 + 
!        C(4,I) * H^3 / 6.
!    where H=X-TAU(I).  The routine PPVALU may be used to
!    evaluate F or its derivatives from TAU, C, L=N-1,
!    and K=4.
!
!    Input, integer ( kind = 4 ) IBCBEG, IBCEND, boundary condition indicators.
!    IBCBEG = 0 means no boundary condition at TAU(1) is given.
!    In this case, the "not-a-knot condition" is used.  That
!    is, the jump in the third derivative across TAU(2) is
!    forced to zero.  Thus the first and the second cubic
!    polynomial pieces are made to coincide.
!    IBCBEG = 1 means the slope at TAU(1) is to equal the
!    input value C(2,1).
!    IBCBEG = 2 means the second derivative at TAU(1) is
!    to equal C(2,1).
!    IBCEND = 0, 1, or 2 has analogous meaning concerning the
!    boundary condition at TAU(N), with the additional
!    information taken from C(2,N).
!
		implicit none

		integer ( kind = 4 ) n

		real ( kind = fp_kind ) c(4,n)
		real ( kind = fp_kind ) divdf1
		real ( kind = fp_kind ) divdf3
		real ( kind = fp_kind ) dtau
		real ( kind = fp_kind ) g
		integer ( kind = 4 ) i
		integer ( kind = 4 ) ibcbeg
		integer ( kind = 4 ) ibcend
		real ( kind = fp_kind ) tau(n)
!
!  C(3,*) and C(4,*) are used initially for temporary storage.
!
!  Store first differences of the TAU sequence in C(3,*).
!
!  Store first divided difference of data in C(4,*).
!
		do i = 2, n
			c(3,i) = tau(i) - tau(i-1)
		end do

		do i = 2, n 
			c(4,i) = ( c(1,i) - c(1,i-1) ) / ( tau(i) - tau(i-1) )
		end do
!
!  Construct the first equation from the boundary condition
!  at the left endpoint, of the form:
!
!    C(4,1) * S(1) + C(3,1) * S(2) = C(2,1)
!
!  IBCBEG = 0: Not-a-knot
!
		if ( ibcbeg == 0 ) then

			if ( n <= 2 ) then
				c(4,1) = 1.0D+00
				c(3,1) = 1.0D+00
				c(2,1) = 2.0D+00 * c(4,2)
				go to 120
			end if

			c(4,1) = c(3,3)
			c(3,1) = c(3,2) + c(3,3)
			c(2,1) = ( ( c(3,2) + 2.0D+00 * c(3,1) ) * c(4,2) * c(3,3) &
			  + c(3,2)**2 * c(4,3) ) / c(3,1)
!
!  IBCBEG = 1: derivative specified.
!
		else if ( ibcbeg == 1 ) then

			c(4,1) = 1.0D+00
			c(3,1) = 0.0D+00

			if ( n == 2 ) then
				go to 120
			end if
!
!  Second derivative prescribed at left end.
!
		else

			c(4,1) = 2.0D+00
			c(3,1) = 1.0D+00
			c(2,1) = 3.0D+00 * c(4,2) - c(3,2) / 2.0D+00 * c(2,1)

			if ( n == 2 ) then
				go to 120
			end if

		end if
!
!  If there are interior knots, generate the corresponding
!  equations and carry out the forward pass of Gauss elimination,
!  after which the I-th equation reads:
!
!    C(4,I) * S(I) + C(3,I) * S(I+1) = C(2,I).
!
		do i = 2, n-1
			g = -c(3,i+1) / c(4,i-1)
			c(2,i) = g * c(2,i-1) + 3.0D+00 * ( c(3,i) * c(4,i+1) + c(3,i+1) * c(4,i) )
			c(4,i) = g * c(3,i-1) + 2.0D+00 * ( c(3,i) + c(3,i+1))
		end do
!
!  Construct the last equation from the second boundary condition, of
!  the form
!
!    -G * C(4,N-1) * S(N-1) + C(4,N) * S(N) = C(2,N)
!
!  If slope is prescribed at right end, one can go directly to
!  back-substitution, since the C array happens to be set up just
!  right for it at this point.
!
		if ( ibcend == 1 ) then
			go to 160
		end if

		if ( 1 < ibcend ) then
			go to 110
		end if
 
90		continue
!
!  Not-a-knot and 3 <= N, and either 3 < N or also not-a-knot
!  at left end point.
!
		if ( n /= 3 .or. ibcbeg /= 0 ) then
			g = c(3,n-1) + c(3,n)
			c(2,n) = ( ( c(3,n) + 2.0D+00 * g ) * c(4,n) * c(3,n-1) + c(3,n)**2 &
			  * ( c(1,n-1) - c(1,n-2) ) / c(3,n-1) ) / g
			g = - g / c(4,n-1)
			c(4,n) = c(3,n-1)
			c(4,n) = c(4,n) + g * c(3,n-1)
			c(2,n) = ( g * c(2,n-1) + c(2,n) ) / c(4,n)
			go to 160
		end if
!
!  N = 3 and not-a-knot also at left.
!
100		continue
 
		c(2,n) = 2.0D+00 * c(4,n)
		c(4,n) = 1.0D+00
		g = -1.0D+00 / c(4,n-1)
		c(4,n) = c(4,n) - c(3,n-1) / c(4,n-1)
		c(2,n) = ( g * c(2,n-1) + c(2,n) ) / c(4,n)
		go to 160
!
!  IBCEND = 2: Second derivative prescribed at right endpoint.
!
110		continue
 
		c(2,n) = 3.0D+00 * c(4,n) + c(3,n) / 2.0D+00 * c(2,n)
		c(4,n) = 2.0D+00
		g = -1.0D+00 / c(4,n-1)
		c(4,n) = c(4,n) - c(3,n-1) / c(4,n-1)
		c(2,n) = ( g * c(2,n-1) + c(2,n) ) / c(4,n)
		go to 160
!
!  N = 2.
!
120		continue
  
		if ( ibcend == 2  ) then

			c(2,n) = 3.0D+00 * c(4,n) + c(3,n) / 2.0D+00 * c(2,n)
			c(4,n) = 2.0D+00
			g = -1.0D+00 / c(4,n-1)
			c(4,n) = c(4,n) - c(3,n-1) / c(4,n-1)
			c(2,n) = ( g * c(2,n-1) + c(2,n) ) / c(4,n)

		else if ( ibcend == 0 .and. ibcbeg /= 0 ) then

			c(2,n) = 2.0D+00 * c(4,n)
			c(4,n) = 1.0D+00
			g = -1.0D+00 / c(4,n-1)
			c(4,n) = c(4,n) - c(3,n-1) / c(4,n-1)
			c(2,n) = ( g * c(2,n-1) + c(2,n) ) / c(4,n)

		else if ( ibcend == 0 .and. ibcbeg == 0 ) then

			c(2,n) = c(4,n)

		end if
!
!  Back solve the upper triangular system 
!
!    C(4,I) * S(I) + C(3,I) * S(I+1) = B(I)
!
!  for the slopes C(2,I), given that S(N) is already known.
!
160		continue
 
		do i = n-1, 1, -1
			c(2,i) = ( c(2,i) - c(3,i) * c(2,i+1) ) / c(4,i)
		end do
!
!  Generate cubic coefficients in each interval, that is, the
!  derivatives at its left endpoint, from value and slope at its
!  endpoints.
!
		do i = 2, n
			dtau = c(3,i)
			divdf1 = ( c(1,i) - c(1,i-1) ) / dtau
			divdf3 = c(2,i-1) + c(2,i) - 2.0D+00 * divdf1
			c(3,i-1) = 2.0D+00 * ( divdf1 - c(2,i-1) - divdf3 ) / dtau
			c(4,i-1) = 6.0D+00 * divdf3 / dtau**2
		end do
 
		return
	end subroutine



	function ppvalu ( break, coef, l, k, x, jderiv )

!*****************************************************************************80
!
!! PPVALU evaluates a piecewise polynomial function or its derivative.
!
!  Discussion:
!
!    PPVALU calculates the value at X of the JDERIV-th derivative of
!    the piecewise polynomial function F from its piecewise
!    polynomial representation.
!
!    The interval index I, appropriate for X, is found through a
!    call to INTERV.  The formula for the JDERIV-th derivative
!    of F is then evaluated by nested multiplication.
!
!    The J-th derivative of F is given by:
!      (d^J) F(X) = 
!        COEF(J+1,I) + H * (
!        COEF(J+2,I) + H * (
!        ...
!        COEF(K-1,I) + H * (
!        COEF(K,  I) / (K-J-1) ) / (K-J-2) ... ) / 2 ) / 1
!    with
!      H = X - BREAK(I)
!    and
!      I = max ( 1, max ( J, BREAK(J) <= X, 1 <= J <= L ) ).
!
!  Modified:
!
!    16 February 2007
!
!  Author:
!
!    Original FORTRAN77 version by Carl de Boor.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Carl de Boor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) BREAK(L+1), real COEF(*), integer L, the
!    piecewise polynomial representation of the function F to be evaluated.
!
!    Input, integer ( kind = 4 ) K, the order of the polynomial pieces that 
!    make up the function F.  The usual value for K is 4, signifying a 
!    piecewise cubic polynomial.
!
!    Input, real ( kind = 8 ) X, the point at which to evaluate F or
!    of its derivatives.
!
!    Input, integer ( kind = 4 ) JDERIV, the order of the derivative to be
!    evaluated.  If JDERIV is 0, then F itself is evaluated,
!    which is actually the most common case.  It is assumed
!    that JDERIV is zero or positive.
!
!    Output, real ( kind = 8 ) PPVALU, the value of the JDERIV-th
!    derivative of F at X.
!
		implicit none

		integer ( kind = 4 ) k
		integer ( kind = 4 ) l

		real ( kind = fp_kind ) break(l+1)
		real ( kind = fp_kind ) coef(k,l)
		real ( kind = fp_kind ) fmmjdr
		real ( kind = fp_kind ) h
		integer ( kind = 4 ) i
		integer ( kind = 4 ) jderiv
		integer ( kind = 4 ) m
		integer ( kind = 4 ) ndummy
		real ( kind = fp_kind ) ppvalu
		real ( kind = fp_kind ) value
		real ( kind = fp_kind ) x

		value = 0.0D+00

		fmmjdr = k - jderiv
!
!  Derivatives of order K or higher are identically zero.
!
		if ( k <= jderiv ) then
			return
		end if
!
!  Find the index I of the largest breakpoint to the left of X.
!
		call interv ( break, l+1, x, i, ndummy )
!
!  Evaluate the JDERIV-th derivative of the I-th polynomial piece at X.
!
		h = x - break(i)
		m = k
 
		do

			value = ( value / fmmjdr ) * h + coef(m,i)
			m = m - 1
			fmmjdr = fmmjdr - 1.0D+00

			if ( fmmjdr <= 0.0D+00 ) then
				exit
			end if

		end do

		ppvalu = value

		return
	end function



	subroutine interv ( xt, lxt, x, left, mflag )

!*****************************************************************************80
!
!! INTERV brackets a real value in an ascending vector of values.
!
!  Discussion:
!
!    The XT array is a set of increasing values.  The goal of the routine
!    is to determine the largest index I so that 
!
!      XT(I) < XT(LXT)  and  XT(I) <= X.
!
!    The routine is designed to be efficient in the common situation
!    that it is called repeatedly, with X taken from an increasing
!    or decreasing sequence.
!
!    This will happen when a piecewise polynomial is to be graphed.
!    The first guess for LEFT is therefore taken to be the value
!    returned at the previous call and stored in the local variable ILO.
!
!    A first check ascertains that ILO < LXT.  This is necessary
!    since the present call may have nothing to do with the previous
!    call.  Then, if 
!      XT(ILO) <= X < XT(ILO+1), 
!    we set LEFT = ILO and are done after just three comparisons.
!
!    Otherwise, we repeatedly double the difference ISTEP = IHI - ILO
!    while also moving ILO and IHI in the direction of X, until
!      XT(ILO) <= X < XT(IHI)
!    after which we use bisection to get, in addition, ILO + 1 = IHI.
!    The value LEFT = ILO is then returned.
!
!    Thanks to Daniel Gloger for pointing out an important modification
!    to the routine, so that the piecewise polynomial in B-form is
!    left-continuous at the right endpoint of the basic interval,
!    17 April 2014.
!
!  Modified:
!
!    17 April 2014
!
!  Author:
!
!    Original FORTRAN77 version by Carl de Boor.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Carl de Boor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XT(LXT), a nondecreasing sequence of values.
!
!    Input, integer ( kind = 4 ) LXT, the dimension of XT.
!
!    Input, real ( kind = 8 ) X, the point whose location with 
!    respect to the sequence XT is to be determined.
!
!    Output, integer ( kind = 4 ) LEFT, the index of the bracketing value:
!      1     if             X  <  XT(1)
!      I     if   XT(I)  <= X  < XT(I+1)
!      I     if   XT(I)  <  X == XT(I+1) == XT(LXT)
!
!    Output, integer ( kind = 4 ) MFLAG, indicates whether X lies within the
!    range of the data.
!    -1:            X  <  XT(1)
!     0: XT(I)   <= X  < XT(I+1)
!    +1: XT(LXT) <  X
!
		implicit none

		integer ( kind = 4 ) lxt

		integer ( kind = 4 ) left
		integer ( kind = 4 ) mflag
		integer ( kind = 4 ) ihi
		integer ( kind = 4 ), save :: ilo = 1
		integer ( kind = 4 ) istep
		integer ( kind = 4 ) middle
		real ( kind = fp_kind ) x
		real ( kind = fp_kind ) xt(lxt)

		ihi = ilo + 1

		if ( lxt <= ihi ) then

			if ( xt(lxt) <= x ) then
				go to 110
			end if

			if ( lxt <= 1 ) then
				mflag = -1
				left = 1
				return
			end if

			ilo = lxt - 1
			ihi = lxt

		end if

		if ( xt(ihi) <= x ) then
			go to 20
		end if

		if ( xt(ilo) <= x ) then
			mflag = 0
			left = ilo
			return
		end if
!
!  Now X < XT(ILO).  Decrease ILO to capture X.
!
		istep = 1

10		continue

		ihi = ilo
		ilo = ihi - istep

		if ( 1 < ilo ) then
			if ( xt(ilo) <= x ) then
				go to 50
			end if
			istep = istep * 2
			go to 10
		end if

		ilo = 1

		if ( x < xt(1) ) then
			mflag = -1
			left = 1
			return
		end if

		go to 50
!
!  Now XT(IHI) <= X.  Increase IHI to capture X.
!
20		continue

		istep = 1

30		continue

		ilo = ihi
		ihi = ilo + istep

		if ( ihi < lxt ) then

			if ( x < xt(ihi) ) then
				go to 50
			end if

			istep = istep * 2
			go to 30

		end if

		if ( xt(lxt) <= x ) then
			go to 110
		end if
!
!  Now XT(ILO) < = X < XT(IHI).  Narrow the interval.
!
		ihi = lxt

50		continue

		do

			middle = ( ilo + ihi ) / 2

			if ( middle == ilo ) then
				mflag = 0
				left = ilo
				return
			end if
			!
			!  It is assumed that MIDDLE = ILO in case IHI = ILO+1.
			!
			if ( xt(middle) <= x ) then
				ilo = middle
			else
				ihi = middle
			end if

		end do
!
!  Set output and return.
!
110		continue

		mflag = 1

		if ( x == xt(lxt) ) then
			mflag = 0
		end if

		do left = lxt - 1, 1, -1
			if ( xt(left) < xt(lxt) ) then
				return
			end if
		end do

		return
  end subroutine
  
  !
  ! returns the integral of the curve described by points
  ! calculated by the trapezoidal rule
  ! https://en.wikipedia.org/wiki/Trapezoidal_rule
  !
  ! uses a double precision accumulator but input and
  ! output is with precision "fp_kind"
  !
  ! INPUT
  !   integer(4) n = number of points
  !   real(fp_kind) points(2,n) = list of n tuples (x_i,y_i)
  !
  ! RETURNS
  !   real(fp_kind) = the numerical integral for points
  !
  function trapz(n, points)
  
    implicit none
    
    integer(4), intent(in) :: n ! number of data points
    real(fp_kind), intent(in) :: points(2, n) ! data points (x_i, y_i)
    
    real(fp_kind) :: trapz ! the integral
    
    integer(4) :: i
    real(fp_kind) :: ym, dx
    real(8) :: s
    
    s = 0.0D+0
    
    if (n > 1) then
      do i=2, n
        ym = 0.5_fp_kind * (points(2,i) + points(2,i-1))
        dx = points(1,i) - points(1,i-1)
        s = s + real(ym * dx, kind=8)
      end do
    end if
    
    trapz = real(s, kind=fp_kind)
  
  end function
  

end module