
!     ..................................................................
!
!  The rest of this file is taken from 
!  fftpack from http://www.psc.edu/~burkardt/src/fftpack/fftpack.html
!  and adapted for double precision. It is not formatted as the rest of
!  the CP-PAW source code
!                                                J. Schimpl 13.12.2001
!     ******************************************************************
subroutine cfftb ( n, c, wsave )
!
!*******************************************************************************
!
!! CFFTB computes the backward complex discrete Fourier transform.
!
!
!  Discussion:
!
!    This process is sometimes called Fourier synthesis.
!
!    CFFTB computes a complex periodic sequence from its Fourier coefficients.
!
!    A call of CFFTF followed by a call of CFFTB will multiply the
!    sequence by N.  In other words, the transforms are not normalized.
!
!    The array WSAVE must be initialized by CFFTI.
!
!    The transform is defined by:
!
!      C_out(J) = sum ( 1 <= K <= N ) 
!        C_in(K) * exp ( sqrt ( - 1 ) * ( J - 1 ) * ( K - 1 ) * 2 * PI / N )
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!    P N Swarztrauber, 
!    Vectorizing the FFT's, 
!    in Parallel Computations (G. Rodrigue, editor), 
!    Academic Press, 1982, pages 51-83.
!
!    B L Buzbee, 
!    The SLATEC Common Math Library, 
!    in Sources and Development of Mathematical Software (W. Cowell, editor),
!    Prentice Hall, 1984, pages 302-318.
!
!  Parameters:
!
!    Input, integer N, the length of the sequence to be transformed.  
!    The method is more efficient when N is the product of small primes.
!
!    Input/output, complex C(N).
!    On input, C contains the sequence of Fourier coefficients.
!    On output, C contains the sequence of data values that correspond
!    to the input coefficients.
!
!    Input, real WSAVE(4*N+15).  The array must be initialized by calling 
!    CFFTI.  A different WSAVE array must be used for each different
!    value of N.  
!
  integer(4) n
!
  complex(8) c(n)
  real(8) wsave(4*n+15)
!
  if ( n <= 1 ) then
    return
  end if

  call cfftb1 ( n, c, wsave(1), wsave(2*n+1), wsave(4*n+1) )

  return
end
subroutine cfftb1 ( n, c, ch, wa, ifac )
!
!*******************************************************************************
!
!! CFFTB1 is a lower-level routine used by CFFTB.
!
!
!  Modified:
!
!    12 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
!    Input, integer N, the length of the sequence to be transformed.  
!
!    Input/output, complex C(N).
!    On input, C contains the sequence of Fourier coefficients.
!    On output, C contains the sequence of data values that correspond
!    to the input coefficients.
!
!    Input, complex CH(N).
!
!    Input, real WA(2*N).
!
!    Input, integer IFAC(15).
!    IFAC(1) = N, the number that was factored.
!    IFAC(2) = NF, the number of factors.
!    IFAC(3:2+NF), the factors.
!
  integer(4) n
!
  complex(8) c(n)
  complex(8) ch(n)
  integer(4) idl1
  integer(4) ido
  integer(4) ifac(15)
  integer(4) ip
  integer(4) iw
  integer(4) ix2
  integer(4) ix3
  integer(4) ix4
  integer(4) k1
  integer(4) l1
  integer(4) l2
  integer(4) na
  integer(4) nac
  integer(4) nf
  real(8)wa(2*n)
!
  nf = ifac(2)
  na = 0
  l1 = 1
  iw = 1

  do k1 = 1, nf

    ip = ifac(k1+2)
    l2 = ip * l1
    ido = n / l2
    idl1 = 2 * ido * l1

    if ( ip == 4 ) then

      ix2 = iw + 2 * ido
      ix3 = ix2 + 2 * ido

      if ( na == 0 ) then
        call passb4 ( 2*ido, l1, c, ch, wa(iw), wa(ix2), wa(ix3) )
      else
        call passb4 ( 2*ido, l1, ch, c, wa(iw), wa(ix2), wa(ix3) )
      end if

      na = 1 - na

    else if ( ip == 2 ) then

      if ( na == 0 ) then
        call passb2 ( 2*ido, l1, c, ch, wa(iw) )
      else
        call passb2 ( 2*ido, l1, ch, c, wa(iw) )
      end if

      na = 1 - na

    else if ( ip == 3 ) then

      ix2 = iw + 2 * ido
 
      if ( na == 0 ) then
        call passb3 ( 2*ido, l1, c, ch, wa(iw), wa(ix2) )
      else
        call passb3 ( 2*ido, l1, ch, c, wa(iw), wa(ix2) )
      end if

      na = 1 - na

    else if ( ip == 5 ) then

      ix2 = iw + 2 * ido
      ix3 = ix2 + 2 * ido
      ix4 = ix3 + 2 * ido

      if ( na == 0 ) then
        call passb5 ( 2*ido, l1, c, ch, wa(iw), wa(ix2), wa(ix3), wa(ix4) )
      else
        call passb5 ( 2*ido, l1, ch, c, wa(iw), wa(ix2), wa(ix3), wa(ix4) )
      end if

      na = 1 - na

    else

      if ( na == 0 ) then
        call passb ( nac, 2*ido, ip, l1, idl1, c, c, c, ch, ch, wa(iw) )
      else
        call passb ( nac, 2*ido, ip, l1, idl1, ch, ch, ch, c, c, wa(iw) )
      end if

      if ( nac /= 0 ) then
        na = 1 - na
      end if

    end if

    l1 = l2
    iw = iw + ( ip - 1 ) * 2 * ido

  end do

  if ( na /= 0 ) then
    c(1:n) = ch(1:n)
  end if

  return
end
subroutine cfftb_2d ( ldf, n, f, wsave )
!
!*******************************************************************************
!
!! CFFTB_2D computes a backward two dimensional complex fast Fourier transform.
!
!
!  Discussion:
!
!    The routine computes the backward two dimensional fast Fourier transform,
!    of a complex N by N matrix of data.
!
!    The output is unscaled, that is, a call to CFFTB_2D followed by a call 
!    to CFFTF_2D will return the original data multiplied by N*N.
!
!    For some applications it is desirable to have the transform scaled so
!    the center of the N by N frequency square corresponds to zero
!    frequency.  The user can do this replacing the original input data
!    F(I,J) by F(I,J) * (-1.)**(I+J),  I,J =0,...,N-1.
!
!    Before calling CFFTF_2D or CFFTB_2D, it is necessary to initialize
!    the array WSAVE by calling CFFTI.
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!  Modified:
!
!    12 March 2001
!
!  Parameters:
!
!    Input, integer(4) LDF, the leading dimension of the matrix.
!
!    Input, integer(4) N, the number of rows and columns in the matrix.
!
!    Input/output, complex F(LDF,N),
!    On input, an N by N array of complex values to be transformed.
!    On output, the transformed values.
!
!    Input, real WSAVE(4*N+15), a work array whose values depend on N,
!    and which must be initialized by calling CFFTI.
!
  integer(4) ldf
  integer(4) n
!
  complex(8) f(ldf,n)
  integer(4) i
  real(8)wsave(4*n+15)
!
!  Row transforms:
!
  f(1:n,1:n) = transpose ( f(1:n,1:n) )

  do i = 1, n
    call cfftb ( n, f(1,i), wsave )
  end do

  f(1:n,1:n) = transpose ( f(1:n,1:n) )
!
!  Column transforms:
!
  do i = 1, n
    call cfftb ( n, f(1,i), wsave )
  end do

  return
end
subroutine cfftf ( n, c, wsave )
!
!*******************************************************************************
!
!! CFFTF computes the forward complex discrete Fourier transform.
!
!
!  Discussion:
!
!    This process is sometimes called Fourier analysis.
!
!    CFFTF computes the Fourier coefficients of a complex periodic sequence.
!
!    The transform is not normalized.  To obtain a normalized transform,
!    the output must be divided by N.  Otherwise a call of CFFTF
!    followed by a call of CFFTB will multiply the sequence by N.
!
!    The array WSAVE must be initialized by calling CFFTI.
!
!    The transform is defined by:
!
!      C_out(J) = sum ( 1 <= K <= N ) 
!        C_in(K) * exp ( - sqrt ( -1 ) * ( J - 1 ) * ( K - 1 ) * 2 * PI / N )
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!    P N Swarztrauber, 
!    Vectorizing the FFT's, 
!    in Parallel Computations (G. Rodrigue, editor), 
!    Academic Press, 1982, pages 51-83.
!
!    B L Buzbee, 
!    The SLATEC Common Math Library, 
!    in Sources and Development of Mathematical Software (W. Cowell, editor),
!    Prentice Hall, 1984, pages 302-318.
!
!  Parameters:
!
!    Input, integer(4) N, the length of the sequence to be transformed.  
!    The method is more efficient when N is the product of small primes.
!
!    Input/output, complex C(N).
!    On input, the data sequence to be transformed.
!    On output, the Fourier coefficients.
!
!    Input, real WSAVE(4*N+15).  The array must be initialized by calling 
!    CFFTI.  A different WSAVE array must be used for each different
!    value of N. 
!
  integer(4) n
!
  complex(8) c(n)
  real(8)wsave(4*n+15)
!
  if ( n <= 1 ) then
    return
  end if

  call cfftf1 ( n, c, wsave(1), wsave(2*n+1), wsave(4*n+1) )

  return
end
subroutine cfftf1 ( n, c, ch, wa, ifac )
!
!*******************************************************************************
!
!! CFFTF1 is a lower level routine used by CFFTF.
!
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
!    Input, integer(4) N, the length of the sequence to be transformed.  
!
!    Input/output, complex C(N).
!    On input, the data sequence to be transformed.
!    On output, the Fourier coefficients.
!
!    Input, complex CH(N).
!
!    Input, real WA(2*N).
!
!    Input, integer(4) IFAC(15).
!    IFAC(1) = N, the number that was factored.
!    IFAC(2) = NF, the number of factors.
!    IFAC(3:2+NF), the factors.
!
  integer(4) n
!
  complex(8) c(n)
  complex(8) ch(n)
  integer(4) idl1
  integer(4) ido
  integer(4) ifac(15)
  integer(4) ip
  integer(4) iw
  integer(4) ix2
  integer(4) ix3
  integer(4) ix4
  integer(4) k1
  integer(4) l1
  integer(4) l2
  integer(4) na
  integer(4) nac
  integer(4) nf
  real(8)wa(2*n)
!
  nf = ifac(2)
  na = 0
  l1 = 1
  iw = 1

  do k1 = 1, nf

    ip = ifac(k1+2)
    l2 = ip * l1
    ido = n / l2
    idl1 = 2 * ido * l1

    if ( ip == 4 ) then

      ix2 = iw + 2 * ido
      ix3 = ix2 + 2 * ido
 
      if ( na == 0 ) then
        call passf4 ( 2*ido, l1, c, ch, wa(iw), wa(ix2), wa(ix3) )
      else
        call passf4 ( 2*ido, l1, ch, c, wa(iw), wa(ix2), wa(ix3) )
      end if

      na = 1 - na

    else if ( ip == 2 ) then

      if ( na == 0 ) then
        call passf2 ( 2*ido, l1, c, ch, wa(iw) )
      else
        call passf2 ( 2*ido, l1, ch, c, wa(iw) )
      end if

      na = 1 - na

    else if ( ip == 3 ) then

      ix2 = iw + 2 * ido

      if ( na == 0 ) then
        call passf3 ( 2*ido, l1, c, ch, wa(iw), wa(ix2) )
      else
        call passf3 ( 2*ido, l1, ch, c, wa(iw), wa(ix2) )
      end if

      na = 1 - na

    else if ( ip == 5 ) then

      ix2 = iw + 2 * ido
      ix3 = ix2 + 2 * ido
      ix4 = ix3 + 2 * ido

      if ( na == 0 ) then
        call passf5 ( 2*ido, l1, c, ch, wa(iw), wa(ix2), wa(ix3), wa(ix4) )
      else
        call passf5 ( 2*ido, l1, ch, c, wa(iw), wa(ix2), wa(ix3), wa(ix4) )
      end if

      na = 1 - na

    else

      if ( na == 0 ) then
        call passf ( nac, 2*ido, ip, l1, idl1, c, c, c, ch, ch, wa(iw) )
      else
        call passf ( nac, 2*ido, ip, l1, idl1, ch, ch, ch, c, c, wa(iw) )
      end if

      if ( nac /= 0 ) then
        na = 1 - na
      end if

    end if

    l1 = l2
    iw = iw + ( ip - 1 ) * 2 * ido

  end do

  if ( na /= 0 ) then
    c(1:n) = ch(1:n)
  end if

  return
end
subroutine cfftf_2d ( ldf, n, f, wsave )
!
!*******************************************************************************
!
!! CFFTF_2D computes a two dimensional complex fast Fourier transform.
!
!
!  Discussion:
!
!    The routine computes the forward two dimensional fast Fourier transform,
!    of a complex N by N matrix of data.
!
!    The output is unscaled, that is, a call to CFFTF_2D,
!    followed by a call to CFFTB_2D will return the original data
!    multiplied by N*N.
!
!    For some applications it is desirable to have the transform scaled so
!    the center of the N by N frequency square corresponds to zero
!    frequency.  The user can do this replacing the original input data
!    F(I,J) by F(I,J) *(-1.)**(I+J),  I,J =0,...,N-1.
!
!    Before calling CFFTF_2D or CFFTB_2D, it is necessary to initialize
!    the array WSAVE by calling CFFTI.
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!  Modified:
!
!    12 March 2001
!
!  Parameters:
!
!    Input, integer(4) LDF, the leading dimension of the matrix.
!
!    Input, integer(4) N, the number of rows and columns in the matrix.
!
!    Input/output, complex F(LDF,N),
!    On input, an N by N array of complex values to be transformed.
!    On output, the transformed values.
!
!    Input, real WSAVE(4*N+15), a work array whose values depend on N,
!    and which must be initialized by calling CFFTI.
!
  integer(4) ldf
  integer(4) n
!
  complex(8) f(ldf,n)
  integer(4) i
  real(8)wsave(4*n+15)
!
!  Row transforms:
!
  f(1:n,1:n) = transpose ( f(1:n,1:n) )

  do i = 1, n
    call cfftf ( n, f(1,i), wsave )
  end do

  f(1:n,1:n) = transpose ( f(1:n,1:n) )
!
!  Column transforms:
!
  do i = 1, n
    call cfftf ( n, f(1,i), wsave )
  end do

  return
end
subroutine cffti ( n, wsave )
!
!*******************************************************************************
!
!! CFFTI initializes WSAVE, used in CFFTF and CFFTB. 
!
!
!  Discussion:
!
!    The prime factorization of N together with a tabulation of the 
!    trigonometric functions are computed and stored in WSAVE.
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!    P N Swarztrauber, 
!    Vectorizing the FFT's, 
!    in Parallel Computations (G. Rodrigue, editor), 
!    Academic Press, 1982, pages 51-83.
!
!    B L Buzbee, 
!    The SLATEC Common Math Library, 
!    in Sources and Development of Mathematical Software (W. Cowell, editor),
!    Prentice Hall, 1984, pages 302-318.
!
!  Parameters:
!
!    Input, integer(4) N, the length of the sequence to be transformed.
!
!    Output, real WSAVE(4*N+15), contains data, dependent on the value
!    of N, which is necessary for the CFFTF or CFFTB routines.  
!
  integer(4) n
!
  real(8)wsave(4*n+15)
!
  if ( n <= 1 ) then
    return
  end if

  call cffti1 ( n, wsave(2*n+1), wsave(4*n+1) )

  return
end
subroutine cffti1 ( n, wa, ifac )
!
!*******************************************************************************
!
!! CFFTI1 is a lower level routine used by CFFTI.
!
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
!    Input, integer(4) N, the length of the sequence to be transformed.
!
!    Input, real WA(2*N).
!
!    Input, integer(4) IFAC(15).
!    IFAC(1) = N, the number that was factored.
!    IFAC(2) = NF, the number of factors.
!    IFAC(3:2+NF), the factors.
!
  integer(4) n
!
  real(8) ::arg
  real(8) ::argh
  real(8) argld
  real(8) fi
  integer(4) i
  integer(4) i1
  integer(4) ib
  integer(4) ido
  integer(4) ifac(15)
  integer(4) ii
  integer(4) ip
  integer(4) j
  integer(4) k1
  integer(4) l1
  integer(4) l2
  integer(4) ld
  integer(4) nf
  real(8) pimach
  real(8) wa(2*n)
!
  call i_factor ( n, ifac )

  nf = ifac(2)

  argh = 2.0D+00 * pimach() / dble ( n )
  i = 2
  l1 = 1

  do k1 = 1, nf

    ip = ifac(k1+2)
    ld = 0
    l2 = l1 * ip
    ido = n / l2

    do j = 1, ip-1

      i1 = i
      wa(i-1) = 1.0D+00
      wa(i) = 0.0D+00
      ld = ld + l1
      fi = 0.0D+00
      argld = dble ( ld ) * argh

      do ii = 4, 2*ido+2, 2
        i = i + 2
        fi = fi + 1.0D+00
        arg = fi * argld
        wa(i-1) = cos ( arg )
        wa(i) = sin ( arg )
      end do

      if ( ip > 5 ) then
        wa(i1-1) = wa(i-1)
        wa(i1) = wa(i)
      end if

    end do

    l1 = l2

  end do

  return
end
subroutine cosqb ( n, x, wsave )
!
!*******************************************************************************
!
!! COSQB computes the fast cosine transform of quarter wave data. 
!
!
!  Discussion:
!
!    COSQB computes a sequence from its representation in terms of a cosine 
!    series with odd wave numbers.
!
!    The transform is defined by:
!
!      X_out(I) = sum ( 1 <= K <= N ) 
!
!        4 * X_in(K) * cos ( ( 2 * K - 1 ) * ( I - 1 ) * PI / ( 2 * N ) )
!
!    COSQB is the unnormalized inverse of COSQF since a call of COSQB
!    followed by a call of COSQF will multiply the input sequence X by 4*N.
!
!    The array WSAVE must be initialized by calling COSQI.
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!    P N Swarztrauber, 
!    Vectorizing the FFT's, 
!    in Parallel Computations (G. Rodrigue, editor), 
!    Academic Press, 1982, pages 51-83.
!
!    B L Buzbee, 
!    The SLATEC Common Math Library, 
!    in Sources and Development of Mathematical Software (W. Cowell, editor),
!    Prentice Hall, 1984, pages 302-318.
!
!  Parameters:
!
!    Input, integer(4) N, the length of the array X.  The method is
!    more efficient when N is the product of small primes.
!
!    Input/output, real X(N).
!    On input, the cosine series coefficients.
!    On output, the corresponding data vector.
!
!    Input, real WSAVE(3*N+15), contains data, depending on N, and 
!    required by the algorithm.  The WSAVE array must be initialized by 
!    calling COSQI.  A different WSAVE array must be used for each different
!    value of N. 
!
  integer(4) n
!
  real(8), parameter :: tsqrt2 = 2.82842712474619D+00
  real(8) wsave(3*n+15)
  real(8) x(n)
  real(8) x1
!
  if ( n < 2 ) then
    x(1) = 4.0D+00 * x(1)
  else if ( n == 2 ) then
    x1 = 4.0D+00 * ( x(1) + x(2) )
    x(2) = tsqrt2 * ( x(1) - x(2) )
    x(1) = x1
  else
    call cosqb1 ( n, x, wsave(1), wsave(n+1) )
  end if

  return
end
subroutine cosqb1 ( n, x, w, xh )
!
!*******************************************************************************
!
!! COSQB1 is a lower level routine used by COSQB.
!
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
!    Input, integer(4) N, the length of the array.  
!
!    Input/output, real X(N).
!    On input, the cosine series coefficients.
!    On output, the corresponding data vector.
!
!    Input, real W(N).
!
!    Input, real XH(2*N+15).
!
  integer(4) n
!
  integer(4) i
  integer(4) k
  integer(4) kc
  integer(4) ns2
  real(8) w(n)
  real(8) x(n)
  real(8) xh(2*n+15)
  real(8) xim1
!
  ns2 = ( n + 1 ) / 2

  do i = 3, n, 2
    xim1 = x(i-1) + x(i)
    x(i) = x(i) - x(i-1)
    x(i-1) = xim1
  end do

  x(1) = x(1) + x(1)

  if ( mod ( n, 2 ) == 0 ) then
    x(n) = 2.0D+00 * x(n)
  end if

  call rfftb ( n, x, xh )

  do k = 2, ns2
    kc = n + 2 - k
    xh(k) = w(k-1) * x(kc) + w(kc-1) * x(k)
    xh(kc) = w(k-1) * x(k) - w(kc-1) * x(kc)
  end do

  if ( mod ( n, 2 ) == 0 ) then
    x(ns2+1) = w(ns2) * ( x(ns2+1) + x(ns2+1) )
  end if

  do k = 2, ns2
    kc = n + 2 - k
    x(k) = xh(k) + xh(kc)
    x(kc) = xh(k) - xh(kc)
  end do

  x(1) = 2.0D+00 * x(1)

  return
end
subroutine cosqf ( n, x, wsave )
!
!*******************************************************************************
!
!! COSQF computes the fast cosine transform of quarter wave data. 
!
!
!  Discussion:
!
!    COSQF computes the coefficients in a cosine series representation 
!    with only odd wave numbers. 
!
!    COSQF is the unnormalized inverse of COSQB since a call of COSQF
!    followed by a call of COSQB will multiply the input sequence X
!    by 4*N.
!
!    The array WSAVE must be initialized by calling COSQI.
!
!    The transform is defined by:
!
!      X_out(I) = X_in(1) + sum ( 2 <= K <= N )
!
!        2 * X_in(K) * cos ( ( 2 * I - 1 ) * ( K - 1 ) * PI / ( 2 * N ) )
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!    P N Swarztrauber, 
!    Vectorizing the FFT's, 
!    in Parallel Computations (G. Rodrigue, editor), 
!    Academic Press, 1982, pages 51-83.
!
!    B L Buzbee, 
!    The SLATEC Common Math Library, 
!    in Sources and Development of Mathematical Software (W. Cowell, editor),
!    Prentice Hall, 1984, pages 302-318.
!
!  Parameters:
!
!    Input, integer(4) N, the length of the array X.  The method is
!    more efficient when N is the product of small primes.
!
!    Input/output, real(8) X(N).
!    On input, the data to be transformed.
!    On output, the transformed data.
!
!    Input, real(8) WSAVE(3*N+15), contains data, depending on N, and 
!    required by the algorithm.  The WSAVE array must be initialized by 
!    calling COSQI.  A different WSAVE array must be used for each different
!    value of N. 
!
  integer(4) n
!
  real(8), parameter :: sqrt2 = 1.4142135623731D+00
  real(8) tsqx
  real(8) wsave(3*n+15)
  real(8) x(n)
!
  if ( n < 2 ) then

  else if ( n == 2 ) then
    tsqx = sqrt2 * x(2)
    x(2) = x(1) - tsqx
    x(1) = x(1) + tsqx
  else
    call cosqf1 ( n, x, wsave(1), wsave(n+1) )
  end if

  return
end
subroutine cosqf1 ( n, x, w, xh )
!
!*******************************************************************************
!
!! COSQF1 is a lower level routine used by COSQF.
!
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
!    Input, integer(4) N, the length of the array to be transformed.  
!
!    Input/output, real(8) X(N).
!    On input, the data to be transformed.
!    On output, the transformed data.
!
!    Input, real(8) W(N).
!
!    Input, real(8) XH(2*N+15).
!
  integer(4) n
!
  integer(4) i
  integer(4) k
  integer(4) kc
  integer(4) ns2
  real(8) w(n)
  real(8) x(n)
  real(8) xh(2*n+15)
  real(8) xim1
!
  ns2 = ( n + 1 ) / 2

  do k = 2, ns2
    kc = n + 2 - k
    xh(k) = x(k) + x(kc)
    xh(kc) = x(k) - x(kc)
  end do

  if ( mod ( n, 2 ) == 0 ) then
    xh(ns2+1) = x(ns2+1) + x(ns2+1)
  end if

  do k = 2, ns2
    kc = n+2-k
    x(k) = w(k-1) * xh(kc) + w(kc-1) * xh(k)
    x(kc) = w(k-1) * xh(k) - w(kc-1) * xh(kc)
  end do

  if ( mod ( n, 2 ) == 0 ) then
    x(ns2+1) = w(ns2) * xh(ns2+1)
  end if

  call rfftf ( n, x, xh )

  do i = 3, n, 2
    xim1 = x(i-1) - x(i)
    x(i) = x(i-1) + x(i)
    x(i-1) = xim1
  end do

  return
end
subroutine cosqi ( n, wsave )
!
!*******************************************************************************
!
!! COSQI initializes WSAVE, used in COSQF and COSQB. 
!
!
!  Discussion:
!
!    The prime factorization of N together with a tabulation of the 
!    trigonometric functions are computed and stored in WSAVE.
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!    P N Swarztrauber, 
!    Vectorizing the FFT's, 
!    in Parallel Computations (G. Rodrigue, editor), 
!    Academic Press, 1982, pages 51-83.
!
!    B L Buzbee, 
!    The SLATEC Common Math Library, 
!    in Sources and Development of Mathematical Software (W. Cowell, editor),
!    Prentice Hall, 1984, pages 302-318.
!
!  Parameters:
!
!    Input, integer(4) N, the length of the array to be transformed.  The method 
!    is more efficient when N is the product of small primes.
!
!    Output, real(8) WSAVE(3*N+15), contains data, depending on N, and 
!    required by the COSQB and COSQF algorithms.  
!
  integer(4) n
!
  real(8) dt
  integer(4) k
  real(8) pimach
  real(8) wsave(3*n+15)
!
  dt = 0.5D+00 * pimach() / dble ( n )

  do k = 1, n
    wsave(k) = cos ( dble ( k ) * dt )
  end do

  call rffti ( n, wsave(n+1) )

  return
end
subroutine cost ( n, x, wsave )
!
!*******************************************************************************
!
!! COST computes the discrete Fourier cosine transform of an even sequence. 
!
!
!  Discussion:
!
!    COST is the unnormalized inverse of itself since a call of COST
!    followed by another call of COST will multiply the input sequence
!    X by 2*(N-1). 
!
!    The array WSAVE must be initialized by calling COSTI.
!
!    The transform is defined by:
!
!      X_out(I) = X_in(1) + (-1) **(I-1) * X_in(N) + sum ( 2 <= K <= N-1 )
!
!        2 * X_in(K) * cos ( ( K - 1 ) * ( I - 1 ) * PI / ( N - 1 ) )
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!    P N Swarztrauber, 
!    Vectorizing the FFT's, 
!    in Parallel Computations (G. Rodrigue, editor), 
!    Academic Press, 1982, pages 51-83.
!
!    B L Buzbee, 
!    The SLATEC Common Math Library, 
!    in Sources and Development of Mathematical Software (W. Cowell, editor),
!    Prentice Hall, 1984, pages 302-318.
!
!  Parameters:
!
!    Input, integer(4) N, the length of the sequence to be transformed.  The 
!    method is more efficient when N-1 is the product of small primes.
!
!    Input/output, real(8) X(N).
!    On input, the sequence to be transformed.
!    On output, the transformed sequence.
!
!    Input, real(8) WSAVE(3*N+15).
!    The WSAVE array must be initialized by calling COSTI.  A different 
!    array must be used for each different value of N. 
!
  integer(4) n
!
  real(8) c1
  integer(4) i
  integer(4) k
  integer(4) kc
  integer(4) ns2
  real(8) t1
  real(8) t2
  real(8) tx2
  real(8) wsave(3*n+15)
  real(8) x(n)
  real(8) x1h
  real(8) x1p3
  real(8) xi
  real(8) xim2
!
  ns2 = n / 2

  if ( n <= 1 ) then
    return
  end if

  if ( n == 2 ) then
    x1h = x(1) + x(2)
    x(2) = x(1) - x(2)
    x(1) = x1h
    return
  end if

  if ( n == 3 ) then
    x1p3 = x(1) + x(3)
    tx2 = x(2) + x(2)
    x(2) = x(1) - x(3)
    x(1) = x1p3 + tx2
    x(3) = x1p3 - tx2
    return
  end if

  c1 = x(1) - x(n)
  x(1) = x(1) + x(n)

  do k = 2, ns2
    kc = n + 1 - k
    t1 = x(k) + x(kc)
    t2 = x(k) - x(kc)
    c1 = c1 + wsave(kc) * t2
    t2 = wsave(k) * t2
    x(k) = t1 - t2
    x(kc) = t1 + t2
  end do

  if ( mod ( n, 2 ) /= 0 ) then
    x(ns2+1) = x(ns2+1) + x(ns2+1)
  end if

  call rfftf ( n-1, x, wsave(n+1) )

  xim2 = x(2)
  x(2) = c1

  do i = 4, n, 2
    xi = x(i)
    x(i) = x(i-2) - x(i-1)
    x(i-1) = xim2
    xim2 = xi
  end do

  if ( mod ( n, 2 ) /= 0 ) then
    x(n) = xim2
  end if

  return
end
subroutine costi ( n, wsave )
!
!*******************************************************************************
!
!! COSTI initializes WSAVE, used in COST.
!
!
!  Discussion:
!
!    The prime factorization of N together with a tabulation of the 
!    trigonometric functions are computed and stored in WSAVE.
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!    P N Swarztrauber, 
!    Vectorizing the FFT's, 
!    in Parallel Computations (G. Rodrigue, editor), 
!    Academic Press, 1982, pages 51-83.
!
!    B L Buzbee, 
!    The SLATEC Common Math Library, 
!    in Sources and Development of Mathematical Software (W. Cowell, editor),
!    Prentice Hall, 1984, pages 302-318.
!
!  Parameters:
!
!    Input, integer(4) N, the length of the sequence to be transformed.  The 
!    method is more efficient when N-1 is the product of small primes.
!
!    Output, real WSAVE(3*N+15), contains data, depending on N, and 
!    required by the COST algorithm.
!
  integer(4) n
!
  real(8) dt
  integer(4) k
  real(8) pimach
  real(8) wsave(3*n+15)
!
  if ( n <= 3 ) then
    return
  end if

  dt = pimach ( ) / dble ( n - 1 )

  do k = 2, ( n / 2 )
    wsave(k)     = 2.0D+00 * sin ( dble ( k - 1 ) * dt )
    wsave(n+1-k) = 2.0D+00 * cos ( dble ( k - 1 ) * dt )
  end do

  call rffti ( n-1, wsave(n+1) )

  return
end
subroutine cvec_random ( alo, ahi, n, a )
!
!*******************************************************************************
!
!! CVEC_RANDOM returns a random complex vector in a given range.
!
!
!  Modified:
!
!    08 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ALO, AHI, the range allowed for the entries.
!
!    Input, integer(4) N, the number of entries in the vector.
!
!    Output, complex A(N), the vector of randomly chosen values.
!
  integer(4) n
!
  complex(8) a(n)
  real(8) ahi
  real(8) ai
  real(8) alo
  real(8) ar
  integer(4) i
!
  do i = 1, n

    call r_random ( alo, ahi, ar )
    call r_random ( alo, ahi, ai )

    a(i) = cmplx ( ar, ai )

  end do

  return
end
subroutine ezfftb ( n, r, azero, a, b, wsave )
!
!*******************************************************************************
!
!! EZFFTB computes a real periodic sequence from its Fourier coefficients.
!
!
!  Discussion:
!
!    This process is sometimes called Fourier synthesis.
!
!    EZFFTB is a simplified but slower version of RFFTB.
!
!    The transform is defined by: 
!
!      R(I) = AZERO + sum ( 1 <= K <= N/2 )
!
!          A(K) * cos ( K * ( I - 1 ) * 2 * PI / N ) 
!        + B(K) * sin ( K * ( I - 1 ) * 2 * PI / N )
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!    P N Swarztrauber, 
!    Vectorizing the FFT's, 
!    in Parallel Computations (G. Rodrigue, editor), 
!    Academic Press, 1982, pages 51-83.
!
!    B L Buzbee, 
!    The SLATEC Common Math Library, 
!    in Sources and Development of Mathematical Software (W. Cowell, editor),
!    Prentice Hall, 1984, pages 302-318.
!
!  Parameters:
!
!    Input, integer(4) N, the length of the output array.  The 
!    method is more efficient when N is the product of small primes.
!
!    Output, real R(N), the reconstructed data sequence.
!
!    Input, real AZERO, the constant Fourier coefficient.
!
!    Input, real A(N/2), B(N/2), the Fourier coefficients.
!
!    Input, real WSAVE(3*N+15), a work array.  The WSAVE array must be
!    initialized by calling EZFFFTI.  A different WSAVE array must be used 
!    for each different value of N. 
!
  integer(4) n
!
  real(8) a(n/2)
  real(8) azero
  real(8) b(n/2)
  integer(4) i
  integer(4) ns2
  real(8) r(n)
  real(8) wsave(3*n+15)
!
  if ( n < 2 ) then

    r(1) = azero

  else if ( n == 2 ) then

    r(1) = azero + a(1)
    r(2) = azero - a(1)

  else

    ns2 = ( n - 1 ) / 2

    do i = 1, ns2
      r(2*i) = 0.5D+00 * a(i)
      r(2*i+1) = -0.5D+00 * b(i)
    end do

    r(1) = azero

    if ( mod ( n, 2 ) == 0 ) then
      r(n) = a(ns2+1)
    end if

    call rfftb ( n, r, wsave(n+1) )

  end if

  return
end
subroutine ezfftf ( n, r, azero, a, b, wsave )
!
!*******************************************************************************
!
!! EZFFTF computes the Fourier coefficients of a real periodic sequence.
!
!
!  Discussion:
!
!    This process is sometimes called Fourier analysis.
!
!    EZFFTF is a simplified but slower version of RFFTF.
!
!    The transform is defined by:
!
!      AZERO = sum ( 1 <= I <= N ) R(I) / N,
!
!    and, for K = 1 to (N-1)/2,
!
!      A(K) = sum ( 1 <= I <= N )
!        ( 2 / N ) * R(I) * cos ( K * ( I - 1 ) * 2 * PI / N )
!
!    and, if N is even, then
!
!      A(N/2) = sum ( 1 <= I <= N ) (-1) **(I-1) * R(I) / N
!
!    For K = 1 to (N-1)/2,
!
!      B(K) = sum ( 1 <= I <= N )
!        ( 2 / N ) * R(I) * sin ( K * ( I - 1 ) * 2 * PI / N )
!
!    and, if N is even, then
!
!      B(N/2) = 0.
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!    P N Swarztrauber, 
!    Vectorizing the FFT's, 
!    in Parallel Computations (G. Rodrigue, editor), 
!    Academic Press, 1982, pages 51-83.
!
!    B L Buzbee, 
!    The SLATEC Common Math Library, 
!    in Sources and Development of Mathematical Software (W. Cowell, editor),
!    Prentice Hall, 1984, pages 302-318.
!
!  Parameters:
!
!    Input, integer(4) N, the length of the array to be transformed.  The 
!    method is more efficient when N is the product of small primes.
!
!    Input, real R(N), the sequence to be transformed.
!
!    Input, real WSAVE(3*N+15), a work array.  The WSAVE array must be
!    initialized by calling EZFFTI.  A different WSAVE array must be used 
!    for each different value of N. 
!
!    Output, real AZERO, the constant Fourier coefficient.
!
!    Output, real A(N/2), B(N/2), the Fourier coefficients.
!
  integer(4) n
!
  real(8) a(n/2)
  real(8) azero
  real(8) b(n/2)
  real(8) cf
  integer(4) i
  integer(4) ns2
  real(8) r(n)
  real(8) wsave(3*n+15)
!
  if ( n < 2 ) then

    azero = r(1)

  else if ( n == 2 ) then

    azero = 0.5D+00 * ( r(1) + r(2) )
    a(1) = 0.5D+00 * ( r(1) - r(2) )

  else

    wsave(1:n) = r(1:n)

    call rfftf ( n, wsave(1), wsave(n+1) )

    cf = 2.0D+00 / dble ( n )
    azero = 0.5D+00 * cf * wsave(1)
    ns2 = ( n + 1 ) / 2

    do i = 1, ns2-1
      a(i) = cf * wsave(2*i)
      b(i) = -cf * wsave(2*i+1)
    end do

    if ( mod ( n, 2 ) /= 1 ) then
      a(ns2) = 0.5D+00 * cf * wsave(n)
      b(ns2) = 0.0D+00
    end if

  end if

  return
end
subroutine ezffti ( n, wsave )
!
!*******************************************************************************
!
!! EZFFTI initializes WSAVE, used in EZFFTF and EZFFTB. 
!
!
!  Discussion:
!
!    The prime factorization of N together with a tabulation of the 
!    trigonometric functions are computed and stored in WSAVE.
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!    P N Swarztrauber, 
!    Vectorizing the FFT's, 
!    in Parallel Computations (G. Rodrigue, editor), 
!    Academic Press, 1982, pages 51-83.
!
!    B L Buzbee, 
!    The SLATEC Common Math Library, 
!    in Sources and Development of Mathematical Software (W. Cowell, editor),
!    Prentice Hall, 1984, pages 302-318.
!
!  Parameters:
!
!    Input, integer(4) N, the length of the array to be transformed.  The 
!    method is more efficient when N is the product of small primes.
!
!    Output, real WSAVE(3*N+15), contains data, dependent on the value
!    of N, which is necessary for the EZFFTF or EZFFTB routines.  
!
  integer(4) n
!
  real(8) wsave(3*n+15)
!
  if ( n <= 1 ) then
    return
  end if

  call ezffti1 ( n, wsave(2*n+1), wsave(3*n+1) )

  return
end
subroutine ezffti1 ( n, wa, ifac )
!
!*******************************************************************************
!
!! EZFFTI1 is a lower level routine used by EZFFTI.
!
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
!    Input, integer(4) N, the length of the array to be transformed. 
!
!    Output, real WA(N).
!
!    Input, integer(4) IFAC(15).
!    IFAC(1) = N, the number that was factored.
!    IFAC(2) = NF, the number of factors.
!    IFAC(3:2+NF), the factors.
!
  integer(4) n
!
  real(8) arg1
  real(8) argh
  real(8) ch1
  real(8) ch1h
  real(8) dch1
  real(8) dsh1
  integer(4) i
  integer(4) ib
  integer(4) ido
  integer(4) ifac(15)
  integer(4) ii
  integer(4) ip
  integer(4) is
  integer(4) j
  integer(4) k1
  integer(4) l1
  integer(4) l2
  integer(4) nf
  real(8) pimach
  real(8) sh1
  real(8) wa(n)
!
  call i_factor ( n, ifac )

  nf = ifac(2)

  argh = 2.0D+00 * pimach() / dble ( n )
  is = 0
  l1 = 1

  do k1 = 1, nf-1

    ip = ifac(k1+2)
    l2 = l1 * ip
    ido = n / l2
    arg1 = dble ( l1 ) * argh
    ch1 = 1.0D+00
    sh1 = 0.0D+00
    dch1 = cos ( arg1 )
    dsh1 = sin ( arg1 )

    do j = 1, ip-1

      ch1h = dch1 * ch1 - dsh1 * sh1
      sh1  = dch1 * sh1 + dsh1 * ch1
      ch1 = ch1h
      i = is + 2
      wa(i-1) = ch1
      wa(i) = sh1

      do ii = 5, ido, 2
        i = i + 2
        wa(i-1) = ch1 * wa(i-3) - sh1 * wa(i-2)
        wa(i)   = ch1 * wa(i-2) + sh1 * wa(i-3)
      end do

      is = is + ido

    end do

    l1 = l2

  end do

  return
end
subroutine i_factor ( n, ifac )
!
!*******************************************************************************
!
!! I_FACTOR factors an integer.
!
!
!  Modified:
!
!    14 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
!    Input, integer(4) N, the number to be factored.
!
!    Output, integer(4) IFAC(15).
!    IFAC(1) = N, the number that was factored.
!    IFAC(2) = NF, the number of factors.
!    IFAC(3:2+NF), the factors.
!
  integer(4) i
  integer(4) ib
  integer(4) ifac(15)
  integer(4) j
  integer(4) n
  integer(4) nf
  integer(4) nl
  integer(4) nq
  integer(4) nr
  integer(4) ntry
!
  ifac(1) = n

  nf = 0
  nl = n

  if ( n == 0 ) then
    nf = 1
    ifac(2) = nf
    ifac(2+nf) = 0
    return
  end if

  if ( n < 1 ) then
    nf = nf + 1
    ifac(2+nf) = -1
    nl = - n
  end if

  if ( nl == 1 ) then
    nf = nf + 1
    ifac(2) = nf
    ifac(2+nf) = 1
    return
  end if

  j = 0

  do while ( nl > 1 )

    j = j + 1
!
!  Choose a trial divisor, NTRY.
!
    if ( j == 1 ) then
      ntry = 4
    else if ( j == 2 ) then
      ntry = 2
    else if ( j == 3 ) then
      ntry = 3
    else if ( j == 4 ) then
      ntry = 5
    else
      ntry = ntry + 2
    end if
!
!  Divide by the divisor as many times as possible.
!
    do

      nq = nl / ntry
      nr = nl - ntry * nq

      if ( nr /= 0 ) then
        exit
      end if

      nl = nq
      nf = nf + 1
!
!  Make sure factors of 2 appear in the front of the list.
!
      if ( ntry /= 2 ) then

        ifac(2+nf) = ntry

      else

        do i = nf, 2, -1
          ifac(i+2) = ifac(i+1)
        end do
        ifac(3) = 2

      end if

    end do

  end do

  ifac(2) = nf

  return
end
subroutine passb ( nac, ido, ip, l1, idl1, cc, c1, c2, ch, ch2, wa )
!
!*******************************************************************************
!
!! PASSB is a lower level routine used by CFFTB1.
!
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
  integer(4) idl1
  integer(4) ido
  integer(4) ip
  integer(4) l1
!
  real(8) c1(ido,l1,ip)
  real(8) c2(idl1,ip)
  real(8) cc(ido,ip,l1)
  real(8) ch(ido,l1,ip)
  real(8) ch2(idl1,ip)
  integer(4) i
  integer(4) idij
  integer(4) idj
  integer(4) idl
  integer(4) idlj
  integer(4) idp
  integer(4) ik
  integer(4) inc
  integer(4) ipph
  integer(4) j
  integer(4) jc
  integer(4) k
  integer(4) l
  integer(4) lc
  integer(4) nac
  integer(4) nt
  real(8) wa(*)
  real(8) wai
  real(8) war
!
  nt = ip * idl1
  ipph = ( ip + 1 ) / 2
  idp = ip * ido

  if ( ido >= l1 ) then

    do j = 2, ipph
      jc = ip + 2 - j
      do k = 1, l1
        ch(1:ido,k,j)  = cc(1:ido,j,k) + cc(1:ido,jc,k)
        ch(1:ido,k,jc) = cc(1:ido,j,k) - cc(1:ido,jc,k)
      end do
    end do

    ch(1:ido,1:l1,1) = cc(1:ido,1,1:l1)

  else

    do j = 2, ipph
      jc = ip + 2 - j
      do i = 1, ido
        ch(i,1:l1,j)  = cc(i,j,1:l1) + cc(i,jc,1:l1)
        ch(i,1:l1,jc) = cc(i,j,1:l1) - cc(i,jc,1:l1)
      end do
    end do

    ch(1:ido,1:l1,1) = cc(1:ido,1,1:l1)

  end if

  idl = 2 - ido
  inc = 0

  do l = 2, ipph

    lc = ip + 2 - l
    idl = idl + ido

    do ik = 1, idl1
      c2(ik,l) = ch2(ik,1) + wa(idl-1) * ch2(ik,2)
      c2(ik,lc) =            wa(idl)   * ch2(ik,ip)
    end do

    idlj = idl
    inc = inc + ido

    do j = 3, ipph

      jc = ip + 2 - j
      idlj = idlj + inc
      if ( idlj > idp ) then
        idlj = idlj - idp
      end if

      war = wa(idlj-1)
      wai = wa(idlj)

      do ik = 1, idl1
        c2(ik,l)  = c2(ik,l)  + war * ch2(ik,j)
        c2(ik,lc) = c2(ik,lc) + wai * ch2(ik,jc)
      end do

    end do

  end do

  do j = 2, ipph
    ch2(1:idl1,1) = ch2(1:idl1,1) + ch2(1:idl1,j)
  end do

  do j = 2, ipph
    jc = ip + 2 - j
    do ik = 2, idl1, 2
      ch2(ik-1,j)  = c2(ik-1,j) - c2(ik,jc)
      ch2(ik-1,jc) = c2(ik-1,j) + c2(ik,jc)
      ch2(ik,j)    = c2(ik,j)   + c2(ik-1,jc)
      ch2(ik,jc)   = c2(ik,j)   - c2(ik-1,jc)
    end do
  end do

  nac = 1

  if ( ido == 2 ) then
    return
  end if

  nac = 0
  c2(1:idl1,1) = ch2(1:idl1,1)
  c1(1:2,1:l1,2:ip) = ch(1:2,1:l1,2:ip)

  if ( ( ido / 2 ) <= l1 ) then

    idij = 0
    do j = 2, ip
      idij = idij + 2
      do i = 4, ido, 2
        idij = idij + 2
        c1(i-1,1:l1,j) = wa(idij-1) * ch(i-1,1:l1,j) - wa(idij) * ch(i,1:l1,j)
        c1(i,1:l1,j)   = wa(idij-1) * ch(i,1:l1,j)   + wa(idij) * ch(i-1,1:l1,j)
      end do
    end do

  else

    idj = 2 - ido

    do j = 2, ip
      idj = idj + ido
      do k = 1, l1
        idij = idj
        do i = 4, ido, 2
          idij = idij + 2
          c1(i-1,k,j) = wa(idij-1) * ch(i-1,k,j) - wa(idij) * ch(i,k,j)
          c1(i,k,j)   = wa(idij-1) * ch(i,k,j)   + wa(idij) * ch(i-1,k,j)
        end do
      end do
    end do

  end if

  return
end
subroutine passb2 ( ido, l1, cc, ch, wa1 )
!
!*******************************************************************************
!
!! PASSB2 is a lower level routine used by CFFTB1.
!
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
  integer(4) ido
  integer(4) l1
!
  real(8) cc(ido,2,l1)
  real(8) ch(ido,l1,2)
  integer(4) i
  integer(4) k
  real(8) ti2
  real(8) tr2
  real(8) wa1(ido)
!
  if ( ido <= 2 ) then

    ch(1,1:l1,1) = cc(1,1,1:l1) + cc(1,2,1:l1)
    ch(1,1:l1,2) = cc(1,1,1:l1) - cc(1,2,1:l1)
    ch(2,1:l1,1) = cc(2,1,1:l1) + cc(2,2,1:l1)
    ch(2,1:l1,2) = cc(2,1,1:l1) - cc(2,2,1:l1)

  else

    do k = 1, l1
      do i = 2, ido, 2

        ch(i-1,k,1) = cc(i-1,1,k) + cc(i-1,2,k)
        tr2         = cc(i-1,1,k) - cc(i-1,2,k)
        ch(i,k,1)   = cc(i,1,k)   + cc(i,2,k)
        ti2         = cc(i,1,k)   - cc(i,2,k)

        ch(i,k,2)   = wa1(i-1) * ti2 + wa1(i) * tr2
        ch(i-1,k,2) = wa1(i-1) * tr2 - wa1(i) * ti2

      end do
    end do

  end if

  return
end
subroutine passb3 ( ido, l1, cc, ch, wa1, wa2 )
!
!*******************************************************************************
!
!! PASSB3 is a lower level routine used by CFFTB1.
!
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
  integer(4) ido
  integer(4) l1
!
  real(8) cc(ido,3,l1)
  real(8) ch(ido,l1,3)
  real(8) ci2
  real(8) ci3
  real(8) cr2
  real(8) cr3
  real(8) di2
  real(8) di3
  real(8) dr2
  real(8) dr3
  integer(4) i
  integer(4) k
  real(8), parameter :: taui = 0.866025403784439D+00
  real(8), parameter :: taur = -0.5D+00
  real(8) ti2
  real(8) tr2
  real(8) wa1(ido)
  real(8) wa2(ido)
!
  if ( ido == 2 ) then

    do k = 1, l1

      tr2 = cc(1,2,k) + cc(1,3,k)
      cr2 = cc(1,1,k) + taur * tr2
      ch(1,k,1) = cc(1,1,k) + tr2

      ti2 = cc(2,2,k) + cc(2,3,k)
      ci2 = cc(2,1,k) + taur * ti2
      ch(2,k,1) = cc(2,1,k) + ti2

      cr3 = taui * ( cc(1,2,k) - cc(1,3,k) )
      ci3 = taui * ( cc(2,2,k) - cc(2,3,k) )

      ch(1,k,2) = cr2 - ci3
      ch(1,k,3) = cr2 + ci3
      ch(2,k,2) = ci2 + cr3
      ch(2,k,3) = ci2 - cr3

    end do

  else

    do k = 1, l1
      do i = 2, ido, 2

        tr2 = cc(i-1,2,k) + cc(i-1,3,k)
        cr2 = cc(i-1,1,k) + taur * tr2
        ch(i-1,k,1) = cc(i-1,1,k) + tr2

        ti2 = cc(i,2,k) + cc(i,3,k)
        ci2 = cc(i,1,k) + taur * ti2
        ch(i,k,1) = cc(i,1,k) + ti2

        cr3 = taui * ( cc(i-1,2,k) - cc(i-1,3,k) )
        ci3 = taui * ( cc(i,2,k) - cc(i,3,k) )

        dr2 = cr2 - ci3
        dr3 = cr2 + ci3
        di2 = ci2 + cr3
        di3 = ci2 - cr3

        ch(i,k,2)   = wa1(i-1) * di2 + wa1(i) * dr2
        ch(i-1,k,2) = wa1(i-1) * dr2 - wa1(i) * di2
        ch(i,k,3)   = wa2(i-1) * di3 + wa2(i) * dr3
        ch(i-1,k,3) = wa2(i-1) * dr3 - wa2(i) * di3

      end do
    end do

  end if

  return
end
subroutine passb4 ( ido, l1, cc, ch, wa1, wa2, wa3 )
!
!*******************************************************************************
!
!! PASSB4 is a lower level routine used by CFFTB1.
!
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
  integer(4) ido
  integer(4) l1
!
  real(8) cc(ido,4,l1)
  real(8) ch(ido,l1,4)
  real(8) ci1
  real(8) ci2
  real(8) ci3
  real(8) ci4
  real(8) cr1
  real(8) cr2
  real(8) cr3
  real(8) cr4
  integer(4) i
  integer(4) k
  real(8) ti1
  real(8) ti2
  real(8) ti3
  real(8) ti4
  real(8) tr1
  real(8) tr2
  real(8) tr3
  real(8) tr4
  real(8) wa1(ido)
  real(8) wa2(ido)
  real(8) wa3(ido)
!
  if ( ido == 2 ) then

    do k = 1, l1

      ti1 = cc(2,1,k) - cc(2,3,k)
      ti2 = cc(2,1,k) + cc(2,3,k)
      tr4 = cc(2,4,k) - cc(2,2,k)
      ti3 = cc(2,2,k) + cc(2,4,k)
      tr1 = cc(1,1,k) - cc(1,3,k)
      tr2 = cc(1,1,k) + cc(1,3,k)
      ti4 = cc(1,2,k) - cc(1,4,k)
      tr3 = cc(1,2,k) + cc(1,4,k)

      ch(1,k,1) = tr2 + tr3
      ch(1,k,3) = tr2 - tr3
      ch(2,k,1) = ti2 + ti3
      ch(2,k,3) = ti2 - ti3
      ch(1,k,2) = tr1 + tr4
      ch(1,k,4) = tr1 - tr4
      ch(2,k,2) = ti1 + ti4
      ch(2,k,4) = ti1 - ti4

    end do

  else

    do k = 1, l1
      do i = 2, ido, 2

        ti1 = cc(i,1,k) - cc(i,3,k)
        ti2 = cc(i,1,k) + cc(i,3,k)
        ti3 = cc(i,2,k) + cc(i,4,k)
        tr4 = cc(i,4,k) - cc(i,2,k)

        tr1 = cc(i-1,1,k) - cc(i-1,3,k)
        tr2 = cc(i-1,1,k) + cc(i-1,3,k)
        ti4 = cc(i-1,2,k) - cc(i-1,4,k)
        tr3 = cc(i-1,2,k) + cc(i-1,4,k)

        ch(i-1,k,1) = tr2 + tr3
        cr3 = tr2 - tr3
        ch(i,k,1) = ti2 + ti3
        ci3 = ti2 - ti3

        cr2 = tr1 + tr4
        cr4 = tr1 - tr4
        ci2 = ti1 + ti4
        ci4 = ti1 - ti4

        ch(i-1,k,2) = wa1(i-1) * cr2 - wa1(i) * ci2
        ch(i,k,2)   = wa1(i-1) * ci2 + wa1(i) * cr2
        ch(i-1,k,3) = wa2(i-1) * cr3 - wa2(i) * ci3
        ch(i,k,3)   = wa2(i-1) * ci3 + wa2(i) * cr3
        ch(i-1,k,4) = wa3(i-1) * cr4 - wa3(i) * ci4
        ch(i,k,4)   = wa3(i-1) * ci4 + wa3(i) * cr4

      end do
    end do

  end if

  return
end
subroutine passb5 ( ido, l1, cc, ch, wa1, wa2, wa3, wa4 )
!
!*******************************************************************************
!
!! PASSB5 is a lower level routine used by CFFTB1.
!
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
  integer(4) ido
  integer(4) l1
!
  real(8) cc(ido,5,l1)
  real(8) ch(ido,l1,5)
  real(8) ci2
  real(8) ci3
  real(8) ci4
  real(8) ci5
  real(8) cr2
  real(8) cr3
  real(8) cr4
  real(8) cr5
  real(8) di2
  real(8) di3
  real(8) di4
  real(8) di5
  real(8) dr2
  real(8) dr3
  real(8) dr4
  real(8) dr5
  integer(4) i
  integer(4) k
  real(8), parameter :: ti11 = 0.951056516295154D+00
  real(8), parameter :: ti12 = 0.587785252292473D+00
  real(8) ti2
  real(8) ti3
  real(8) ti4
  real(8) ti5
  real(8), parameter :: tr11 = 0.309016994374947D+00
  real(8), parameter :: tr12 = -0.809016994374947D+00
  real(8) tr2
  real(8) tr3
  real(8) tr4
  real(8) tr5
  real(8) wa1(ido)
  real(8) wa2(ido)
  real(8) wa3(ido)
  real(8) wa4(ido)
!
  if ( ido == 2 ) then

    do k = 1, l1

      ti5 = cc(2,2,k) - cc(2,5,k)
      ti2 = cc(2,2,k) + cc(2,5,k)
      ti4 = cc(2,3,k) - cc(2,4,k)
      ti3 = cc(2,3,k) + cc(2,4,k)
      tr5 = cc(1,2,k) - cc(1,5,k)
      tr2 = cc(1,2,k) + cc(1,5,k)
      tr4 = cc(1,3,k) - cc(1,4,k)
      tr3 = cc(1,3,k) + cc(1,4,k)

      ch(1,k,1) = cc(1,1,k) + tr2 + tr3
      ch(2,k,1) = cc(2,1,k) + ti2 + ti3

      cr2 = cc(1,1,k) + tr11 * tr2 + tr12 * tr3
      ci2 = cc(2,1,k) + tr11 * ti2 + tr12 * ti3
      cr3 = cc(1,1,k) + tr12 * tr2 + tr11 * tr3
      ci3 = cc(2,1,k) + tr12 * ti2 + tr11 * ti3

      cr5 = ti11 * tr5 + ti12 * tr4
      ci5 = ti11 * ti5 + ti12 * ti4
      cr4 = ti12 * tr5 - ti11 * tr4
      ci4 = ti12 * ti5 - ti11 * ti4

      ch(1,k,2) = cr2 - ci5
      ch(1,k,5) = cr2 + ci5
      ch(2,k,2) = ci2 + cr5
      ch(2,k,3) = ci3 + cr4
      ch(1,k,3) = cr3 - ci4
      ch(1,k,4) = cr3 + ci4
      ch(2,k,4) = ci3 - cr4
      ch(2,k,5) = ci2 - cr5

    end do

  else

    do k = 1, l1
      do i = 2, ido, 2

        ti5 = cc(i,2,k) - cc(i,5,k)
        ti2 = cc(i,2,k) + cc(i,5,k)
        ti4 = cc(i,3,k) - cc(i,4,k)
        ti3 = cc(i,3,k) + cc(i,4,k)

        tr5 = cc(i-1,2,k) - cc(i-1,5,k)
        tr2 = cc(i-1,2,k) + cc(i-1,5,k)
        tr4 = cc(i-1,3,k) - cc(i-1,4,k)
        tr3 = cc(i-1,3,k) + cc(i-1,4,k)

        ch(i-1,k,1) = cc(i-1,1,k) + tr2 + tr3
        ch(i,k,1)   = cc(i,1,k)   + ti2 + ti3

        cr2 = cc(i-1,1,k) + tr11 * tr2 + tr12 * tr3
        ci2 = cc(i,1,k)   + tr11 * ti2 + tr12 * ti3
        cr3 = cc(i-1,1,k) + tr12 * tr2 + tr11 * tr3
        ci3 = cc(i,1,k)   + tr12 * ti2 + tr11 * ti3

        cr5 = ti11 * tr5 + ti12 * tr4
        ci5 = ti11 * ti5 + ti12 * ti4
        cr4 = ti12 * tr5 - ti11 * tr4
        ci4 = ti12 * ti5 - ti11 * ti4

        dr3 = cr3 - ci4
        dr4 = cr3 + ci4
        di3 = ci3 + cr4
        di4 = ci3 - cr4
        dr5 = cr2 + ci5
        dr2 = cr2 - ci5
        di5 = ci2 - cr5
        di2 = ci2 + cr5

        ch(i-1,k,2) = wa1(i-1) * dr2 - wa1(i) * di2
        ch(i,k,2)   = wa1(i-1) * di2 + wa1(i) * dr2
        ch(i-1,k,3) = wa2(i-1) * dr3 - wa2(i) * di3
        ch(i,k,3)   = wa2(i-1) * di3 + wa2(i) * dr3
        ch(i-1,k,4) = wa3(i-1) * dr4 - wa3(i) * di4
        ch(i,k,4)   = wa3(i-1) * di4 + wa3(i) * dr4
        ch(i-1,k,5) = wa4(i-1) * dr5 - wa4(i) * di5
        ch(i,k,5)   = wa4(i-1) * di5 + wa4(i) * dr5

      end do
    end do

  end if

  return
end
subroutine passf ( nac, ido, ip, l1, idl1, cc, c1, c2, ch, ch2, wa )
!
!*******************************************************************************
!
!! PASSF is a lower level routine used by CFFTF1.
!
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
  integer(4) idl1
  integer(4) ido
  integer(4) ip
  integer(4) l1
!
  real(8) c1(ido,l1,ip)
  real(8) c2(idl1,ip)
  real(8) cc(ido,ip,l1)
  real(8) ch(ido,l1,ip)
  real(8) ch2(idl1,ip)
  integer(4) i
  integer(4) idij
  integer(4) idj
  integer(4) idl
  integer(4) idlj
  integer(4) idp
  integer(4) ik
  integer(4) inc
  integer(4) ipph
  integer(4) j
  integer(4) jc
  integer(4) k
  integer(4) l
  integer(4) lc
  integer(4) nac
  integer(4) nt
  real(8) wa(*)
  real(8) wai
  real(8) war
!
  nt = ip * idl1
  ipph = (ip+1) / 2
  idp = ip * ido

  if ( ido >= l1 ) then

    do j = 2, ipph
      jc = ip + 2 - j
      ch(1:ido,1:l1,j)  = cc(1:ido,j,1:l1) + cc(1:ido,jc,1:l1)
      ch(1:ido,1:l1,jc) = cc(1:ido,j,1:l1) - cc(1:ido,jc,1:l1)
    end do

    ch(1:ido,1:l1,1) = cc(1:ido,1,1:l1)

  else

    do j = 2, ipph
      jc = ip + 2 - j
      ch(1:ido,1:l1,j)  = cc(1:ido,j,1:l1) + cc(1:ido,jc,1:l1)
      ch(1:ido,1:l1,jc) = cc(1:ido,j,1:l1) - cc(1:ido,jc,1:l1)
    end do

    ch(1:ido,1:l1,1) = cc(1:ido,1,1:l1)

  end if

  idl = 2 - ido
  inc = 0

  do l = 2, ipph

    lc = ip + 2 - l
    idl = idl + ido

    do ik = 1, idl1
      c2(ik,l)  = ch2(ik,1) + wa(idl-1) * ch2(ik,2)
      c2(ik,lc) =           - wa(idl)   * ch2(ik,ip)
    end do

    idlj = idl
    inc = inc + ido

    do j = 3, ipph

      jc = ip + 2 - j

      idlj = idlj + inc
      if ( idlj > idp ) then
        idlj = idlj - idp
      end if

      war = wa(idlj-1)
      wai = wa(idlj)

      do ik = 1, idl1
        c2(ik,l)  = c2(ik,l)  + war * ch2(ik,j)
        c2(ik,lc) = c2(ik,lc) - wai * ch2(ik,jc)
      end do

    end do

  end do

  do j = 2, ipph
    ch2(1:idl1,1) = ch2(1:idl1,1) + ch2(1:idl1,j)
  end do

  do j = 2, ipph
    jc = ip + 2 - j
    do ik = 2, idl1, 2
      ch2(ik-1,j)  = c2(ik-1,j) - c2(ik,jc)
      ch2(ik-1,jc) = c2(ik-1,j) + c2(ik,jc)
      ch2(ik,j)    = c2(ik,j)   + c2(ik-1,jc)
      ch2(ik,jc)   = c2(ik,j)   - c2(ik-1,jc)
    end do
  end do

  if ( ido == 2 ) then
    nac = 1
    return
  end if

  nac = 0

  c2(1:idl1,1)    = ch2(1:idl1,1)
  c1(1,1:l1,2:ip) = ch(1,1:l1,2:ip)
  c1(2,1:l1,2:ip) = ch(2,1:l1,2:ip)

  if ( ( ido / 2 ) <= l1 ) then

    idij = 0
    do j = 2, ip
      idij = idij + 2
      do i = 4, ido, 2
        idij = idij + 2
        c1(i-1,1:l1,j) = wa(idij-1) * ch(i-1,1:l1,j) + wa(idij) * ch(i,1:l1,j)
        c1(i,1:l1,j)   = wa(idij-1) * ch(i,1:l1,j)   - wa(idij) * ch(i-1,1:l1,j)
      end do
    end do

  else

    idj = 2 - ido

    do j = 2, ip
      idj = idj + ido
      do k = 1, l1
        idij = idj
        do i = 4, ido, 2
          idij = idij + 2
          c1(i-1,k,j) = wa(idij-1) * ch(i-1,k,j) + wa(idij) * ch(i,k,j)
          c1(i,k,j)   = wa(idij-1) * ch(i,k,j)   - wa(idij) * ch(i-1,k,j)
        end do
      end do
    end do

  end if

  return
end
subroutine passf2 ( ido, l1, cc, ch, wa1 )
!
!*******************************************************************************
!
!! PASSF2 is a lower level routine used by CFFTF1.
!
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
  integer(4) ido
  integer(4) l1
!
  real(8) cc(ido,2,l1)
  real(8) ch(ido,l1,2)
  integer(4) i
  integer(4) k
  real(8) ti2
  real(8) tr2
  real(8) wa1(ido)
!
  if ( ido <= 2 ) then

    ch(1,1:l1,1) = cc(1,1,1:l1) + cc(1,2,1:l1)
    ch(1,1:l1,2) = cc(1,1,1:l1) - cc(1,2,1:l1)
    ch(2,1:l1,1) = cc(2,1,1:l1) + cc(2,2,1:l1)
    ch(2,1:l1,2) = cc(2,1,1:l1) - cc(2,2,1:l1)

  else

    do k = 1, l1
      do i = 2, ido, 2

        ch(i-1,k,1) = cc(i-1,1,k) + cc(i-1,2,k)
        tr2         = cc(i-1,1,k) - cc(i-1,2,k)

        ch(i,k,1) = cc(i,1,k) + cc(i,2,k)
        ti2       = cc(i,1,k) - cc(i,2,k)

        ch(i,k,2)   = wa1(i-1) * ti2 - wa1(i) * tr2
        ch(i-1,k,2) = wa1(i-1) * tr2 + wa1(i) * ti2

      end do
    end do

  end if

  return
end
subroutine passf3 ( ido, l1, cc, ch, wa1, wa2 )
!
!*******************************************************************************
!
!! PASSF3 is a lower level routine used by CFFTF1.
!
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
  integer(4) ido
  integer(4) l1
!
  real(8) cc(ido,3,l1)
  real(8) ch(ido,l1,3)
  real(8) ci2
  real(8) ci3
  real(8) cr2
  real(8) cr3
  real(8) di2
  real(8) di3
  real(8) dr2
  real(8) dr3
  integer(4) i
  integer(4) k
  real(8), parameter :: taui = -0.866025403784439D+00
  real(8), parameter :: taur = -0.5D+00
  real(8) ti2
  real(8) tr2
  real(8) wa1(ido)
  real(8) wa2(ido)
!
  if ( ido == 2 ) then

    do k = 1, l1

      tr2 = cc(1,2,k) + cc(1,3,k)
      cr2 = cc(1,1,k) + taur * tr2
      ch(1,k,1) = cc(1,1,k) + tr2

      ti2 = cc(2,2,k) + cc(2,3,k)
      ci2 = cc(2,1,k) + taur * ti2
      ch(2,k,1) = cc(2,1,k) + ti2

      cr3 = taui * ( cc(1,2,k) - cc(1,3,k) )
      ci3 = taui * ( cc(2,2,k) - cc(2,3,k) )

      ch(1,k,2) = cr2 - ci3
      ch(1,k,3) = cr2 + ci3
      ch(2,k,2) = ci2 + cr3
      ch(2,k,3) = ci2 - cr3

    end do

  else

    do k = 1, l1
      do i = 2, ido, 2

        tr2 = cc(i-1,2,k) + cc(i-1,3,k)
        cr2 = cc(i-1,1,k) + taur * tr2
        ch(i-1,k,1) = cc(i-1,1,k) + tr2

        ti2 = cc(i,2,k) + cc(i,3,k)
        ci2 = cc(i,1,k) + taur * ti2
        ch(i,k,1) = cc(i,1,k) + ti2

        cr3 = taui * ( cc(i-1,2,k) - cc(i-1,3,k) )
        ci3 = taui * ( cc(i,2,k)   - cc(i,3,k) )

        dr2 = cr2 - ci3
        dr3 = cr2 + ci3
        di2 = ci2 + cr3
        di3 = ci2 - cr3

        ch(i,k,2)   = wa1(i-1) * di2 - wa1(i) * dr2
        ch(i-1,k,2) = wa1(i-1) * dr2 + wa1(i) * di2
        ch(i,k,3)   = wa2(i-1) * di3 - wa2(i) * dr3
        ch(i-1,k,3) = wa2(i-1) * dr3 + wa2(i) * di3

      end do
    end do

  end if

  return
end
subroutine passf4 ( ido, l1, cc, ch, wa1, wa2, wa3 )
!
!*******************************************************************************
!
!! PASSF4 is a lower level routine used by CFFTF1.
!
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
  integer(4) ido
  integer(4) l1
!
  real(8) cc(ido,4,l1)
  real(8) ch(ido,l1,4)
  real(8) ci1
  real(8) ci2
  real(8) ci3
  real(8) ci4
  real(8) cr1
  real(8) cr2
  real(8) cr3
  real(8) cr4
  integer(4) i
  integer(4) k
  real(8) ti1
  real(8) ti2
  real(8) ti3
  real(8) ti4
  real(8) tr1
  real(8) tr2
  real(8) tr3
  real(8) tr4
  real(8) wa1(ido)
  real(8) wa2(ido)
  real(8) wa3(ido)
!
  if ( ido == 2 ) then

    do k = 1, l1

      ti1 = cc(2,1,k) - cc(2,3,k)
      ti2 = cc(2,1,k) + cc(2,3,k)
      tr4 = cc(2,2,k) - cc(2,4,k)
      ti3 = cc(2,2,k) + cc(2,4,k)
      tr1 = cc(1,1,k) - cc(1,3,k)
      tr2 = cc(1,1,k) + cc(1,3,k)
      ti4 = cc(1,4,k) - cc(1,2,k)
      tr3 = cc(1,2,k) + cc(1,4,k)

      ch(1,k,1) = tr2 + tr3
      ch(1,k,3) = tr2 - tr3
      ch(2,k,1) = ti2 + ti3
      ch(2,k,3) = ti2 - ti3
      ch(1,k,2) = tr1 + tr4
      ch(1,k,4) = tr1 - tr4
      ch(2,k,2) = ti1 + ti4
      ch(2,k,4) = ti1 - ti4

    end do

  else

    do k = 1, l1
      do i = 2, ido, 2

        ti1 = cc(i,1,k)   - cc(i,3,k)
        ti2 = cc(i,1,k)   + cc(i,3,k)
        ti3 = cc(i,2,k)   + cc(i,4,k)
        tr4 = cc(i,2,k)   - cc(i,4,k)
        tr1 = cc(i-1,1,k) - cc(i-1,3,k)
        tr2 = cc(i-1,1,k) + cc(i-1,3,k)
        ti4 = cc(i-1,4,k) - cc(i-1,2,k)
        tr3 = cc(i-1,2,k) + cc(i-1,4,k)

        ch(i-1,k,1) = tr2 + tr3
        cr3         = tr2 - tr3
        ch(i,k,1)   = ti2 + ti3
        ci3         = ti2 - ti3

        cr2 = tr1 + tr4
        cr4 = tr1 - tr4
        ci2 = ti1 + ti4
        ci4 = ti1 - ti4

        ch(i-1,k,2) = wa1(i-1) * cr2 + wa1(i) * ci2
        ch(i,k,2)   = wa1(i-1) * ci2 - wa1(i) * cr2
        ch(i-1,k,3) = wa2(i-1) * cr3 + wa2(i) * ci3
        ch(i,k,3)   = wa2(i-1) * ci3 - wa2(i) * cr3
        ch(i-1,k,4) = wa3(i-1) * cr4 + wa3(i) * ci4
        ch(i,k,4)   = wa3(i-1) * ci4 - wa3(i) * cr4

      end do
    end do

  end if

  return
end
subroutine passf5 ( ido, l1, cc, ch, wa1, wa2, wa3, wa4 )
!
!*******************************************************************************
!
!! PASSF5 is a lower level routine used by CFFTF1.
!
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
  integer(4) ido
  integer(4) l1
!
  real(8) cc(ido,5,l1)
  real(8) ch(ido,l1,5)
  real(8) ci2
  real(8) ci3
  real(8) ci4
  real(8) ci5
  real(8) cr2
  real(8) cr3
  real(8) cr4
  real(8) cr5
  real(8) di2
  real(8) di3
  real(8) di4
  real(8) di5
  real(8) dr2
  real(8) dr3
  real(8) dr4
  real(8) dr5
  integer(4) i
  integer(4) k
  real(8), parameter :: ti11 = -0.951056516295154D+00
  real(8), parameter :: ti12 = -0.587785252292473D+00
  real(8) ti2
  real(8) ti3
  real(8) ti4
  real(8) ti5
  real(8) tr2
  real(8) tr3
  real(8) tr4
  real(8) tr5
  real(8), parameter :: tr11 =  0.309016994374947D+00
  real(8), parameter :: tr12 = -0.809016994374947D+00
  real(8) wa1(ido)
  real(8) wa2(ido)
  real(8) wa3(ido)
  real(8) wa4(ido)
!
  if ( ido == 2 ) then

    do k = 1, l1

      ti5 = cc(2,2,k) - cc(2,5,k)
      ti2 = cc(2,2,k) + cc(2,5,k)
      ti4 = cc(2,3,k) - cc(2,4,k)
      ti3 = cc(2,3,k) + cc(2,4,k)
      tr5 = cc(1,2,k) - cc(1,5,k)
      tr2 = cc(1,2,k) + cc(1,5,k)
      tr4 = cc(1,3,k) - cc(1,4,k)
      tr3 = cc(1,3,k) + cc(1,4,k)

      ch(1,k,1) = cc(1,1,k) + tr2 + tr3
      ch(2,k,1) = cc(2,1,k) + ti2 + ti3

      cr2 = cc(1,1,k) + tr11 * tr2 + tr12 * tr3
      ci2 = cc(2,1,k) + tr11 * ti2 + tr12 * ti3
      cr3 = cc(1,1,k) + tr12 * tr2 + tr11 * tr3
      ci3 = cc(2,1,k) + tr12 * ti2 + tr11 * ti3

      cr5 = ti11 * tr5 + ti12 * tr4
      ci5 = ti11 * ti5 + ti12 * ti4
      cr4 = ti12 * tr5 - ti11 * tr4
      ci4 = ti12 * ti5 - ti11 * ti4

      ch(1,k,2) = cr2 - ci5
      ch(1,k,5) = cr2 + ci5
      ch(2,k,2) = ci2 + cr5
      ch(2,k,3) = ci3 + cr4
      ch(1,k,3) = cr3 - ci4
      ch(1,k,4) = cr3 + ci4
      ch(2,k,4) = ci3 - cr4
      ch(2,k,5) = ci2 - cr5

    end do

  else

    do k = 1, l1
      do i = 2, ido, 2

        ti5 = cc(i,2,k) - cc(i,5,k)
        ti2 = cc(i,2,k) + cc(i,5,k)
        ti4 = cc(i,3,k) - cc(i,4,k)
        ti3 = cc(i,3,k) + cc(i,4,k)

        tr5 = cc(i-1,2,k) - cc(i-1,5,k)
        tr2 = cc(i-1,2,k) + cc(i-1,5,k)
        tr4 = cc(i-1,3,k) - cc(i-1,4,k)
        tr3 = cc(i-1,3,k) + cc(i-1,4,k)

        ch(i-1,k,1) = cc(i-1,1,k) + tr2 + tr3
        ch(i,k,1)   = cc(i,1,k)   + ti2 + ti3

        cr2 = cc(i-1,1,k) + tr11 * tr2 + tr12 * tr3
        ci2 = cc(i,1,k)   + tr11 * ti2 + tr12 * ti3
        cr3 = cc(i-1,1,k) + tr12 * tr2 + tr11 * tr3
        ci3 = cc(i,1,k)   + tr12 * ti2 + tr11 * ti3

        cr5 = ti11 * tr5 + ti12 * tr4
        ci5 = ti11 * ti5 + ti12 * ti4
        cr4 = ti12 * tr5 - ti11 * tr4
        ci4 = ti12 * ti5 - ti11 * ti4

        dr3 = cr3 - ci4
        dr4 = cr3 + ci4
        di3 = ci3 + cr4
        di4 = ci3 - cr4
        dr5 = cr2 + ci5
        dr2 = cr2 - ci5
        di5 = ci2 - cr5
        di2 = ci2 + cr5

        ch(i-1,k,2) = wa1(i-1) * dr2 + wa1(i) * di2
        ch(i,k,2)   = wa1(i-1) * di2 - wa1(i) * dr2
        ch(i-1,k,3) = wa2(i-1) * dr3 + wa2(i) * di3
        ch(i,k,3)   = wa2(i-1) * di3 - wa2(i) * dr3
        ch(i-1,k,4) = wa3(i-1) * dr4 + wa3(i) * di4
        ch(i,k,4)   = wa3(i-1) * di4 - wa3(i) * dr4
        ch(i-1,k,5) = wa4(i-1) * dr5 + wa4(i) * di5
        ch(i,k,5)   = wa4(i-1) * di5 - wa4(i) * dr5

      end do
    end do

  end if

  return
end
function pimach ()
!
!*******************************************************************************
!
!! PIMACH returns the value of pi.
!
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
!    Output, real(8) PIMACH, the value of PI.
!
  real(8) pimach
!
  pimach = 4.0D+00 * atan ( 1.0D+00 )

  return
end
subroutine r_random ( rlo, rhi, r )
!
!*******************************************************************************
!
!! R_RANDOM returns a random real(8) in a given range.
!
!
!  Discussion:
!
!    Calls to the FORTRAN 90 random number generator should go through
!    this routine, to guarantee that the random number seed has been set.
!
!  Modified:
!
!    05 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real(8) RLO, RHI, the minimum and maximum values.
!
!    Output, real(8) R, the randomly chosen value.
!
  real(8) r
  real(8) rhi
  real(8) rlo
  integer(4), save :: seed = 0
  logical, save :: seeded = .false.
  real(8) t
  real(8) uniform_01_sample
!
!  Make sure the random number generator has been seeded.
!
  if ( .not. seeded ) then
    call random_initialize ( seed )
    seeded = .true.
  end if
!
!  Pick T, a random number in (0,1).
!
! call random_number ( harvest = t )
!
  t = uniform_01_sample ( seed )
!
!  Set R in ( RLO, RHI ).
!
  r = ( 1.0D+00 - t ) * rlo + t * rhi

  return
end
subroutine r_swap ( x, y )
!
!*******************************************************************************
!
!! R_SWAP swaps two real values.
!
!
!  Modified:
!
!    01 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real X, Y.  On output, the values of X and
!    Y have been interchanged.
!
  real(8) x
  real(8) y
  real(8) z
!
  z = x
  x = y
  y = z

  return
end
subroutine radb2 ( ido, l1, cc, ch, wa1 )
!
!*******************************************************************************
!
!! RADB2 is a lower level routine used by RFFTB1.
!
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
  integer(4) ido
  integer(4) l1
!
  real(8) cc(ido,2,l1)
  real(8) ch(ido,l1,2)
  integer(4) i
  integer(4) ic
  integer(4) k
  real(8) ti2
  real(8) tr2
  real(8) wa1(ido)
!
  ch(1,1:l1,1) = cc(1,1,1:l1) + cc(ido,2,1:l1)
  ch(1,1:l1,2) = cc(1,1,1:l1) - cc(ido,2,1:l1)

  if ( ido < 2 ) then
    return
  end if

  if ( ido > 2 ) then

    do k = 1, l1
      do i = 3, ido, 2

        ic = ido + 2 - i

        ch(i-1,k,1) = cc(i-1,1,k) + cc(ic-1,2,k)
        tr2         = cc(i-1,1,k) - cc(ic-1,2,k)
        ch(i,k,1)   = cc(i,1,k)   - cc(ic,2,k)
        ti2         = cc(i,1,k)   + cc(ic,2,k)

        ch(i-1,k,2) = wa1(i-2) * tr2 - wa1(i-1) * ti2
        ch(i,k,2)   = wa1(i-2) * ti2 + wa1(i-1) * tr2

      end do
    end do

    if ( mod ( ido, 2 ) == 1 ) then
      return
    end if

  end if

  ch(ido,1:l1,1) =    cc(ido,1,1:l1) + cc(ido,1,1:l1)
  ch(ido,1:l1,2) = -( cc(1,2,1:l1)   + cc(1,2,1:l1) )

  return
end
subroutine radb3 ( ido, l1, cc, ch, wa1, wa2 )
!
!*******************************************************************************
!
!! RADB3 is a lower level routine used by RFFTB1.
!
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
  integer(4) ido
  integer(4) l1
!
  real(8) cc(ido,3,l1)
  real(8) ch(ido,l1,3)
  real(8) ci2
  real(8) ci3
  real(8) cr2
  real(8) cr3
  real(8) di2
  real(8) di3
  real(8) dr2
  real(8) dr3
  integer(4) i
  integer(4) ic
  integer(4) k
  real(8), parameter :: taui =  0.866025403784439D+00
  real(8), parameter :: taur = -0.5D+00
  real(8) ti2
  real(8) tr2
  real(8) wa1(ido)
  real(8) wa2(ido)
!
  do k = 1, l1

    tr2 = cc(ido,2,k) + cc(ido,2,k)
    cr2 = cc(1,1,k) + taur * tr2
    ch(1,k,1) = cc(1,1,k) + tr2
    ci3 = taui * ( cc(1,3,k) + cc(1,3,k) )

    ch(1,k,2) = cr2 - ci3
    ch(1,k,3) = cr2 + ci3

  end do

  if ( ido == 1 ) then
    return
  end if

  do k = 1, l1
    do i = 3, ido, 2

      ic = ido + 2 - i

      tr2 = cc(i-1,3,k) + cc(ic-1,2,k)
      cr2 = cc(i-1,1,k) + taur * tr2
      ch(i-1,k,1) = cc(i-1,1,k) + tr2

      ti2 = cc(i,3,k) - cc(ic,2,k)
      ci2 = cc(i,1,k) + taur * ti2
      ch(i,k,1) = cc(i,1,k) + ti2

      cr3 = taui * ( cc(i-1,3,k) - cc(ic-1,2,k) )
      ci3 = taui * ( cc(i,3,k)   + cc(ic,2,k) )

      dr2 = cr2 - ci3
      dr3 = cr2 + ci3
      di2 = ci2 + cr3
      di3 = ci2 - cr3

      ch(i-1,k,2) = wa1(i-2) * dr2 - wa1(i-1) * di2
      ch(i,k,2)   = wa1(i-2) * di2 + wa1(i-1) * dr2
      ch(i-1,k,3) = wa2(i-2) * dr3 - wa2(i-1) * di3
      ch(i,k,3)   = wa2(i-2) * di3 + wa2(i-1) * dr3

    end do
  end do

  return
end
subroutine radb4 ( ido, l1, cc, ch, wa1, wa2, wa3 )
!
!*******************************************************************************
!
!! RADB4 is a lower level routine used by RFFTB1.
!
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
  integer(4) ido
  integer(4) l1
!
  real(8) cc(ido,4,l1)
  real(8) ch(ido,l1,4)
  real(8) ci2
  real(8) ci3
  real(8) ci4
  real(8) cr2
  real(8) cr3
  real(8) cr4
  integer(4) i
  integer(4) ic
  integer(4) k
  real(8), parameter :: sqrt2 = 1.414213562373095D+00
  real(8) ti1
  real(8) ti2
  real(8) ti3
  real(8) ti4
  real(8) tr1
  real(8) tr2
  real(8) tr3
  real(8) tr4
  real(8) wa1(ido)
  real(8) wa2(ido)
  real(8) wa3(ido)
!
  do k = 1, l1

    tr1 = cc(1,1,k) - cc(ido,4,k)
    tr2 = cc(1,1,k) + cc(ido,4,k)
    tr3 = cc(ido,2,k) + cc(ido,2,k)
    tr4 = cc(1,3,k) + cc(1,3,k)

    ch(1,k,1) = tr2 + tr3
    ch(1,k,2) = tr1 - tr4
    ch(1,k,3) = tr2 - tr3
    ch(1,k,4) = tr1 + tr4

  end do

  if ( ido < 2 ) then
    return
  end if

  if ( ido > 2 ) then

    do k = 1, l1
      do i = 3, ido, 2

        ic = ido + 2 - i

        ti1 = cc(i,1,k) + cc(ic,4,k)
        ti2 = cc(i,1,k) - cc(ic,4,k)
        ti3 = cc(i,3,k) - cc(ic,2,k)
        tr4 = cc(i,3,k) + cc(ic,2,k)

        tr1 = cc(i-1,1,k) - cc(ic-1,4,k)
        tr2 = cc(i-1,1,k) + cc(ic-1,4,k)
        ti4 = cc(i-1,3,k) - cc(ic-1,2,k)
        tr3 = cc(i-1,3,k) + cc(ic-1,2,k)

        ch(i-1,k,1) = tr2 + tr3
        cr3         = tr2 - tr3
        ch(i,k,1)   = ti2 + ti3
        ci3         = ti2 - ti3

        cr2 = tr1 - tr4
        cr4 = tr1 + tr4
        ci2 = ti1 + ti4
        ci4 = ti1 - ti4

        ch(i-1,k,2) = wa1(i-2) * cr2 - wa1(i-1) * ci2
        ch(i,k,2)   = wa1(i-2) * ci2 + wa1(i-1) * cr2
        ch(i-1,k,3) = wa2(i-2) * cr3 - wa2(i-1) * ci3
        ch(i,k,3)   = wa2(i-2) * ci3 + wa2(i-1) * cr3
        ch(i-1,k,4) = wa3(i-2) * cr4 - wa3(i-1) * ci4
        ch(i,k,4)   = wa3(i-2) * ci4 + wa3(i-1) * cr4

      end do
    end do

    if ( mod ( ido, 2 ) == 1 ) then
      return
    end if

  end if

  do k = 1, l1

    ti1 = cc(1,2,k)   + cc(1,4,k)
    ti2 = cc(1,4,k)   - cc(1,2,k)
    tr1 = cc(ido,1,k) - cc(ido,3,k)
    tr2 = cc(ido,1,k) + cc(ido,3,k)

    ch(ido,k,1) = tr2 + tr2
    ch(ido,k,2) = sqrt2 * ( tr1 - ti1 )
    ch(ido,k,3) = ti2 + ti2
    ch(ido,k,4) = -sqrt2 * ( tr1 + ti1 )

  end do

  return
end
subroutine radb5 ( ido, l1, cc, ch, wa1, wa2, wa3, wa4 )
!
!*******************************************************************************
!
!! RADB5 is a lower level routine used by RFFTB1.
!
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
  integer(4) ido
  integer(4) l1
!
  real(8) cc(ido,5,l1)
  real(8) ch(ido,l1,5)
  real(8) ci2
  real(8) ci3
  real(8) ci4
  real(8) ci5
  real(8) cr2
  real(8) cr3
  real(8) cr4
  real(8) cr5
  real(8) di2
  real(8) di3
  real(8) di4
  real(8) di5
  real(8) dr2
  real(8) dr3
  real(8) dr4
  real(8) dr5
  integer(4) i
  integer(4) ic
  integer(4) k
  real(8), parameter :: ti11 =  0.951056516295154D+00
  real(8), parameter :: ti12 =  0.587785252292473D+00
  real(8) ti2
  real(8) ti3
  real(8) ti4
  real(8) ti5
  real(8), parameter :: tr11 =  0.309016994374947D+00
  real(8), parameter :: tr12 = -0.809016994374947D+00
  real(8) tr2
  real(8) tr3
  real(8) tr4
  real(8) tr5
  real(8) wa1(ido)
  real(8) wa2(ido)
  real(8) wa3(ido)
  real(8) wa4(ido)
!
  do k = 1, l1

    ti5 = cc(1,3,k) + cc(1,3,k)
    ti4 = cc(1,5,k) + cc(1,5,k)
    tr2 = cc(ido,2,k) + cc(ido,2,k)
    tr3 = cc(ido,4,k) + cc(ido,4,k)

    ch(1,k,1) = cc(1,1,k) + tr2 + tr3
    cr2 = cc(1,1,k) + tr11 * tr2 + tr12 * tr3
    cr3 = cc(1,1,k) + tr12 * tr2 + tr11 * tr3
    ci5 = ti11 * ti5 + ti12 * ti4
    ci4 = ti12 * ti5 - ti11 * ti4

    ch(1,k,2) = cr2 - ci5
    ch(1,k,3) = cr3 - ci4
    ch(1,k,4) = cr3 + ci4
    ch(1,k,5) = cr2 + ci5

  end do

  if ( ido == 1 ) then
    return
  end if

  do k = 1, l1
    do i = 3, ido, 2

      ic = ido + 2 - i

      ti5 = cc(i,3,k) + cc(ic,2,k)
      ti2 = cc(i,3,k) - cc(ic,2,k)
      ti4 = cc(i,5,k) + cc(ic,4,k)
      ti3 = cc(i,5,k) - cc(ic,4,k)
      tr5 = cc(i-1,3,k) - cc(ic-1,2,k)
      tr2 = cc(i-1,3,k) + cc(ic-1,2,k)
      tr4 = cc(i-1,5,k) - cc(ic-1,4,k)
      tr3 = cc(i-1,5,k) + cc(ic-1,4,k)

      ch(i-1,k,1) = cc(i-1,1,k) + tr2 + tr3
      ch(i,k,1)   = cc(i,1,k) + ti2 + ti3

      cr2 = cc(i-1,1,k) + tr11 * tr2 + tr12 * tr3
      ci2 = cc(i,1,k)   + tr11 * ti2 + tr12 * ti3
      cr3 = cc(i-1,1,k) + tr12 * tr2 + tr11 * tr3
      ci3 = cc(i,1,k)   + tr12 * ti2 + tr11 * ti3

      cr5 = ti11 * tr5 + ti12 * tr4
      ci5 = ti11 * ti5 + ti12 * ti4
      cr4 = ti12 * tr5 - ti11 * tr4
      ci4 = ti12 * ti5 - ti11 * ti4

      dr3 = cr3 - ci4
      dr4 = cr3 + ci4
      di3 = ci3 + cr4
      di4 = ci3 - cr4
      dr5 = cr2 + ci5
      dr2 = cr2 - ci5
      di5 = ci2 - cr5
      di2 = ci2 + cr5

      ch(i-1,k,2) = wa1(i-2) * dr2 - wa1(i-1) * di2
      ch(i,k,2)   = wa1(i-2) * di2 + wa1(i-1) * dr2
      ch(i-1,k,3) = wa2(i-2) * dr3 - wa2(i-1) * di3
      ch(i,k,3)   = wa2(i-2) * di3 + wa2(i-1) * dr3
      ch(i-1,k,4) = wa3(i-2) * dr4 - wa3(i-1) * di4
      ch(i,k,4)   = wa3(i-2) * di4 + wa3(i-1) * dr4
      ch(i-1,k,5) = wa4(i-2) * dr5 - wa4(i-1) * di5
      ch(i,k,5)   = wa4(i-2) * di5 + wa4(i-1) * dr5

    end do
  end do

  return
end
subroutine radbg ( ido, ip, l1, idl1, cc, c1, c2, ch, ch2, wa )
!
!*******************************************************************************
!
!! RADBG is a lower level routine used by RFFTB1.
!
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
  integer(4) idl1
  integer(4) ido
  integer(4) ip
  integer(4) l1
!
  real(8) ai1
  real(8) ai2
  real(8) ar1
  real(8) ar1h
  real(8) ar2
  real(8) ar2h
  real(8) arg
  real(8) c1(ido,l1,ip)
  real(8) c2(idl1,ip)
  real(8) cc(ido,ip,l1)
  real(8) ch(ido,l1,ip)
  real(8) ch2(idl1,ip)
  real(8) dc2
  real(8) dcp
  real(8) ds2
  real(8) dsp
  integer(4) i
  integer(4) ic
  integer(4) idij
  integer(4) ik
  integer(4) ipph
  integer(4) is
  integer(4) j
  integer(4) j2
  integer(4) jc
  integer(4) k
  integer(4) l
  integer(4) lc
  integer(4) nbd
  real(8) pimach
  real(8) wa(*)
!
  arg = 2.0D+00 * pimach() / dble ( ip )
  dcp = cos ( arg )
  dsp = sin ( arg )
  nbd = ( ido - 1 ) / 2
  ipph = ( ip + 1 ) / 2
  ch(1:ido,1:l1,1) = cc(1:ido,1,1:l1)

  do j = 2, ipph
    jc = ip + 2 - j
    j2 = j + j
    ch(1,1:l1,j) =  cc(ido,j2-2,1:l1) + cc(ido,j2-2,1:l1)
    ch(1,1:l1,jc) = cc(1,j2-1,1:l1)   + cc(1,j2-1,1:l1)
  end do

  if ( ido /= 1 ) then

    if ( nbd >= l1 ) then

      do j = 2, ipph
        jc = ip + 2 - j
        do k = 1, l1
          do i = 3, ido, 2
            ic = ido + 2 - i
            ch(i-1,k,j)  = cc(i-1,2*j-1,k) + cc(ic-1,2*j-2,k)
            ch(i-1,k,jc) = cc(i-1,2*j-1,k) - cc(ic-1,2*j-2,k)
            ch(i,k,j)    = cc(i,2*j-1,k)   - cc(ic,2*j-2,k)
            ch(i,k,jc)   = cc(i,2*j-1,k)   + cc(ic,2*j-2,k)
          end do
        end do
      end do

    else

      do j = 2, ipph
        jc = ip + 2 - j
        do i = 3, ido, 2
          ic = ido + 2 - i
          ch(i-1,1:l1,j)  = cc(i-1,2*j-1,1:l1) + cc(ic-1,2*j-2,1:l1)
          ch(i-1,1:l1,jc) = cc(i-1,2*j-1,1:l1) - cc(ic-1,2*j-2,1:l1)
          ch(i,1:l1,j)    = cc(i,2*j-1,1:l1)   - cc(ic,2*j-2,1:l1)
          ch(i,1:l1,jc)   = cc(i,2*j-1,1:l1)   + cc(ic,2*j-2,1:l1)
        end do
      end do

    end if

  end if

  ar1 = 1.0D+00
  ai1 = 0.0D+00

  do l = 2, ipph

    lc = ip + 2 - l
    ar1h = dcp * ar1 - dsp * ai1
    ai1 =  dcp * ai1 + dsp * ar1
    ar1 = ar1h

    do ik = 1, idl1
      c2(ik,l)  = ch2(ik,1) + ar1 * ch2(ik,2)
      c2(ik,lc) =             ai1 * ch2(ik,ip)
    end do

    dc2 = ar1
    ds2 = ai1
    ar2 = ar1
    ai2 = ai1

    do j = 3, ipph

      jc = ip + 2 - j
      ar2h = dc2 * ar2 - ds2 * ai2
      ai2  = dc2 * ai2 + ds2 * ar2
      ar2 = ar2h

      do ik = 1, idl1
        c2(ik,l)  = c2(ik,l)  + ar2 * ch2(ik,j)
        c2(ik,lc) = c2(ik,lc) + ai2 * ch2(ik,jc)
      end do

    end do

  end do

  do j = 2, ipph
    ch2(1:idl1,1) = ch2(1:idl1,1) + ch2(1:idl1,j)
  end do

  do j = 2, ipph
    jc = ip + 2 - j
    ch(1,1:l1,j)  = c1(1,1:l1,j) - c1(1,1:l1,jc)
    ch(1,1:l1,jc) = c1(1,1:l1,j) + c1(1,1:l1,jc)
  end do

  if ( ido /= 1 ) then

    if ( nbd >= l1 ) then

      do j = 2, ipph
        jc = ip + 2 - j
        do k = 1, l1
          do i = 3, ido, 2
            ch(i-1,k,j)  = c1(i-1,k,j) - c1(i,k,jc)
            ch(i-1,k,jc) = c1(i-1,k,j) + c1(i,k,jc)
            ch(i,k,j)    = c1(i,k,j)   + c1(i-1,k,jc)
            ch(i,k,jc)   = c1(i,k,j)   - c1(i-1,k,jc)
          end do
        end do
      end do

    else

      do j = 2, ipph
        jc = ip + 2 - j
        do i = 3, ido, 2
          ch(i-1,1:l1,j)  = c1(i-1,1:l1,j) - c1(i,1:l1,jc)
          ch(i-1,1:l1,jc) = c1(i-1,1:l1,j) + c1(i,1:l1,jc)
          ch(i,1:l1,j)    = c1(i,1:l1,j)   + c1(i-1,1:l1,jc)
          ch(i,1:l1,jc)   = c1(i,1:l1,j)   - c1(i-1,1:l1,jc)
        end do
      end do

    end if

  end if

  if ( ido == 1 ) then
    return
  end if

  c2(1:idl1,1) = ch2(1:idl1,1)
  c1(1,1:l1,2:ip) = ch(1,1:l1,2:ip)

  if ( nbd <= l1 ) then

    is = -ido

    do j = 2, ip
      is = is + ido
      idij = is
      do i = 3, ido, 2
        idij = idij + 2
        c1(i-1,1:l1,j) = wa(idij-1) * ch(i-1,1:l1,j) - wa(idij) * ch(i,1:l1,j)
        c1(i,1:l1,j)   = wa(idij-1) * ch(i,1:l1,j)   + wa(idij) * ch(i-1,1:l1,j)
      end do
    end do

  else

    is = -ido
    do j = 2, ip
      is = is + ido
      do k = 1, l1
        idij = is
        do i = 3, ido, 2
          idij = idij + 2
          c1(i-1,k,j) = wa(idij-1) * ch(i-1,k,j) - wa(idij) * ch(i,k,j)
          c1(i,k,j)   = wa(idij-1) * ch(i,k,j)   + wa(idij) * ch(i-1,k,j)
        end do
      end do
    end do

  end if

  return
end
subroutine radf2 ( ido, l1, cc, ch, wa1 )
!
!*******************************************************************************
!
!! RADF2 is a lower level routine used by RFFTF1.
!
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
  integer(4) ido
  integer(4) l1
!
  real(8) cc(ido,l1,2)
  real(8) ch(ido,2,l1)
  integer(4) i
  integer(4) ic
  integer(4) k
  real(8) ti2
  real(8) tr2
  real(8) wa1(ido)
!
  ch(1,1,1:l1)   = cc(1,1:l1,1) + cc(1,1:l1,2)
  ch(ido,2,1:l1) = cc(1,1:l1,1) - cc(1,1:l1,2)

  if ( ido < 2 ) then
    return
  end if

  if ( ido > 2 ) then

    do k = 1, l1
      do i = 3, ido, 2

        ic = ido + 2 - i

        tr2 = wa1(i-2) * cc(i-1,k,2) + wa1(i-1) * cc(i,k,2)
        ti2 = wa1(i-2) * cc(i,k,2)   - wa1(i-1) * cc(i-1,k,2)

        ch(i,1,k) = cc(i,k,1) + ti2
        ch(ic,2,k) = ti2 - cc(i,k,1)
        ch(i-1,1,k) = cc(i-1,k,1) + tr2
        ch(ic-1,2,k) = cc(i-1,k,1) - tr2

      end do
    end do

    if ( mod ( ido, 2 ) == 1 ) then
      return
    end if

  end if

  ch(1,2,1:l1) = -cc(ido,1:l1,2)
  ch(ido,1,1:l1) = cc(ido,1:l1,1)

  return
end
subroutine radf3 ( ido, l1, cc, ch, wa1, wa2 )
!
!*******************************************************************************
!
!! RADF3 is a lower level routine used by RFFTF1.
!
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
  integer(4) ido
  integer(4) l1
!
  real(8) cc(ido,l1,3)
  real(8) ch(ido,3,l1)
  real(8) ci2
  real(8) cr2
  real(8) di2
  real(8) di3
  real(8) dr2
  real(8) dr3
  integer(4) i
  integer(4) ic
  integer(4) k
  real(8), parameter :: taui = 0.866025403784439D+00
  real(8), parameter :: taur = -0.5D+00
  real(8) ti2
  real(8) ti3
  real(8) tr2
  real(8) tr3
  real(8) wa1(ido)
  real(8) wa2(ido)
!
  do k = 1, l1
    cr2 = cc(1,k,2) + cc(1,k,3)
    ch(1,1,k) = cc(1,k,1) + cr2
    ch(1,3,k) = taui * ( cc(1,k,3) - cc(1,k,2) )
    ch(ido,2,k) = cc(1,k,1) + taur * cr2
  end do

  if ( ido == 1 ) then
    return
  end if

  do k = 1, l1
    do i = 3, ido, 2

      ic = ido + 2 - i

      dr2 = wa1(i-2) * cc(i-1,k,2) + wa1(i-1) * cc(i,k,2)
      di2 = wa1(i-2) * cc(i,k,2)   - wa1(i-1) * cc(i-1,k,2)
      dr3 = wa2(i-2) * cc(i-1,k,3) + wa2(i-1) * cc(i,k,3)
      di3 = wa2(i-2) * cc(i,k,3)   - wa2(i-1) * cc(i-1,k,3)

      cr2 = dr2 + dr3
      ci2 = di2 + di3

      ch(i-1,1,k) = cc(i-1,k,1) + cr2
      ch(i,1,k)   = cc(i,k,1) + ci2

      tr2 = cc(i-1,k,1) + taur * cr2
      ti2 = cc(i,k,1) + taur * ci2
      tr3 = taui * ( di2 - di3 )
      ti3 = taui * ( dr3 - dr2 )

      ch(i-1,3,k) = tr2 + tr3
      ch(ic-1,2,k) = tr2 - tr3
      ch(i,3,k) = ti2 + ti3
      ch(ic,2,k) = ti3 - ti2

    end do
  end do

  return
end
subroutine radf4 ( ido, l1, cc, ch, wa1, wa2, wa3 )
!
!*******************************************************************************
!
!! RADF4 is a lower level routine used by RFFTF1.
!
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
  integer(4) ido
  integer(4) l1
!
  real(8) cc(ido,l1,4)
  real(8) ch(ido,4,l1)
  real(8) ci2
  real(8) ci3
  real(8) ci4
  real(8) cr2
  real(8) cr3
  real(8) cr4
  real(8), parameter :: hsqt2 = 0.7071067811865475D+00
  integer(4) i
  integer(4) ic
  integer(4) k
  real(8) ti1
  real(8) ti2
  real(8) ti3
  real(8) ti4
  real(8) tr1
  real(8) tr2
  real(8) tr3
  real(8) tr4
  real(8) wa1(ido)
  real(8) wa2(ido)
  real(8) wa3(ido)
!
  do k = 1, l1
    tr1 = cc(1,k,2) + cc(1,k,4)
    tr2 = cc(1,k,1) + cc(1,k,3)
    ch(1,1,k) = tr1 + tr2
    ch(ido,4,k) = tr2 - tr1
    ch(ido,2,k) = cc(1,k,1) - cc(1,k,3)
    ch(1,3,k) = cc(1,k,4) - cc(1,k,2)
  end do

  if ( ido < 2 ) then
    return
  end if

  if ( ido > 2 ) then

    do k = 1, l1
      do i = 3, ido, 2

        ic = ido + 2 - i

        cr2 = wa1(i-2) * cc(i-1,k,2) + wa1(i-1) * cc(i,k,2)
        ci2 = wa1(i-2) * cc(i,k,2)   - wa1(i-1) * cc(i-1,k,2)
        cr3 = wa2(i-2) * cc(i-1,k,3) + wa2(i-1) * cc(i,k,3)
        ci3 = wa2(i-2) * cc(i,k,3)   - wa2(i-1) * cc(i-1,k,3)
        cr4 = wa3(i-2) * cc(i-1,k,4) + wa3(i-1) * cc(i,k,4)
        ci4 = wa3(i-2) * cc(i,k,4)   - wa3(i-1) * cc(i-1,k,4)

        tr1 = cr2+cr4
        tr4 = cr4-cr2
        ti1 = ci2+ci4
        ti4 = ci2-ci4
        ti2 = cc(i,k,1) + ci3
        ti3 = cc(i,k,1) - ci3
        tr2 = cc(i-1,k,1) + cr3
        tr3 = cc(i-1,k,1) - cr3

        ch(i-1,1,k)  = tr1 + tr2
        ch(ic-1,4,k) = tr2 - tr1
        ch(i,1,k)    = ti1 + ti2
        ch(ic,4,k)   = ti1 - ti2
        ch(i-1,3,k)  = ti4 + tr3
        ch(ic-1,2,k) = tr3 - ti4
        ch(i,3,k)    = tr4 + ti3
        ch(ic,2,k)   = tr4 - ti3

      end do
    end do

    if ( mod ( ido, 2 ) == 1 ) then
      return
    end if

  end if

  do k = 1, l1

    ti1 = -hsqt2 * ( cc(ido,k,2) + cc(ido,k,4) )
    tr1 =  hsqt2 * ( cc(ido,k,2) - cc(ido,k,4) )

    ch(ido,1,k) = tr1 + cc(ido,k,1)
    ch(ido,3,k) = cc(ido,k,1) - tr1

    ch(1,2,k) = ti1 - cc(ido,k,3)
    ch(1,4,k) = ti1 + cc(ido,k,3)

  end do

  return
end
subroutine radf5 ( ido, l1, cc, ch, wa1, wa2, wa3, wa4 )
!
!*******************************************************************************
!
!! RADF5 is a lower level routine used by RFFTF1.
!
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
  integer(4) ido
  integer(4) l1
!
  real(8) cc(ido,l1,5)
  real(8) ch(ido,5,l1)
  real(8) ci2
  real(8) ci3
  real(8) ci4
  real(8) ci5
  real(8) cr2
  real(8) cr3
  real(8) cr4
  real(8) cr5
  real(8) di2
  real(8) di3
  real(8) di4
  real(8) di5
  real(8) dr2
  real(8) dr3
  real(8) dr4
  real(8) dr5
  integer(4) i
  integer(4) ic
  integer(4) k
  real(8), parameter :: ti11 =  0.951056516295154D+00
  real(8), parameter :: ti12 =  0.587785252292473D+00
  real(8) ti2
  real(8) ti3
  real(8) ti4
  real(8) ti5
  real(8), parameter :: tr11 =  0.309016994374947D+00
  real(8), parameter :: tr12 = -0.809016994374947D+00
  real(8) tr2
  real(8) tr3
  real(8) tr4
  real(8) tr5
  real(8) wa1(ido)
  real(8) wa2(ido)
  real(8) wa3(ido)
  real(8) wa4(ido)
!
  do k = 1, l1

    cr2 = cc(1,k,5) + cc(1,k,2)
    ci5 = cc(1,k,5) - cc(1,k,2)
    cr3 = cc(1,k,4) + cc(1,k,3)
    ci4 = cc(1,k,4) - cc(1,k,3)

    ch(1,1,k)   = cc(1,k,1) + cr2 + cr3
    ch(ido,2,k) = cc(1,k,1) + tr11 * cr2 + tr12 * cr3
    ch(1,3,k)   = ti11 * ci5 + ti12 * ci4
    ch(ido,4,k) = cc(1,k,1) + tr12 * cr2 + tr11 * cr3
    ch(1,5,k)   = ti12 * ci5 - ti11 * ci4

  end do

  if ( ido == 1 ) then
    return
  end if

  do k = 1, l1
    do i = 3, ido, 2

      ic = ido + 2 - i

      dr2 = wa1(i-2) * cc(i-1,k,2) + wa1(i-1) * cc(i,k,2)
      di2 = wa1(i-2) * cc(i,k,2)   - wa1(i-1) * cc(i-1,k,2)
      dr3 = wa2(i-2) * cc(i-1,k,3) + wa2(i-1) * cc(i,k,3)
      di3 = wa2(i-2) * cc(i,k,3)   - wa2(i-1) * cc(i-1,k,3)
      dr4 = wa3(i-2) * cc(i-1,k,4) + wa3(i-1) * cc(i,k,4)
      di4 = wa3(i-2) * cc(i,k,4)   - wa3(i-1) * cc(i-1,k,4)
      dr5 = wa4(i-2) * cc(i-1,k,5) + wa4(i-1) * cc(i,k,5)
      di5 = wa4(i-2) * cc(i,k,5)   - wa4(i-1) * cc(i-1,k,5)

      cr2 = dr2 + dr5
      ci5 = dr5 - dr2
      cr5 = di2 - di5
      ci2 = di2 + di5
      cr3 = dr3 + dr4
      ci4 = dr4 - dr3
      cr4 = di3 - di4
      ci3 = di3 + di4

      ch(i-1,1,k) = cc(i-1,k,1) + cr2 + cr3
      ch(i,1,k)   = cc(i,k,1)   + ci2 + ci3

      tr2 = cc(i-1,k,1) + tr11 * cr2 + tr12 * cr3
      ti2 = cc(i,k,1)   + tr11 * ci2 + tr12 * ci3
      tr3 = cc(i-1,k,1) + tr12 * cr2 + tr11 * cr3
      ti3 = cc(i,k,1)   + tr12 * ci2 + tr11 * ci3

      tr5 = ti11 * cr5 + ti12 * cr4
      ti5 = ti11 * ci5 + ti12 * ci4
      tr4 = ti12 * cr5 - ti11 * cr4
      ti4 = ti12 * ci5 - ti11 * ci4

      ch(i-1,3,k)  = tr2 + tr5
      ch(ic-1,2,k) = tr2 - tr5
      ch(i,3,k)    = ti2 + ti5
      ch(ic,2,k)   = ti5 - ti2
      ch(i-1,5,k)  = tr3 + tr4
      ch(ic-1,4,k) = tr3 - tr4
      ch(i,5,k)    = ti3 + ti4
      ch(ic,4,k)   = ti4 - ti3

    end do
  end do

  return
end
subroutine radfg ( ido, ip, l1, idl1, cc, c1, c2, ch, ch2, wa )
!
!*******************************************************************************
!
!! RADFG is a lower level routine used by RFFTF1.
!
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
  integer(4) idl1
  integer(4) ido
  integer(4) ip
  integer(4) l1
!
  real(8) ai1
  real(8) ai2
  real(8) ar1
  real(8) ar1h
  real(8) ar2
  real(8) ar2h
  real(8) arg
  real(8) c1(ido,l1,ip)
  real(8) c2(idl1,ip)
  real(8) cc(ido,ip,l1)
  real(8) ch(ido,l1,ip)
  real(8) ch2(idl1,ip)
  real(8) dc2
  real(8) dcp
  real(8) ds2
  real(8) dsp
  integer(4) i
  integer(4) ic
  integer(4) idij
  integer(4) ik
  integer(4) ipph
  integer(4) is
  integer(4) j
  integer(4) j2
  integer(4) jc
  integer(4) k
  integer(4) l
  integer(4) lc
  integer(4) nbd
  real(8) pimach
  real(8) wa(*)
!
  arg = 2.0D+00 * pimach() / dble ( ip )
  dcp = cos ( arg )
  dsp = sin ( arg )
  ipph = ( ip + 1 ) / 2
  nbd = ( ido - 1 ) / 2

  if ( ido == 1 ) then

    c2(1:idl1,1) = ch2(1:idl1,1)

  else

    ch2(1:idl1,1) = c2(1:idl1,1)
    ch(1,1:l1,2:ip) = c1(1,1:l1,2:ip)

    if ( nbd <= l1 ) then

      is = -ido
      do j = 2, ip
        is = is + ido
        idij = is
        do i = 3, ido, 2
          idij = idij + 2
          do k = 1, l1
            ch(i-1,k,j) = wa(idij-1) * c1(i-1,k,j) + wa(idij) * c1(i,k,j)
            ch(i,k,j)   = wa(idij-1) * c1(i,k,j)   - wa(idij) * c1(i-1,k,j)
          end do
        end do
      end do

    else

      is = -ido

      do j = 2, ip
        is = is + ido
        do k = 1, l1
          idij = is
          do i = 3, ido, 2
            idij = idij + 2
            ch(i-1,k,j) = wa(idij-1) * c1(i-1,k,j) + wa(idij) * c1(i,k,j)
            ch(i,k,j)   = wa(idij-1) * c1(i,k,j)   - wa(idij) * c1(i-1,k,j)
          end do
        end do
      end do

    end if

    if ( nbd >= l1 ) then

      do j = 2, ipph
        jc = ip + 2 - j
        do k = 1, l1
          do i = 3, ido, 2
            c1(i-1,k,j)  = ch(i-1,k,j)  + ch(i-1,k,jc)
            c1(i-1,k,jc) = ch(i,k,j)    - ch(i,k,jc)
            c1(i,k,j)    = ch(i,k,j)    + ch(i,k,jc)
            c1(i,k,jc)   = ch(i-1,k,jc) - ch(i-1,k,j)
          end do
        end do
      end do

    else

      do j = 2, ipph
        jc = ip + 2 - j
        do i = 3, ido, 2
          c1(i-1,1:l1,j)  = ch(i-1,1:l1,j)  + ch(i-1,1:l1,jc)
          c1(i-1,1:l1,jc) = ch(i,1:l1,j)    - ch(i,1:l1,jc)
          c1(i,1:l1,j)    = ch(i,1:l1,j)    + ch(i,1:l1,jc)
          c1(i,1:l1,jc)   = ch(i-1,1:l1,jc) - ch(i-1,1:l1,j)
        end do
      end do

    end if

  end if

  do j = 2, ipph
    jc = ip + 2 - j
    c1(1,1:l1,j)  = ch(1,1:l1,j)  + ch(1,1:l1,jc)
    c1(1,1:l1,jc) = ch(1,1:l1,jc) - ch(1,1:l1,j)
  end do

  ar1 = 1.0D+00
  ai1 = 0.0D+00

  do l = 2, ipph

    lc = ip + 2 - l
    ar1h = dcp * ar1 - dsp * ai1
    ai1 =  dcp * ai1 + dsp * ar1
    ar1 = ar1h

    do ik = 1, idl1
      ch2(ik,l) = c2(ik,1) + ar1 * c2(ik,2)
      ch2(ik,lc) =           ai1 * c2(ik,ip)
    end do

    dc2 = ar1
    ds2 = ai1
    ar2 = ar1
    ai2 = ai1

    do j = 3, ipph

      jc = ip + 2 - j
      ar2h = dc2 * ar2 - ds2 * ai2
      ai2 =  dc2 * ai2 + ds2 * ar2
      ar2 = ar2h

      do ik = 1, idl1
        ch2(ik,l) =  ch2(ik,l)  + ar2 * c2(ik,j)
        ch2(ik,lc) = ch2(ik,lc) + ai2 * c2(ik,jc)
      end do

    end do

  end do

  do j = 2, ipph
    ch2(1:idl1,1) = ch2(1:idl1,1) + c2(1:idl1,j)
  end do

  cc(1:ido,1,1:l1) = ch(1:ido,1:l1,1)

  do j = 2, ipph
    jc = ip + 2 - j
    j2 = j + j
    cc(ido,j2-2,1:l1) = ch(1,1:l1,j)
    cc(1,j2-1,1:l1)   = ch(1,1:l1,jc)
  end do

  if ( ido == 1 ) then
    return
  end if

  if ( nbd >= l1 ) then

    do j = 2, ipph
      jc = ip + 2 - j
      j2 = j + j
      do k = 1, l1
        do i = 3, ido, 2
          ic = ido + 2 - i
          cc(i-1,j2-1,k)  = ch(i-1,k,j) + ch(i-1,k,jc)
          cc(ic-1,j2-2,k) = ch(i-1,k,j) - ch(i-1,k,jc)
          cc(i,j2-1,k)    = ch(i,k,j)   + ch(i,k,jc)
          cc(ic,j2-2,k)   = ch(i,k,jc)  - ch(i,k,j)
        end do
      end do
    end do

  else

    do j = 2, ipph
      jc = ip + 2 - j
      j2 = j + j
      do i = 3, ido, 2
        ic = ido + 2 - i
        cc(i-1,j2-1,1:l1)  = ch(i-1,1:l1,j) + ch(i-1,1:l1,jc)
        cc(ic-1,j2-2,1:l1) = ch(i-1,1:l1,j) - ch(i-1,1:l1,jc)
        cc(i,j2-1,1:l1)    = ch(i,1:l1,j)   + ch(i,1:l1,jc)
        cc(ic,j2-2,1:l1)   = ch(i,1:l1,jc)  - ch(i,1:l1,j)
      end do
    end do

  end if

  return
end
subroutine random_initialize ( seed )
!
!*******************************************************************************
!
!! RANDOM_INITIALIZE initializes the FORTRAN 90 random number seed.
!
!
!  Discussion:
!
!    If you don't initialize the random number generator, its behavior
!    is not specified.  If you initialize it simply by:
!
!      call random_seed
!
!    its behavior is not specified.  On the DEC ALPHA, if that's all you
!    do, the same random number sequence is returned.  In order to actually
!    try to scramble up the random number generator a bit, this routine
!    goes through the tedious process of getting the size of the random
!    number seed, making up values based on the current time, and setting
!    the random number seed.
!
!    And this is the FORTRAN 90 people's idea of convenience?
!
!    And I still get poorly randomized values, somehow, having to do
!    with a bad seed, or something.  I am about ready to go back to
!    using my own damn routine!
!
!  Modified:
!
!    06 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer(4) SEED, a seed value.
!
  integer(4) date_time(8)
  integer(4) i
  integer(4) seed
  integer(4), allocatable :: seed_vector(:)
  integer(4) seed_size
  real(8) t
  integer(4) value
!
!  Initialize the random number seed.
!
  call random_seed
!
!  Determine the size of the random number seed.
!
  call random_seed ( size = seed_size )
!
!  Allocate a seed of the right size.
!
  allocate ( seed_vector(seed_size) )
!
!  Get the current date and time.
!
  call date_and_time ( values = date_time )
!
!  Construct a slightly random value.
!
  seed = 0
  do i = 1, 8
    seed = ieor ( seed, date_time(i) )
  end do
!
!  Make slightly random assignments to SEED_VECTOR.
!
  do i = 1, seed_size
    seed_vector(i) = ieor ( seed, i )
  end do
!
!  Set the random number seed value.
!
  call random_seed ( put = seed_vector(1:seed_size) )
!
!  Free up the seed space.
!
  deallocate ( seed_vector )
!
!  Because EVEN THIS DOESN'T SEEM TO PROPERLY MIX UP THE RANDOM
!  NUMBERS, call the random number routine a bunch of times.
!
  do i = 1, 100
    call random_number ( harvest = t )
  end do
!
!  I STILL GET LOUSY RESULTS.  THE HELL WITH IT!
!
  return
end
subroutine rfftb ( n, r, wsave )
!
!*******************************************************************************
!
!! RFFTB computes a real(8) periodic sequence from its Fourier coefficients.
!
!
!  Discussion:
!
!    This process is sometimes called Fourier synthesis.
!
!    The transform is unnormalized.  A call to RFFTF followed by a call to
!    RFFTB will multiply the input sequence by N.
!
!    If N is even, the transform is defined by:
!
!      R_out(I) = R_in(1) + (-1)**(I-1) * R_in(N) + sum ( 2 <= K <= N/2 )
!
!        + 2 * R_in(2*K-2) * cos ( ( K - 1 ) * ( I - 1 ) * 2 * PI / N )
!
!        - 2 * R_in(2*K-1) * sin ( ( K - 1 ) * ( I - 1 ) * 2 * PI / N )
!
!    If N is odd, the transform is defined by:
!
!      R_out(I) = R_in(1) + sum ( 2 <= K <= (N+1)/2 )
!
!        + 2 * R_in(2*K-2) * cos ( ( K - 1 ) * ( I - 1 ) * 2 * PI / N )
!
!        - 2 * R_in(2*K-1) * sin ( ( K - 1 ) * ( I - 1 ) * 2 * PI / N )
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!    P N Swarztrauber, 
!    Vectorizing the FFT's, 
!    in Parallel Computations (G. Rodrigue, editor), 
!    Academic Press, 1982, pages 51-83.
!
!    B L Buzbee, 
!    The SLATEC Common Math Library, 
!    in Sources and Development of Mathematical Software (W. Cowell, editor),
!    Prentice Hall, 1984, pages 302-318.
!
!  Parameters:
!
!    Input, integer(4) N, the length of the array to be transformed.  The 
!    method is more efficient when N is the product of small primes.
!
!    Input/output, real(8) R(N).
!    On input, the sequence to be transformed.
!    On output, the transformed sequence.
!
!    Input, real(8) WSAVE(2*N+15), a work array.  The WSAVE array must be
!    initialized by calling RFFTI.  A different WSAVE array must be used 
!    for each different value of N.
!
  integer(4) n
!
  real(8) r(n)
  real(8) wsave(2*n+15)
!
  if ( n <= 1 ) then
    return
  end if

  call rfftb1 ( n, r, wsave(1), wsave(n+1), wsave(2*n+1) )

  return
end
subroutine rfftb1 ( n, c, ch, wa, ifac )
!
!*******************************************************************************
!
!! RFFTB1 is a lower level routine used by RFFTB.
!
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
!    Input, integer(4) N, the length of the array to be transformed.  
!
!    Input/output, real(8) C(N).
!    On input, the sequence to be transformed.
!    On output, the transformed sequence.
!
!    Input, real(8) CH(N).
!
!    Input, real(8) WA(N).
!
!    Input, integer(4) IFAC(15).
!    IFAC(1) = N, the number that was factored.
!    IFAC(2) = NF, the number of factors.
!    IFAC(3:2+NF), the factors.
!
  integer(4) n
!
  real(8) c(n)
  real(8) ch(n)
  integer(4) idl1
  integer(4) ido
  integer(4) ifac(15)
  integer(4) ip
  integer(4) iw
  integer(4) ix2
  integer(4) ix3
  integer(4) ix4
  integer(4) k1
  integer(4) l1
  integer(4) l2
  integer(4) na
  integer(4) nf
  real(8) wa(n)
!
  nf = ifac(2)
  na = 0
  l1 = 1
  iw = 1

  do k1 = 1, nf

    ip = ifac(k1+2)
    l2 = ip * l1
    ido = n / l2
    idl1 = ido * l1

    if ( ip == 4 ) then

      ix2 = iw + ido
      ix3 = ix2 + ido

      if ( na == 0 ) then
        call radb4 ( ido, l1, c, ch, wa(iw), wa(ix2), wa(ix3) )
      else
        call radb4 ( ido, l1, ch, c, wa(iw), wa(ix2), wa(ix3) )
      end if

      na = 1 - na

    else if ( ip == 2 ) then

      if ( na == 0 ) then
        call radb2 ( ido, l1, c, ch, wa(iw) )
      else
        call radb2 ( ido, l1, ch, c, wa(iw) )
      end if

      na = 1 - na

    else if ( ip == 3 ) then

      ix2 = iw + ido

      if ( na == 0 ) then
        call radb3 ( ido, l1, c, ch, wa(iw), wa(ix2) )
      else
        call radb3 ( ido, l1, ch, c, wa(iw), wa(ix2) )
      end if

      na = 1 - na

    else if ( ip == 5 ) then

      ix2 = iw + ido
      ix3 = ix2 + ido
      ix4 = ix3 + ido

      if ( na == 0 ) then
        call radb5 ( ido, l1, c, ch, wa(iw), wa(ix2), wa(ix3), wa(ix4) )
      else
        call radb5 ( ido, l1, ch, c, wa(iw), wa(ix2), wa(ix3), wa(ix4) )
      end if

      na = 1 - na

    else

      if ( na == 0 ) then
        call radbg ( ido, ip, l1, idl1, c, c, c, ch, ch, wa(iw) )
      else
        call radbg ( ido, ip, l1, idl1, ch, ch, ch, c, c, wa(iw) )
      end if

      if ( ido == 1 ) then
        na = 1 - na
      end if

    end if

    l1 = l2
    iw = iw + ( ip - 1 ) * ido

  end do

  if ( na /= 0 ) then
    c(1:n) = ch(1:n)
  end if

  return
end
subroutine rfftf ( n, r, wsave )
!
!*******************************************************************************
!
!! RFFTF computes the Fourier coefficients of a real(8) periodic sequence.
!
!
!  Discussion:
!
!    This process is sometimes called Fourier analysis.
! 
!    The transform is unnormalized.  A call to RFFTF followed by a call 
!    to RFFTB will multiply the input sequence by N.
!
!    The transform is defined by:
!
!      R_out(1) = sum ( 1 <= I <= N ) R_in(I)
!
!    Letting L = (N+1)/2, then for K = 2,...,L
!
!      R_out(2*K-2) = sum ( 1 <= I <= N )
!
!        R_in(I) * cos ( ( K - 1 ) * ( I - 1 ) * 2 * PI / N )
!
!      R_out(2*K-1) = sum ( 1 <= I <= N )
!
!        -R_in(I) * sin ( ( K - 1 ) * ( I - 1 ) * 2 * PI / N )
!
!    And, if N is even, then:
!
!      R_out(N) = sum ( 1 <= I <= N ) (-1)**(I-1) * R_in(I)
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!    P N Swarztrauber, 
!    Vectorizing the FFT's, 
!    in Parallel Computations (G. Rodrigue, editor), 
!    Academic Press, 1982, pages 51-83.
!
!    B L Buzbee, 
!    The SLATEC Common Math Library, 
!    in Sources and Development of Mathematical Software (W. Cowell, editor),
!    Prentice Hall, 1984, pages 302-318.
!
!  Parameters:
!
!    Input, integer(4) N, the length of the array to be transformed.  The 
!    method is more efficient when N is the product of small primes.
!
!    Input/output, real(8) R(N).
!    On input, the sequence to be transformed.
!    On output, the transformed sequence.
!
!    Input, real(8) WSAVE(2*N+15), a work array.  The WSAVE array must be
!    initialized by calling RFFTI.  A different WSAVE array must be used 
!    for each different value of N.  
!
  integer(4) n
!
  real(8) r(n)
  real(8) wsave(2*n+15)
!
  if ( n <= 1 ) then
    return
  end if

  call rfftf1 ( n, r, wsave(1), wsave(n+1), wsave(2*n+1) )

  return
end
subroutine rfftf1 ( n, c, ch, wa, ifac )
!
!*******************************************************************************
!
!! RFFTF1 is a lower level routine used by RFFTF and SINT.
!
!
!  Modified:
!
!    12 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
!    Input, integer(4) N, the length of the array to be transformed.  
!
!    Input/output, real(8) C(N).
!    On input, the sequence to be transformed.
!    On output, the transformed sequence.
!
!    Input, real(8) CH(N).
!
!    Input, real(8) WA(N).
!
!    Input, integer(4) IFAC(15).
!    IFAC(1) = N, the number that was factored.
!    IFAC(2) = NF, the number of factors.
!    IFAC(3:2+NF), the factors.
!
  integer(4) n
!
  real(8) c(n)
  real(8) ch(n)
  integer(4) idl1
  integer(4) ido
  integer(4) ifac(15)
  integer(4) ip
  integer(4) iw
  integer(4) ix2
  integer(4) ix3
  integer(4) ix4
  integer(4) k1
  integer(4) kh
  integer(4) l1
  integer(4) l2
  integer(4) na
  integer(4) nf
  real(8) wa(n)
!
  nf = ifac(2)
  na = 1
  l2 = n
  iw = n

  do k1 = 1, nf

    kh = nf - k1
    ip = ifac(kh+3)
    l1 = l2 / ip
    ido = n / l2
    idl1 = ido * l1
    iw = iw - ( ip - 1 ) * ido
    na = 1 - na

    if ( ip == 4 ) then

      ix2 = iw + ido
      ix3 = ix2 + ido

      if ( na == 0 ) then
        call radf4 ( ido, l1, c, ch, wa(iw), wa(ix2), wa(ix3) )
      else
        call radf4 ( ido, l1, ch, c, wa(iw), wa(ix2), wa(ix3) )
      end if

    else if ( ip == 2 ) then

      if ( na == 0 ) then
        call radf2 ( ido, l1, c, ch, wa(iw) )
      else
        call radf2 ( ido, l1, ch, c, wa(iw) )
      end if

    else if ( ip == 3 ) then

      ix2 = iw + ido

      if ( na == 0 ) then
        call radf3 ( ido, l1, c, ch, wa(iw), wa(ix2) )
      else
        call radf3 ( ido, l1, ch, c, wa(iw), wa(ix2) )
      end if

    else if ( ip == 5 ) then

      ix2 = iw + ido
      ix3 = ix2 + ido
      ix4 = ix3 + ido

      if ( na == 0 ) then
        call radf5 ( ido, l1, c, ch, wa(iw), wa(ix2), wa(ix3), wa(ix4) )
      else
        call radf5 ( ido, l1, ch, c, wa(iw), wa(ix2), wa(ix3), wa(ix4) )
      end if

    else

      if ( ido == 1 ) then
        na = 1 - na
      end if

      if ( na == 0 ) then
        call radfg ( ido, ip, l1, idl1, c, c, c, ch, ch, wa(iw) )
        na = 1
      else
        call radfg ( ido, ip, l1, idl1, ch, ch, ch, c, c, wa(iw) )
        na = 0
      end if

    end if

    l2 = l1

  end do

  if ( na /= 1 ) then
    c(1:n) = ch(1:n)
  end if

  return
end
subroutine rffti ( n, wsave )
!
!*******************************************************************************
!
!! RFFTI initializes WSAVE, used in RFFTF and RFFTB. 
!
!
!  Discussion:
!
!    The prime factorization of N together with a tabulation of the 
!    trigonometric functions are computed and stored in WSAVE.
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!    P N Swarztrauber, 
!    Vectorizing the FFT's, 
!    in Parallel Computations (G. Rodrigue, editor), 
!    Academic Press, 1982, pages 51-83.
!
!    B L Buzbee, 
!    The SLATEC Common Math Library, 
!    in Sources and Development of Mathematical Software (W. Cowell, editor),
!    Prentice Hall, 1984, pages 302-318.
!
!  Parameters:
!
!    Input, integer(4) N, the length of the sequence to be transformed. 
!
!    Output, real(8) WSAVE(2*N+15), contains data, dependent on the value
!    of N, which is necessary for the RFFTF and RFFTB routines.  
!
  integer(4) n
!
  real(8) wsave(2*n+15)
!
  if ( n <= 1 ) then
    return
  end if

  call rffti1 ( n, wsave(n+1), wsave(2*n+1) )

  return
end
subroutine rffti1 ( n, wa, ifac )
!
!*******************************************************************************
!
!! RFFTI1 is a lower level routine used by RFFTI.
!
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
!    Input, integer(4) N, the length of the sequence to be transformed. 
!
!    Input, real(8) WA(N).
!
!    Input, integer(4) IFAC(15).
!    IFAC(1) = N, the number that was factored.
!    IFAC(2) = NF, the number of factors.
!    IFAC(3:2+NF), the factors.
!
  integer(4) n
!
  real(8) arg
  real(8) argh
  real(8) argld
  real(8) fi
  integer(4) i
  integer(4) ib
  integer(4) ido
  integer(4) ifac(15)
  integer(4) ii
  integer(4) ip
  integer(4) is
  integer(4) j
  integer(4) k1
  integer(4) l1
  integer(4) l2
  integer(4) ld
  integer(4) nf
  real(8) pimach
  real(8) wa(n)
!
  call i_factor ( n, ifac )

  nf = ifac(2)

  argh = 2.0D+00 * pimach() / dble ( n )
  is = 0
  l1 = 1

  do k1 = 1, nf-1

    ip = ifac(k1+2)
    ld = 0
    l2 = l1 * ip
    ido = n / l2

    do j = 1, ip-1

      ld = ld + l1
      i = is
      argld = dble ( ld ) * argh
      fi = 0.0D+00

      do ii = 3, ido, 2
        i = i + 2
        fi = fi + 1.0D+00
        arg = fi * argld
        wa(i-1) = cos ( arg )
        wa(i) = sin ( arg )
      end do

      is = is + ido

    end do

    l1 = l2

  end do

  return
end
subroutine rsftb ( n, r, azero, a, b )
!
!*******************************************************************************
!
!! RSFTB computes a "slow" backward Fourier transform of real(8) data.
!
!
!  Modified:
!
!    13 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer(4) N, the number of data values. 
!
!    Output, real(8) R(N), the reconstructed data sequence.
!
!    Input, real(8) AZERO, the constant Fourier coefficient.
!
!    Input, real(8) A(N/2), B(N/2), the Fourier coefficients.
!
  integer(4) n
!
  real(8) a(n/2)
  real(8) azero
  real(8) b(n/2)
  integer(4) i
  integer(4) k
  real(8), parameter :: pi = 3.14159265358979323846264338327950288419716939937510D+00
  real(8) r(n)
  real(8) theta
!
  r(1:n) = azero
  do i = 1, n
    do k = 1, n/2
      theta = dble ( k * ( i - 1 ) * 2 ) * pi / dble ( n )
      r(i) = r(i) + a(k) * cos ( theta ) + b(k) * sin ( theta )
    end do
  end do

  return
end
subroutine rsftf ( n, r, azero, a, b )
!
!*******************************************************************************
!
!! RSFTF computes a "slow" forward Fourier transform of real(8) data.
!
!
!  Modified:
!
!    13 March 2001
!
!  Parameters:
!
!    Input, integer(4) N, the number of data values.
!
!    Input, real(8) R(N), the data to be transformed.
!
!    Output, real(8) AZERO, = sum ( 1 <= I <= N ) R(I) / N.
!
!    Output, real(8) A(N/2), B(N/2), the Fourier coefficients.
!
  integer(4) n
!
  real(8) a(1:n/2)
  real(8) azero
  real(8) b(1:n/2)
  integer(4) i
  integer(4) j
  real(8), parameter :: pi = 3.14159265358979323846264338327950288419716939937510D+00
  real(8) r(n)
  real(8) theta
!
  azero = sum ( r(1:n) ) / dble ( n )

  do i = 1, n / 2

    a(i) = 0.0D+00
    b(i) = 0.0D+00

    do j = 1, n
      theta = dble ( 2 * i * ( j - 1 ) ) * pi / dble ( n )
      a(i) = a(i) + r(j) * cos ( theta )
      b(i) = b(i) + r(j) * sin ( theta )
    end do

    a(i) = a(i) / dble ( n )
    b(i) = b(i) / dble ( n )

    if ( i /= ( n / 2 ) ) then
      a(i) = 2.0D+00 * a(i)
      b(i) = 2.0D+00 * b(i)
    end if

  end do

  return
end
subroutine rvec_random ( alo, ahi, n, a )
!
!*******************************************************************************
!
!! RVEC_RANDOM returns a random real(8) vector in a given range.
!
!
!  Modified:
!
!    04 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real(8) ALO, AHI, the range allowed for the entries.
!
!    Input, integer(4) N, the number of entries in the vector.
!
!    Output, real(8) A(N), the vector of randomly chosen values.
!
  integer(4) n
!
  real(8) a(n)
  real(8) ahi
  real(8) alo
  integer(4) i
!
  do i = 1, n
    call r_random ( alo, ahi, a(i) )
  end do

  return
end
subroutine rvec_reverse ( n, a )
!
!*******************************************************************************
!
!! RVEC_REVERSE reverses the elements of a real(8) vector.
!
!
!  Example:
!
!    Input:
!
!      N = 5, A = ( 11.0, 12.0, 13.0, 14.0, 15.0 ).
!
!    Output:
!
!      A = ( 15.0, 14.0, 13.0, 12.0, 11.0 ).
!
!  Modified:
!
!    06 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer(4) N, the number of entries in the array.
!
!    Input/output, real(8) A(N), the array to be reversed.
!
  integer(4) n
!
  real(8) a(n)
  integer(4) i
!
  do i = 1, n/2
    call r_swap ( a(i), a(n+1-i) )
  end do

  return
end
subroutine sinqb ( n, x, wsave )
!
!*******************************************************************************
!
!! SINQB computes the fast sine transform of quarter wave data. 
!
!
!  Discussion:
!
!    SINQB computes a sequence from its representation in terms of a sine 
!    series with odd wave numbers.
!
!    SINQF is the unnormalized inverse of SINQB since a call of SINQB
!    followed by a call of SINQF will multiply the input sequence X by 4*N.
!
!    The array WSAVE must be initialized by calling SINQI.
!
!    The transform is defined by:
!
!      X_out(I) = sum ( 1 <= K <= N )
!
!        4 * X_in(K) * sin ( ( 2 * K - 1 ) * I * PI / ( 2 * N ) )
!
!  Modified:
!
!    12 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!    P N Swarztrauber, 
!    Vectorizing the FFT's, 
!    in Parallel Computations (G. Rodrigue, editor), 
!    Academic Press, 1982, pages 51-83.
!
!    B L Buzbee, 
!    The SLATEC Common Math Library, 
!    in Sources and Development of Mathematical Software (W. Cowell, editor),
!    Prentice Hall, 1984, pages 302-318.
!
!  Parameters:
!
!    Input, integer(4) N, the length of the array to be transformed.  The 
!    method is more efficient when N is the product of small primes.
!
!    Input/output, real(8) X(N).
!    On input, the sequence to be transformed.
!    On output, the transformed sequence.
!
!    Input, real(8) WSAVE(3*N+15), a work array.  The WSAVE array must be
!    initialized by calling SINQI.  A different WSAVE array must be used 
!    for each different value of N. 
!
  integer(4) n
!
  integer(4) k
  real(8) wsave(3*n+15)
  real(8) x(n)
!
  if ( n < 1 ) then
    return
  end if

  if ( n == 1 ) then
    x(1) = 4.0D+00 * x(1)
    return
  end if

  x(2:n:2) = -x(2:n:2)

  call cosqb ( n, x, wsave )
!
!  Reverse the X vector.
!
  call rvec_reverse ( n, x )

  return
end
subroutine sinqf ( n, x, wsave )
!
!*******************************************************************************
!
!! SINQF computes the fast sine transform of quarter wave data. 
!
!
!  Discussion:
!
!    SINQF computes the coefficients in a sine series representation with 
!    only odd wave numbers. 
!
!    SINQB is the unnormalized inverse of SINQF since a call of SINQF
!    followed by a call of SINQB will multiply the input sequence X by 4*N.
!
!    The array WSAVE, which is used by SINQF, must be initialized by 
!    calling SINQI.
!
!    The transform is defined by:
!
!      X_out(I) = (-1)**(I-1) * X_in(N) + sum ( 1 <= K <= N-1 )
!        2 * X_in(K) * sin ( ( 2 * I - 1 ) * K * PI / ( 2 * N ) )
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!    P N Swarztrauber, 
!    Vectorizing the FFT's, 
!    in Parallel Computations (G. Rodrigue, editor), 
!    Academic Press, 1982, pages 51-83.
!
!    B L Buzbee, 
!    The SLATEC Common Math Library, 
!    in Sources and Development of Mathematical Software (W. Cowell, editor),
!    Prentice Hall, 1984, pages 302-318.
!
!  Parameters:
!
!    Input, integer(4) N, the length of the array to be transformed.  The 
!    method is more efficient when N is the product of small primes.
!
!    Input/output, real(8) X(N).
!    On input, the sequence to be transformed.
!    On output, the transformed sequence.
!
!    Input, real(8) WSAVE(3*N+15), a work array.  The WSAVE array must be
!    initialized by calling SINQI.  A different WSAVE array must be used 
!    for each different value of N. 
!
  integer(4) n
!
  integer(4) k
  real(8) wsave(3*n+15)
  real(8) x(n)
!
  if ( n <= 1 ) then
    return
  end if
!
!  Reverse the X vector.
!
  call rvec_reverse ( n, x )

  call cosqf ( n, x, wsave )

  x(2:n:2) = -x(2:n:2)

  return
end
subroutine sinqi ( n, wsave )
!
!*******************************************************************************
!
!! SINQI initializes WSAVE, used in SINQF and SINQB. 
!
!
!  Discussion:
!
!    The prime factorization of N together with a tabulation of the 
!    trigonometric functions are computed and stored in WSAVE.
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!    P N Swarztrauber, 
!    Vectorizing the FFT's, 
!    in Parallel Computations (G. Rodrigue, editor), 
!    Academic Press, 1982, pages 51-83.
!
!    B L Buzbee, 
!    The SLATEC Common Math Library, 
!    in Sources and Development of Mathematical Software (W. Cowell, editor),
!    Prentice Hall, 1984, pages 302-318.
!
!  Parameters:
!
!    Input, integer(4) N, the length of the array to be transformed.  
!
!    Output, real(8) WSAVE(3*N+15), contains data, dependent on the value
!    of N, which is necessary for the SINQF or SINQB routines.  
!
  integer(4) n
!
  real(8) wsave(3*n+15)
!
  call cosqi ( n, wsave )

  return
end
subroutine sint ( n, x, wsave )
!
!*******************************************************************************
!
!! SINT computes the discrete Fourier sine transform of an odd sequence. 
!
!
!  Discussion:
!
!    SINT is the unnormalized inverse of itself since a call of SINT
!    followed by another call of SINT will multiply the input sequence
!    X by 2*(N+1).
!
!    The array WSAVE must be initialized by calling SINTI.
!
!    The transform is defined by:
!
!      X_out(I) = sum ( 1 <= K <= N ) 
!        2 * X_in(K) * sin ( K * I * PI / ( N + 1 ) )
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!    P N Swarztrauber, 
!    Vectorizing the FFT's, 
!    in Parallel Computations (G. Rodrigue, editor), 
!    Academic Press, 1982, pages 51-83.
!
!    B L Buzbee, 
!    The SLATEC Common Math Library, 
!    in Sources and Development of Mathematical Software (W. Cowell, editor),
!    Prentice Hall, 1984, pages 302-318.
!
!  Parameters:
!
!    Input, integer(4) N, the length of the sequence to be transformed. 
!    The method is most efficient when N+1 is the product of small primes. 
!
!    Input/output, real(8) X(N).
!    On input, the sequence to be transformed.
!    On output, the transformed sequence.
!
!    Input, real(8) WSAVE((5*N+30)/2), a work array.  The WSAVE array must be
!    initialized by calling SINTI.  A different WSAVE array must be used 
!    for each different value of N. 
!
  integer(4) n
!
  integer(4) iw1
  integer(4) iw2
  integer(4) iw3
  real(8) wsave((5*n+30)/2)
  real(8) x(n)
!
  write ( *, * ) 'DEBUG: entered SINT'
  iw1 = n / 2 + 1
  iw2 = iw1 + n + 1
  iw3 = iw2 + n + 1

  call sint1 ( n, x, wsave(1), wsave(iw1), wsave(iw2), wsave(iw3) )
  write ( *, * ) 'DEBUG: leaving SINT'
  return
end
subroutine sint1 ( n, war, was, xh, x, ifac )
!
!*******************************************************************************
!
!! SINT1 is a lower level routine used by SINT.
!
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
!    Input, integer(4) N, the length of the sequence to be transformed. 
!
!    Input/output, real(8) WAR(N).
!    On input, the sequence to be transformed.
!    On output, the transformed sequence.
!
!    Input, real(8) WAS(N/2).
!
!    Input, real(8) XH(N).
!
!    Input, real(8) X(N).
!
!    Input, integer(4) IFAC(15).
!    IFAC(1) = N, the number that was factored.
!    IFAC(2) = NF, the number of factors.
!    IFAC(3:2+NF), the factors.
!
  integer(4) n
!
  integer(4) i
  integer(4) ifac(15)
  integer(4) k
  integer(4) kc
  integer(4) ns2
  real(8), parameter :: sqrt3 = 1.73205080756888D+00
  real(8) t1
  real(8) t2
  real(8) war(n)
  real(8) was(n/2)
  real(8) x(n)
  real(8) xh(n)
  real(8) xhold
!
  write ( *, * ) 'DEBUG: entered SINT1'
  xh(1:n) = war(1:n)
  war(1:n) = x(1:n)

  if ( n <= 1 ) then
    xh(1) = 2.0D+00 * xh(1)
    return
  end if

  if ( n == 2 ) then
    xhold = sqrt3 * ( xh(1) + xh(2) )
    xh(2) = sqrt3 * ( xh(1) - xh(2) )
    xh(1) = xhold
    return
  end if

  ns2 = n / 2
  x(1) = 0.0D+00

  do k = 1, n/2
    t1 = xh(k) - xh(n+1-k)
    t2 = was(k) * ( xh(k) + xh(n+1-k) )
    x(k+1) = t1 + t2
    x(n+2-k) = t2 - t1
  end do

  if ( mod ( n, 2 ) /= 0 ) then
    x(n/2+2) = 4.0D+00 * xh(n/2+1)
  end if
  write ( *, * ) 'DEBUG: calling RFFTF1'
  call rfftf1 ( n+1, x, xh, war, ifac )
  write ( *, * ) 'DEBUG: back from RFFTF1'
  xh(1) = 0.5D+00 * x(1)
  do i = 3, n, 2
    xh(i-1) = -x(i)
    xh(i) = xh(i-2) + x(i-1)
  end do

  if ( mod ( n, 2 ) == 0 ) then
    xh(n) = -x(n+1)
  end if

  x(1:n) = war(1:n)
  war(1:n) = xh(1:n)

  write ( *, * ) 'DEBUG: leaving SINT1'
  return
end
subroutine sinti ( n, wsave )
!
!*******************************************************************************
!
!! SINTI initializes WSAVE, used in SINT.
!
!
!  Discussion:
!
!    The prime factorization of N together with a tabulation of the 
!    trigonometric functions are computed and stored in WSAVE.
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!    P N Swarztrauber, 
!    Vectorizing the FFT's, 
!    in Parallel Computations (G. Rodrigue, editor), 
!    Academic Press, 1982, pages 51-83.
!
!    B L Buzbee, 
!    The SLATEC Common Math Library, 
!    in Sources and Development of Mathematical Software (W. Cowell, editor),
!    Prentice Hall, 1984, pages 302-318.
!
!  Parameters:
!
!    Input, integer(4) N, the length of the sequence to be transformed.  
!    The method is most efficient when N+1 is a product of small primes.
!
!    Output, real(8) WSAVE((5*N+30)/2), contains data, dependent on the value
!    of N, which is necessary for the SINT routine.  
!
  integer(4) n
!
  real(8) dt
  integer(4) k
  real(8) pimach
  real(8) wsave((5*n+30)/2)
!
  if ( n <= 1 ) then
    return
  end if

  dt = pimach() / dble ( n + 1 )

  do k = 1, n/2
    wsave(k) = 2.0D+00 * sin ( dble ( k ) * dt )
  end do

  call rffti ( n+1, wsave((n/2)+1) )

  return
end
function uniform_01_sample ( iseed )
!
!*******************************************************************************
!
!! UNIFORM_01_SAMPLE is a portable random number generator.
!
!
!  Formula:
!
!    ISEED = ISEED * (7**5) mod (2**31 - 1)
!    RANDOM = ISEED * / ( 2**31 - 1 )
!
!  Modified:
!
!    01 March 1999
!
!  Parameters:
!
!    Input/output, integer(4) ISEED, the integer(4) "seed" used to generate
!    the output random number, and updated in preparation for the
!    next one.  ISEED should not be zero.
!
!    Output, real(8) UNIFORM_01_SAMPLE, a random value between 0 and 1.
!

!
!  IA = 7**5
!  IB = 2**15
!  IB16 = 2**16
!  IP = 2**31-1
!
  integer(4), parameter :: ia = 16807
  integer(4), parameter :: ib15 = 32768
  integer(4), parameter :: ib16 = 65536
  integer(4), parameter :: ip = 2147483647
!
  integer(4) iprhi
  integer(4) iseed
  integer(4) ixhi
  integer(4) k
  integer(4) leftlo
  integer(4) loxa
  real(8) uniform_01_sample
!
!  Don't let ISEED be 0.
!
  if ( iseed == 0 ) then
    iseed = ip
  end if
!
!  Get the 15 high order bits of ISEED.
!
  ixhi = iseed / ib16
!
!  Get the 16 low bits of ISEED and form the low product.
!
  loxa = ( iseed - ixhi * ib16 ) * ia
!
!  Get the 15 high order bits of the low product.
!
  leftlo = loxa / ib16
!
!  Form the 31 highest bits of the full product.
!
  iprhi = ixhi * ia + leftlo
!
!  Get overflow past the 31st bit of full product.
!
  k = iprhi / ib15
!
!  Assemble all the parts and presubtract IP.  The parentheses are
!  essential.
!
  iseed = ( ( ( loxa - leftlo * ib16 ) - ip ) + ( iprhi - k * ib15 ) * ib16 ) &
    + k
!
!  Add IP back in if necessary.
!
  if ( iseed < 0 ) then
    iseed = iseed + ip
  end if
!
!  Multiply by 1 / (2**31-1).
!
  uniform_01_sample = real( iseed ) * 4.656612875E-10

  return
end























