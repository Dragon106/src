module derivfx_m

    implicit none

contains
!*******
! https://en.wikipedia.org/wiki/Finite_difference_coefficient
!*******
subroutine deriv(dx, fx, dfx)
    integer, parameter :: npts = 5
    real(8), parameter :: cent_coeff(npts)   = (/  1.d0/12,    -2.d0/3,    0.d0,    2.d0/3, -1.d0/12 /)
    real(8), parameter :: forw_coeff(npts-1) = (/ -11.d0/6,     3.d0,     -3.d0/2,  1.d0/3  /)
    real(8), parameter :: back_coeff(npts-1) = (/ -1.d0/3,      3.d0/2,   -3.d0,    11.d0/6 /)

    real(8), intent(in)  :: dx, fx(:)
    real(8), intent(out) :: dfx(size(fx))

    integer :: nhalf = npts/2
    integer :: i, n

    N = size(fx)
    if ( n < npts ) then
        write(*,*) 'The size of fx is smaller than npts!'
        stop
    endif

!-- at left boundaries using forward formula
    forall (i = 1 : nhalf)
        dfx(i) = sum( forw_coeff * fx(i : i+npts-2) ) / dx
    end forall

!-- at right boundaries using backward formula
    forall (i = n-nhalf+1 : n)
        dfx(i) = sum( back_coeff * fx(i-npts+2 : n) ) / dx
    end forall

!-- using central formula at (nhalf+1 : n-nhalf)
    forall (i = nhalf+1 : n-nhalf)
        dfx(i) = sum( cent_coeff * fx(i-nhalf : i+nhalf) ) / dx
    end forall

end subroutine
!-------------

subroutine second_deriv(dx, fx, d2fx)
    integer, parameter :: npts = 5
    real(8), parameter :: cent_coeff(npts)   = (/ -1.d0/12, 4.d0/3, -5.d0/2, 4.d0/3, -1.d0/12 /)
    real(8), parameter :: forw_coeff(npts-1) = (/  2.d0,   -5.d0,    4.d0,  -1.d0 /)
    real(8), parameter :: back_coeff(npts-1) = (/ -1.d0,    4.d0,   -5.d0,   2.d0 /)

    real(8), intent(in)  :: dx, fx(:)
    real(8), intent(out) :: d2fx(size(fx))

    integer :: nhalf = npts/2
    integer :: i, n
    real(8) :: h2

    n = size(fx)
    h2 = dx*dx
    if ( n < npts ) then
        write(*,*) 'The size of fx is smaller than npts!'
        stop
    endif

!-- at left boundaries
    forall (i = 1 : nhalf)
        d2fx(i) = sum( forw_coeff * fx(i : i+npts-2) ) / h2
    end forall

!-- at right boundaries
    forall (i = n-nhalf+1 : n)
        d2fx(i) = sum( back_coeff * fx(i-npts+2 : i) ) / h2
    end forall

!-- using central formula at (nhalf+1 : n-nhalf)
    forall (i = nhalf+1 : n-nhalf)
        d2fx(i) = sum( cent_coeff * fx(i-nhalf : i+nhalf) ) / h2
    end forall

end subroutine
!-------------

end module
!*********
