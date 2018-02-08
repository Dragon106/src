module cubic_spline
	use type_vars
    implicit none
    
contains
!-------

subroutine cspl(x, y, xi, yi)
    use lapack95
    real(dp), intent(in)  :: x(:), y(size(x)), xi(:)
    real(dp), intent(out) :: yi(size(xi))
    
    real(dp) :: A(size(x),size(x)), C(size(x),0:3)
    real(dp) :: h(size(x)-1), d(size(x)-1)
    real(dp) :: w(size(xi))
    integer :: ipiv(size(x))
    
    integer :: i, Nx
    
    Nx = size(x)
    
    h(1:Nx-1) =  x(2:Nx) - x(1:Nx-1)
    d(1:Nx-1) = (y(2:Nx) - y(1:Nx-1)) / h
    
    A = 0.0d0
    A(1,1)   = 2.0d0
    A(nx,nx) = 2.0d0

    C(1,2)  = 0.0d0
    C(nx,2) = 0.0d0

    do i = 2, Nx-1
        A(i,i-1:i+1) = (/ h(i-1), 2*(h(i-1) + h(i)), h(i) /)
        C(i,2) = 3*(d(i) - d(i-1))
    enddo

    call getrf ( A, ipiv )
    call getrs ( A, ipiv , C(:,2:2) )

    do i = 1, Nx-1
        C(i,0) = y(i)
        C(i,1) = d(i) - h(i)/3.0d0 * (C(i+1,2) + 2*C(i,2))
        C(i,3) = (C(i+1,2) - C(i,2)) / 3.0d0 / h(i)
    enddo
!--
    
!~     where( xi == x(Nx) ) yi = y(Nx) 
!~     do i = 1, Nx-1
!~         where( xi >= x(i) .AND. xi < x(i+1) )
!~             w = xi - x(i)
!~             yi = ( (C(i,3)*w + C(i,2))*w + C(i,1) )*w + C(i,0)
!~         end where
!~     enddo
	where( xi == x(1) )  yi = y(1)
    where( xi == x(Nx) ) yi = y(Nx) 
    do i = 1, Nx-1
        where( xi >= x(i) .AND. xi < x(i+1) )
            w = xi - x(i)
            yi = ( (C(i,3)*w + C(i,2))*w + C(i,1) )*w + C(i,0)
        end where
    enddo
 
end subroutine
!--

end module

