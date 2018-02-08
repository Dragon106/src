module linear_interp
	use type_vars
    implicit none
    
contains
!*******

subroutine poly(x, y, xi, yi)
    real(dp), intent(in)  :: x(:), y(size(x)), xi(:)
    real(dp), intent(out) :: yi(size(xi))
    
    integer :: i, Nx
    
    Nx = size(x)

!~     where( xi == x(Nx) ) yi = y(Nx) 
    do i = 1, Nx-1
        where( xi >= x(i) .AND. xi < x(i+1) ) 
            yi = (y(i+1) - y(i))/(x(i+1) - x(i)) * (xi-x(i)) + y(i) 
        end where
    enddo
    
end subroutine
!-------------

end module
!*********
