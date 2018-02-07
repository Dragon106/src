module integ_weigs_m
	use pars
	implicit none
	
contains
!*******

subroutine integ_weigs(x, weig)
	real(dp), intent(in)  :: x(:)
	real(dp), intent(out) :: weig(size(x))
	
	integer :: n
	
	n = size(x) 
	
	if( n < 2 ) stop 'integ_weigs error: must have at least 2 points.'
	
	weig(1) = 0.5_dp * (x(2) - x(1)  )
	weig(n) = 0.5_dp * (x(n) - x(n-1))
	
	if( n > 2 ) & 
		weig(2:n-1) = 0.5_dp * ( x(3:n) - x(1:n-2) )
	
end subroutine
!-------------
	
end module
