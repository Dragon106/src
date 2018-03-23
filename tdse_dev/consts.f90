module consts
	use type_vars
	
	implicit none
	
	complex(wp), parameter :: Im = cmplx(0.0_wp, 1.0_wp)
	real(wp), parameter :: Pi = acos(-1.0_wp)
	real(wp), parameter :: deg2rad = Pi/180.0_wp
	real(wp), parameter :: rad2deg = 180.0_wp/Pi
	real(wp), parameter :: auI = 351.3393_wp
	real(wp), parameter :: auE = 27.21138457_wp
	real(wp), parameter :: auT = 0.024189_wp
	
end module
