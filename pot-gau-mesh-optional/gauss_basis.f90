! module for calculating Cartesian primitive Gaussian functions 
! 
! in Cartesian coordinate 
! G(x,y,z) = x**nx * y**ny * z**nz * exp( -zeta * (x**2 + y**2 + z**2) ) 
! 
! in spherical coordinate 
! G(r,e,p) = r**(nx+ny+nz) * (1-e*e)**((nx+ny)/2) * e**(nz) * sin(p)**ny * cos(p)**nx * exp( -zeta * r**2 )  

module basis_gauss
	use type_vars
	implicit none 

contains
!*******

pure subroutine gaurad(gau, rad, zeta, nx, ny, nz) 

!-- radial dependent part 

	integer,  intent(in)  :: nx, ny, nz
	real(dp), intent(in)  :: zeta, rad(:) 
	real(dp), intent(out) :: gau(size(rad))  

	gau = rad**(nx+ny+nz) * exp(-zeta*rad*rad)  

end subroutine
!-------------

pure subroutine dgaurad(dgau, rad, zeta, nx, ny, nz) 

!-- derivative radial dependent part  

	integer,  intent(in)  :: nx, ny, nz
	real(dp), intent(in)  :: zeta, rad(:)
	real(dp), intent(out) :: dgau(size(rad))   

	dgau = ( (nx+ny+nz) - 2*zeta*rad*rad ) * rad**(nx+ny+nz-1) & 
		   * exp(-zeta*rad*rad) 
	 
end subroutine
!-------------

pure subroutine gaueta(gau, eta, zeta, nx, ny, nz) 

!-- eta ( cos(theta) ) dependent part 

	integer,  intent(in)  :: nx, ny, nz
	real(dp), intent(in)  :: zeta, eta(:) 
	real(dp), intent(out) :: gau(size(eta)) 

	gau = (1.0_dp - eta*eta)**( (nx+ny)/2.0_dp ) * eta**(nz)   
 
end subroutine 
!-------------

pure subroutine dgaueta(dgau, eta, zeta, nx, ny, nz) 

!-- derivative eta ( cos(theta) ) dependent part 

	integer,  intent(in)  :: nx, ny, nz
	real(dp), intent(in)  :: zeta, eta(:)
	real(dp), intent(out) :: dgau(size(eta)) 

	dgau = ( -eta*eta*(nx+ny+nz) + nz ) & 
		   * (1.0_dp - eta*eta)**( (nx+ny)/2.0_dp - 1 ) * eta**(nz-1) 

end subroutine 
!-------------

pure subroutine gauphi(gau, phi, nx, ny, nz, zeta)
 
!-- phi dependent part

	integer,  intent(in)  :: nx, ny, nz 
	real(dp), intent(in)  :: zeta, phi(:) 
	real(dp), intent(out) :: gau(size(phi)) 

	gau = cos(phi)**nx * sin(phi)**ny 

end subroutine 
!-------------

pure subroutine dgauphi(dgau, phi, zeta, nx, ny, nz) 
 
!-- derivative phi dependent part 

	integer,  intent(in)  :: nx, ny, nz 
	real(dp), intent(in)  :: zeta, phi(:) 
	real(dp), intent(out) :: dgau(size(phi)) 

	dgau = ( -nx*sin(phi)**2 + ny*cos(phi)**2 ) & 
		   * cos(phi)**(nx-1) * sin(phi)**(ny-1)

end subroutine 
!-------------

end module
!*********
