module read_pars
	use type_vars
    implicit none 

    integer  :: imesh
    integer  :: Nrad, Neta, Npgau
	real(dp) :: Rmax, Rasym
	  
    integer  :: lmax
    
    real(dp) :: alpha, beta

contains
!*******

subroutine read_input()
    character :: hdr*100 

    open(unit = 1, file = 'pars.in', status = 'old')

!-- read parameters of box
    read (1,*) hdr
    
    read (1,*) hdr
    read (1,*) imesh
    
    read (1,*) hdr
    read (1,*) nrad, neta, rmax, rasym
    
    read (1,*) hdr
    read (1,*) npgau
     
!-- read lmax for Legendre expansion
    read (1,*) hdr
	read (1,*)
    read (1,*) lmax

!-- read parameters for exchange potential    
    read (1,*) hdr
    read (1,*)
    read (1,*) alpha, beta
    
    close(1) 

end subroutine 
!-------------

end module
!*********
