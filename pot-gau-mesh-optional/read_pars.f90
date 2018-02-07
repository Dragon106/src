module read_pars
	use type_vars
    implicit none 

    integer  :: Nrad, Neta, Npgau
    integer  :: imesh
	real(dp) :: Rmax, Rasym
	  
    integer  :: lmax
    
    real(dp) :: alpha, beta
    
    integer  :: iVc

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
    
!-- take V_c into account or not
	read (1,*) hdr
	read (1,*) iVc	
    
    close(1) 

end subroutine 
!-------------

subroutine prt_pars()

	write (*, '("ALPH = ", f11.6)') alpha
	write (*, '("BETA = ", f11.6)') beta
	if (iVc .eq. 0) then
		write (*, '("Vc (Y/N?): ", A)') 'No'
		else
		write (*, '("Vc (Y/N?): ", A)') 'Yes'
	endif
	write (*, '')

end subroutine

end module
!*********
