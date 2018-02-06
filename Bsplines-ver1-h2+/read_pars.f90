module read_pars
    use basis
    use molpot
    use wfs
    implicit none 

    character :: potfln*100

contains
!*******

subroutine read_input(fln)
	character, intent(in) :: fln*(*) 

    character :: sdum*100 

    open(unit = 1, file = fln, status = 'old')

    read (1,*) sdum
    read (1,*) nrad, neta, rmax 

    read (1,*) sdum
    read (1,*) nbr, deg, gam, iline

    read (1,*) sdum
    read (1,*) lmax, mq
    
    read (1,*) sdum
    read (1,*) Emax
    
    read (1,*) sdum
    read (1,*) potfln
    
    read (1,*) sdum
    read (1,*) R
    
    close(1) 

end subroutine 
!-------------

end module 
!*********
