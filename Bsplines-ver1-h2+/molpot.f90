module molpot

    use mesh
    implicit none

    real(8), allocatable :: pot(:,:)
    real(8), allocatable :: grady_v(:,:), gradz_v(:,:)
    
    real(8) :: R

contains
!*******

subroutine cre_pot()

    allocate ( pot(nrad, neta) )
    allocate ( grady_v(nrad, neta), gradz_v(nrad, neta) )

end subroutine 
!-------------

subroutine clr_pot()

    deallocate ( pot ) 
    deallocate ( grady_v, gradz_v )

end subroutine 
!-------------

subroutine wrt_pot(fln) 
    character, intent(in) :: fln*(*) 

    integer :: i, j 

    open(unit = 1, file = fln) 

    do j = 1, neta  
        do i = 1, nrad
            write(1,'(3f)') rad(i), eta(j), pot(i,j)
        enddo
        write(1,*)
    enddo

    close(1)

end subroutine 
!-------------

subroutine read_pot(fln)
    character, intent(in) :: fln*(*) 

    real(8) :: dum
    integer :: i, j

    open(unit = 2, file = fln, status = 'old')

        do i = 1, Nrad    
            do j = 1, Neta 
                read(2,*) dum, dum, pot(i,j) 
            enddo        
            read(2,*)
        enddo

    close(2)

end subroutine 
!-------------

end module 

