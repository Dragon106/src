module hamil

    use consts
    use utils 
    use mesh
    use basis
    use molpot 
    implicit none 
    
    real(8), allocatable :: coeffs(:,:), engs(:) 
 
contains
!-------

subroutine clr_ham()

    deallocate( coeffs, engs )

end subroutine 
!-------------

subroutine ham_diag(ham, Smat, engs)
    use lapack95
    real(8), intent(inout) :: ham(:,:)
    real(8), intent(inout) :: Smat(:,:)
    real(8), intent(out)   :: engs(:)
    
    call sygv(ham, Smat, engs, 1, 'V', 'U')
 
end subroutine 
!-------

subroutine cre_ham() 
	
	allocate( coeffs(nbas,nbas), engs(nbas) ) 

    call calc_ham(coeffs)

end subroutine 
!-------

subroutine calc_ham(ham)
    real(8), intent(out) :: ham(NBAS,NBAS)

    real(8), allocatable :: pote(:,:,:), irad2(:) 
    
    integer :: i1, i2, j1, j2, ji1, ji2, l, ibm  

    allocate( pote(neta,Nbr,Nbr), irad2(nrad) ) 
    
    irad2 = 1.0d0 / rad**2 * wrad(:) 

    pote = 0.0d0
    
    ham = 0.d0

!$omp parallel private(i1,i2, j1,j2, ji1, ji2, l, ibm)

!$omp do 
    do i2 = 1, nbr
    do i1 = ibmin(i2), ibmax(i2)
        do j1 = 1, Neta 
            pote(j1,i1,i2) =  sum( basr(:,i1) * pot(:,j1) * basr(:,i2) * wrad(:) )
        enddo
    enddo
    enddo
!$omp end do    

!$omp do 
    do j2 = 1, nbe
    do j1 = 1, j2 
 
        ji1 = (j1-1)*nbr
        ji2 = (j2-1)*nbr 

        if( j1 == j2 ) then
            l = j1 - 1 + lmin 

            do i2 = 1, nbr
            do i1 = ibmin(i2), i2

! calculate kinetic part             
                ham(i1 + ji1, i2 + ji2) = ham(i1 + ji1, i2 + ji2) &
										- 0.5d0 * ( basr(nrad,i1) * dbasr(nrad,i2) - basr(1,i1) * dbasr(1,i2)) &
                                        + 0.5d0 * sum(dbasr(:,i1) * dbasr(:,i2) * wrad(:) ) 

! calculate centrifugal potential part            
                ham(i1 + ji1, i2 + ji2) = ham(i1 + ji1, i2 + ji2) & 
                                        + sum( basr(:,i1) * irad2(:) * basr(:,i2) ) * (l+1)*l*0.5d0

            enddo
            enddo 
        endif 

! calculate potential part         
        do i2 = 1, nbr
            
            if(j1 == j2) then
                ibm = i2
            else 
                ibm = ibmax(i2)
            endif
            
            do i1 = ibmin(i2), ibm 
                ham(i1 + ji1, i2 + ji2) = ham(i1 + ji1, i2 + ji2) & 
                                        + sum ( weta(:) * base(:,j1) * base(:,j2) * pote(:,i1,i2) )   

            enddo
        enddo 
    enddo
    enddo 
!$omp end do 

!$omp end parallel

    deallocate( pote, irad2 ) 

end subroutine 
!-------

end module 

