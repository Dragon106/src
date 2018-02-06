module basis

    use consts
    use utils
    use mesh 
    implicit none

    integer :: lmin, lmax, mq 

    integer :: nbps, deg   
    real(8) :: gam 
    integer :: iline

    integer :: Nbe, Nbr, Nbas 

    real(8), dimension(:,:), allocatable :: basr, base, dbasr, dbase 
    real(8), dimension(:,:), allocatable :: ovlpr, ovlpe, ovlp, ovlpdr

    integer, dimension(:),  allocatable  :: ibmin, ibmax

contains
!-------

include 'basis_funcs.f90'
!-------

subroutine prt_basis()

    write (*, '(a)') '=== BASIS INFO ==='  
 
    write (*, '("NBR  = ", i11)') nbr
    write (*, '("NBE  = ", i11)') nbe  
    write (*, '("NBAS = ", i11)') nbas 
    write (*, '("MQ   = ", i11)') mq 

    write (*, '')

end subroutine  
!-------

subroutine cre_basis()

    real(8), allocatable :: bps(:)     
    integer :: i 
     
    Nbps = Nbr - deg + 2 + 2
!~     Nbps = Nbr - deg + 2 + 1
    
    allocate ( bps(Nbps) )
    
    call cre_bps_exp(bps, 0.d0, rmax, gam) !!! easy for TuTu control!!!
!   call cre_bps_mix(bps, 0.d0, rmax, iline) !!! Too hard to control!!!
!   call cre_bps_line(bps, 0.d0, rmax) !!! Not good!
    
    call cre_bspline(bps, deg, rad, basr, dbasr)

    deallocate( bps )

    lmin = abs(mq) 
    Nbe = lmax - lmin + 1

    allocate( base(Neta,Nbe), dbase(Neta,Nbe) ) 

    do i = lmin, lmax
        call cre_nplg(i,mq,eta,base(:,i-lmin+1))
    enddo

    Nbas = Nbr*Nbe

    allocate( ibmin(nbr), ibmax(nbr) )

    do i = 1, nbr
        ibmin(i) = max(i - deg + 1, 1  ) 
        ibmax(i) = min(i + deg - 1, nbr)
    enddo 

end subroutine 
!-------

subroutine wrt_basis()

    integer :: i 

    open(unit = 1, file = 'basisr.dat') 
    do i = 1, nrad
        write(1, '(*(e))') rad(i), basr(i,:) 
    enddo
    close(1)

    open(unit = 1, file = 'dbasisr.dat') 
    do i = 1, nrad
        write(1, '(*(e))') rad(i), dbasr(i,:) 
    enddo
    close(1)
    
    open(unit = 1, file = 'basise.dat')
    do i = 1, neta  
        write(1, '(*(e))') eta(i), base(i,:) 
    enddo 
    close(1) 

end subroutine 
!-------

subroutine clr_basis()

    deallocate( basr, dbasr )
    deallocate( base, dbase )
    deallocate( ibmin, ibmax )

end subroutine 
!-------

subroutine cre_overlap()

    integer :: i1, i2, i, j

    allocate( ovlpr(Nbr,Nbr), ovlpe(Nbe,Nbe), ovlpdr(Nbr,Nbr) )
    allocate( ovlp(nbas,nbas) )

    do j = 1, nbr
		do i = 1, j
            ovlpr(i,j) = sum ( wrad(:) * basr(:,j) * basr(:,i) )     
		enddo
    enddo

    do j = 1, nbr
		do i = 1, j
            ovlpdr(i,j) = sum ( wrad(:) * dbasr(:,j) * dbasr(:,i) )     
		enddo
    enddo
    
    do j = 1, nbe 
		do i = 1, j
            ovlpe(i,j) = sum ( weta(:) * base(:,j) * base(:,i) )
		enddo
    enddo 
    
    call mat_symz(ovlpe)
    call mat_symz(ovlpr)  

    do j = 1, nbe
        do i2 = 1, nbr
			do i1 = 1, i2 
                ovlp(i1 + (j-1)*nbr, i2 + (j-1)*nbr) = ovlpr(i1,i2)  
			enddo
        enddo
    enddo
    
    call mat_symz(ovlp)         
        
end subroutine
!-------

subroutine wrt_overlap()

    integer :: i, j

    open(unit = 1, file = 'ovlpe.dat')
    do i = 1, Nbe
        do j = 1, Nbe
            write(1,'(2i,e)') i, j, ovlpe(j,i)
        enddo
        write(1,*)
    enddo     
    close(1)
    
    open(unit = 1, file = 'ovlpr.dat')
    do i = 1, Nbr
        do j = 1, Nbr
            write(1,'(2i,e)') i, j, ovlpr(j,i)
        enddo
        write(1,*)
    enddo     
    close(1)

    open(unit = 1, file = 'ovlpdr.dat')
    do i = 1, Nbr
        do j = 1, Nbr
            write(1,'(2i,e)') i, j, ovlpdr(j,i)
        enddo
        write(1,*)
    enddo     
    close(1)
    
end subroutine
!-------

subroutine clr_overlap()

    deallocate( ovlpr, ovlpe )
    deallocate( ovlp )

end subroutine 
!-------

subroutine cre_bps_exp(bps, rmin, rmax, gam)
    real(8), intent(in) :: rmin, rmax, gam
    real(8), intent(out) :: bps(:)
    integer :: N
    
    integer :: i
    real(8) :: dr, e1, ems
    
!    write(*,*)
    write(*,'(a,f6.2)') 'Exponential bps, gamma = ', gam !(Need writing another place! Think about it! and about him) 
    write(*,*)
    
    N = size(bps)
    
    dr  = rmax - rmin
    e1  = exp( gam/(N-1) )
    ems = exp(gam) - 1    
    
    open(unit=1, file = 'exp-bps.dat')
!    tmp = 1/e1
    do i = 1, N 
        bps(i) = rmin + dr * ( e1**(i-1) - 1 ) / ems
        write(1,'(i,f)') i, bps(i)         
!       tmp = tmp*e1
!       bps = ......     * (tmp - 1) / ...
    enddo
    close(1)
 
end subroutine
!-------

subroutine cre_bps_mix(bps, rmin, rmax, i0) ! P.113
    integer, intent(in) :: i0
    real(8), intent(in) :: rmin, rmax
    real(8), intent(out) :: bps(:)

    integer :: N

    integer :: i
    real(8) :: alpha, beta, r0

    N = size(bps)
!	i0 = ceiling(N/1.5d0) ! Need to check
!   i0 = 80    
    r0 = (rmax * (i0-1) + rmin * (N-i0)) / (2*N - i0 - 1)
    alpha = (r0 - rmin) / (i0-1)**2
    beta = (rmax - r0) / (N - i0)

    open(unit=1, file ='mix-bps.dat')

	do i = 1, N
		if (i < i0) then
            bps(i) = rmin + alpha * (i-1)**2
		else
            bps(i) = r0 + beta * (i - i0)
		endif

        write(1,'(i,f)') i, bps(i) 

	enddo

    close(1) 

end subroutine
!-------

subroutine cre_bps_line(bps, rmin, rmax)
    real(8), intent(in) :: rmin, rmax
    real(8), intent(out) :: bps(:)

    integer :: N

    integer :: i
    real(8) :: drn

    N = size(bps)
!	write(*,*) N

    drn = (rmax - rmin) / (N-1)
	do i = 1, N
        bps(i) = rmin + drn * (i-1)
	enddo
    
end subroutine
!-------

subroutine mat_symz(A)
    real(8), intent(inout) :: A(:,:)

    integer :: sz, i, j
    
    sz = min( size(A,1), size(A,2) )
    
    do j = 1, sz-1
        do i = j+1, sz
            A(i,j) = A(j,i)
        enddo 
    enddo 
    
end subroutine 
!-------

end module 

