module mesh 

    use utils
    implicit none

    integer :: Nrad, Neta   
    real(8) :: rmax 
    
    real(8), allocatable :: rad(:), wrad(:)
    real(8), allocatable :: eta(:), weta(:)
    real(8), allocatable :: w12(:,:)

contains
!*******

subroutine prt_mesh()

    write (*, '(a)' ) '=== MESH INFO ==='

    write (*, '("NRAD = ", i11)') nrad
    write (*, '("NETA = ", i11)') neta 
    write (*, '("RMAX = ", f11.6 )') rmax 

    write (*, '')

end subroutine 
!-------------

subroutine cre_mesh()

    integer :: ir, ie
     
    if( allocated(rad) .or. allocated(wrad) .or. &
        allocated(eta) .or. allocated(weta) ) & 
        call stop_error('Mesh created before!!!')
	
	allocate( rad(nrad), wrad(nrad), eta(neta), weta(neta) )

!-- Gauss-Laguerre and Jacobi  
!   call mesh_gauqua_1(rad, eta, wrad, weta)

!-- Gauss-Legendre (both)
    call mesh_quadrule_1(rad, eta, wrad, weta) 

    if( allocated(w12) ) &
        call stop_error('Mesh weights created before!!!')  

    allocate( w12(nrad,neta) )

    forall ( ir = 1:Nrad, ie = 1:Neta )    
        w12(ir,ie) = rad(ir)**2 * wrad(ir) * weta(ie)
    end forall

end subroutine
!-------------

!-- QUADRULE.F90 -- GAUSS-LEGENDRE
subroutine mesh_quadrule_1(rad, eta, wrad, weta) 

    real(8), intent(out) :: rad(:), eta(:), wrad(:), weta(:)

    real(8), allocatable :: xacs(:), weig(:) 
    integer :: ikind
    integer :: i   

	write(*,'')
    write(*,'(a)') 'Gauss-Legendre for both r and eta!' !(Just writing here temporarily!!!)
    
!-- create radial grid
    allocate( xacs(nrad), weig(nrad) )

    ikind = 1
    call cgqf ( nrad, ikind, 0.D0, 0.D0, 0.d0, rmax, xacs, weig ) 

    rad  = xacs 
    wrad = weig 

	!-- write Rinfo.in for TDSE code!!!!!

    open(unit = 1, file='Rinfo.in') 
        write(1,'(2E24.16)') ( rad(i), weig(i), i = 1, Nrad )
    close(1)   
    
    deallocate( xacs, weig )

!-- create eta grid as Gauss-Legendre (ikind = 1)/Gauss-Jacobi (ikind = 4)
    allocate( xacs(neta), weig(neta) )

!    ikind = 4
!    call cgqf ( Neta, ikind, 0.d0, 0.d0, -1.d0, 1.d0, xacs, weig )
!    eta  = xacs * sqrt(2.0d0 - xacs*xacs)
!    weta = weig * 2.0d0*(1 - xacs*xacs) / sqrt(2.0d0 - xacs*xacs)

    ikind = 1
    call cgqf ( Neta, ikind, 0.d0, 0.d0, -1.d0, 1.d0, xacs, weig )
    eta  = xacs
    weta = weig
    
    !-- write Einfo.in for TDSE!!!!!
    open(unit = 1,file='Einfo.in')
        write(1,'(3E24.16)') (eta(i), 1.0d0 - xacs(i)**2, weig(i), i = 1, Neta)
    close(1)

    deallocate( xacs, weig )

end subroutine
!-------------

!-- GAUQUA.F -- GAUSS-LAGUERRE AND JACOBI
subroutine mesh_gauqua_1(rad,eta,wrad,weta)

    real(8), intent(out) :: rad(:), eta(:), wrad(:), weta(:)

    real(8), allocatable :: xacs(:), weig(:) 
    integer :: itype, ifail
    integer :: i 
    
    write(*,'')
    write(*,'(a)') 'Gauss-Laguerre and Gauss-Jacobi!' !(Just writing here temporarily!!!)

!-- create radial grid as Gauss-Laguerre
    allocate( xacs(nrad), weig(nrad) )

    itype = -3
    call d01bcf(itype,0.d0,1.d0,1.d0,0.d0,nrad,weig,xacs,ifail)

    rad  = xacs * rmax/xacs(nrad)
    wrad = weig * rmax/xacs(nrad)

	!-- write Rinfo.in for TDSE code!!!!!

    open(unit = 1, file='Rinfo.in') 
        write(1,'(3E24.16)') ( xacs(i), rad(i), weig(i), i = 1, Nrad )
    close(1)   

    deallocate(xacs, weig)

!-- create eta grid as Gauss-Jacobi
    allocate( xacs(neta), weig(neta) )

    itype = 1
    call d01bcf(itype,-1.d0,1.d0,0.d0,0.d0,neta,weig,xacs,ifail)

    eta  = xacs * sqrt(2.0d0 - xacs*xacs)
    weta = weig * 2.0d0*(1 - xacs*xacs) / sqrt(2.0d0 - xacs*xacs)

	!-- write Einfo.in for TDSE!!!!!
	
    open(unit = 1,file='Einfo.in')
        write(1,'(4E24.16)') (xacs(i), eta(i), 1.0d0 - xacs(i)**2, weig(i), i = 1, Neta)
    close(1)

    deallocate(xacs, weig)

end subroutine 
!-------------

subroutine clr_mesh()
    
    nrad = 0
    neta = 0
    Rmax = 0.d0
    
    deallocate ( rad, wrad )
    deallocate ( eta, weta )
    deallocate ( w12 )

end subroutine
!-------------

subroutine wrt_mesh()

    integer :: i 

!-- write out file radial points
    open(unit = 1, file ='rad.dat')
    write(1,'(i,2f)') ( i, rad(i), wrad(i), i = 1, nrad )
    close(1)

!-- write out file eta points
    open(unit = 1, file ='eta.dat')
    write(1,'(i,2f)') ( i, eta(i), weta(i), i = 1, neta )
    close(1)

end subroutine
!-------------

end module
!*********

