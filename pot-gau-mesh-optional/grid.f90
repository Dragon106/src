module grid
	use type_vars
	use read_pars
	use consts
	use utils
	implicit none
	 
	real(dp), allocatable :: rad(:), wrad(:)
    real(dp), allocatable :: eta(:), weta(:)
    real(dp), allocatable :: pgau(:)
    
    real(dp), allocatable :: w12(:,:)
    real(dp), allocatable :: wpgau(:)
	
contains
!*******

subroutine create_grids()
	integer :: ir, ie

	allocate ( rad(Nrad), wrad(Nrad), eta(Neta), weta(Neta) )
	allocate ( pgau(npgau), wpgau(npgau) )
	
!-- create Gauss-Legendre for r and eta	
!~ 	call cre_grids_GL(Rmax, rad, wrad, eta, weta)
	
!-- create Gauss-Laguerre for r and Gauss-Jacobi for eta
!~ 	call cre_grids_LJ(Rmax, rad, wrad, eta, weta)
	
	if (imesh .EQ. 1) then
		call cre_grids_GL(Rmax, rad, wrad, eta, weta)
	else
		call cre_grids_LJ(Rmax, rad, wrad, eta, weta)
	endif

	
!-- creta Gauss-Legendre for phi	
	call grid_GL(pgau, wpgau, 0.d0, 2.d0*pi)
!--	
	
	if( allocated(w12) ) &
        call stop_error('Mesh weights created before!!!')
	
	allocate ( w12(Nrad,Neta) )
	
 	forall (ir = 1:nrad, ie=1:neta) 
 		w12(ir,ie) = wrad(ir)*weta(ie)
 	end forall

!-- write in files Rinfo.in and Einfo.in    
	call wrt_grids()
	
end subroutine
!-------------

subroutine prt_grids()

	write (*,'(a)')   '=== GRIDS INFO ==='
    write (*, '("NRAD = ", i11)') nrad
    write (*, '("NETA = ", i11)') neta 
    write (*, '("RMAX = ", f11.6 )') rmax
     
	write (*, '("NPHI = ", i11)') npgau
!~ 	write (*, '("ALPH = ", f11.6)') alpha
!~ 	write (*, '("BETA = ", f11.6)') beta
	
    write (*, '')

end subroutine
!-------------

subroutine wrt_grids()
	integer :: i
	
	open (unit=1, file ='./output/Rinfo.in')
        write(1,'(2f)') (rad(i), wrad(i), i = 1, nrad)
    close(1)	

    open (unit=2, file = './output/Einfo.in')
    	write(2,'(3f)') (eta(i),  sqrt( 1 - eta(i)**2 ), weta(i), i = 1, neta)
    close(2)
    
end subroutine
!-------------

subroutine clr_grids()

	deallocate ( rad, wrad, eta, weta, w12 )
	deallocate ( pgau, wpgau )

end subroutine
!-------------

subroutine cre_grids_GL(Rmax, x, wx, e, we)

    real(dp), intent(in) :: Rmax
    real(dp), intent(out) :: x(:), wx(size(x)), e(:), we(size(e))

!~     integer :: i
	write (*,'(a)') 'Gauss-Legendre for r and eta'
	write (*,'')
	
    call grid_GL(x, wx, 0.d0, Rmax)
		wx = wx * x*x
    
    call grid_GL(e, we, -1.d0, 1.d0)

!~     open (unit=1, file ='Rinfo.in')
!~         write(1,'(2f)') (x(i), wx(i), i = 1, size(x))
!~     close(1)	

!~     open(unit=2, file = 'Einfo.in')
!~     	write(2,'(3f)') (e(i),  sqrt( 1 - e(i)**2 ), we(i), i = 1, size(e))
!~     close(2)

end subroutine
!-------------

subroutine cre_grids_LJ(Rmax, x, wx, e, we)
	real(dp), intent(in) :: Rmax
	real(dp), intent(out) :: x(:), wx(size(x)), e(:), we(size(e))
	
	real(dp), allocatable :: xacs(:), weig(:)
	integer :: Nr, Ne
	integer :: itype, ifail
	
	Nr = size(x)
	Ne = size(e)
	
	write (*,'(a)') 'Gauss-Laguerre and Jacobi'
	write (*,'')
	
!-- create Gauss-Laguerre for r

	allocate ( xacs(nr), weig(nr) ) 
	
    itype = -3
    ifail = 0
    call d01bcf(itype, 0.0d0, 1.0d0, 1.0d0, 0.0d0, nr, weig, xacs, ifail)
		x = xacs / (xacs(nr) / Rmax)
		wx = weig * rad*rad * Rmax/xacs(nr)
	
	deallocate(xacs, weig)
	
!-- create Jacobi grid for eta
	allocate(xacs(ne), weig(ne))

	itype = 1
	call d01bcf(itype, -1.d0, 1.d0, 0.d0, 0.d0, ne, weig, xacs, ifail)

		e  = xacs * sqrt(2.0_dp - xacs*xacs)
		we = weig * 2.0_dp*(1 - xacs*xacs) / sqrt(2.0_dp - xacs*xacs)

	deallocate(xacs, weig)
	
end subroutine
!-------------

subroutine grid_GL(x, w, xi, xf)
	real(dp), intent(in) :: xi, xf
	real(dp), intent(out) :: x(:), w(size(x))
	integer :: N, itype

	N = size(x)
    itype = 1
    call cgqf ( N, itype, 0.0d0, 0.0d0, xi, xf, x, w )

end subroutine
!-------------

subroutine read_grids()
	real(dp), allocatable :: dsn(:)
	integer :: i

	allocate (dsn(neta))
	open(unit=1, file = 'Rinfo.in')
		read(1,*) (rad(i), wrad(i), i = 1, nrad)
	close(1)

	open(unit=2, file = 'Einfo.in')
		read(2, *) (eta(i), dsn(i) , weta(i), i = 1, neta) 
	close(2)

end subroutine
!-------------

end module
!*********
