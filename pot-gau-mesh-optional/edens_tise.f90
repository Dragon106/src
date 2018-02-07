!***********************************************************************
! CALCULATING ELECTRON AND SPIN DENSITY FROM BSPLINE TISE WAVEFUNCTIONS !
!***********************************************************************

module edens_tise
	use type_vars
	use consts
	use grid
	use eledens_info
	use utils
	implicit none

contains
!*******

subroutine calc_dens_tise()

	call read_MOs_wf()
	call dens_tise()
	call write_edens()
	
end subroutine
!-------------

subroutine read_MOs_wf()
	integer :: ir, ie, io
	
	allocate( wf(Nrad, Neta, Norb) )
	allocate( drwf(Nrad, Neta, Norb), dewf(Nrad, Neta, Norb) )
	
	open(unit=1, file = './input-tise/WF.in',  status='old')
	open(unit=2, file = './input-tise/DWF.in', status='old')

	do io = 1, Norb
		do ie = 1, Neta
			do ir = 1, Nrad
!~ 				read(1,*) rad(ir), eta(ie), wf(ir,ie,io)
				read(1,*) wf(ir,ie,io)
				read(2,*) drwf(ir,ie,io), dewf(ir,ie,io)
			enddo
		enddo
		read(1,*)
		read(2,*)
	enddo
	close(1)
	close(2)

end subroutine	
!-------------

subroutine dens_tise()
	integer :: io, ir, ie
	real(dp), allocatable :: keta(:,:)
	
!~ 	real(dp), allocatable :: dsn(:), ngr_seden(:,:), nge_seden(:,:)
	
    den1e(:,:,:) = wf(:,:,:)**2 / (2*pi)
    eden(:,:) = 2 * sum( den1e(:,:,1:Norb-1), 3) + ( 2 - mod(Nele,2) ) * den1e(:,:,Norb)
	seden = sum (den1e, 3)
	
	write(*,'(a)')   '-------------------' 
	write(*,'(a,f)') 'Number of electrons: ', sum( w12*eden ) * 2*pi
		
	allocate ( keta(nrad,neta) )
	
	forall (ir = 1:nrad, ie = 1:neta)
		keta(ir,ie) = - sqrt (1-eta(ie)**2) / rad(ir)
	endforall	
	
	gr_seden = 2*sum(wf*drwf, 3) / (2*pi)
	ge_seden = 2*sum(wf*dewf, 3) * keta / (2*pi)
	
	deallocate (keta, drwf, dewf)
	
!-- calc numerically derivative of electron spin density
!~     allocate ( dsn(neta), ngr_seden(Nrad, Neta), nge_seden(Nrad, Neta) )
    
!~     do ie = 1, neta
!~         ngr_seden(:,ie) = deriv( seden(:,ie), rad(:) )
!~     enddo	
    
!~     dsn = sqrt (1 - eta*eta)
!~     do ir = 1, nrad
!~         nge_seden(ir,:) = - dsn(:) * deriv( seden(ir,:), eta(:)) / rad(ir)
!~     enddo  

!~ !--
!~ 	open(unit = 100, file = './output/ngrad.dat')
!~ 	do ie = 1, neta
!~ 		do ir = 1, nrad
!~ 			write(100, '(4(e,3x))') eden(ir,ie), seden(ir,ie), & 
!~ 									ngr_seden(ir,ie), nge_seden(ir,ie)  
!~ 		enddo
!~ 		write(100, '')
!~ 	enddo

!~ 	close(100)
!~ !--
!~ 	deallocate (dsn, ngr_seden, nge_seden)

end subroutine
!-------------

function deriv(fx, x)

    real(dp), intent(in) :: fx(:), x(:)
    real(dp) :: deriv(size(fx))
    
    integer :: N  
    real(dp) :: h(size(fx) - 1)
    
    N = size(fx)
    
    h(:) = x(2:N) - x(1:N-1)
    
! ** using 3-point derivative using parabola interpolation 
    
    deriv(1) = ( -( h(2)**2 + 2.d0*h(1)*h(2) ) * fx(1) + ( h(2) + h(1) )**2 * fx(2) - h(1)**2 * fx(3) ) &
                    / ( h(2)*h(1) * ( h(2) + h(1) ) )

    deriv(2:N-1) = (- h(2:N-1)**2 * fx(1:N-2) + ( h(2:N-1)**2 - h(1:N-2)**2 )*fx(2:N-1) + h(1:N-2)**2 * fx(3:N) ) &
                    / ( h(2:N-1)*h(1:N-2) * ( h(2:N-1) + h(1:N-2) ) )
                    
    deriv(N) = ( h(N-1)**2 * fx(N-2) - ( h(N-1) + h(N-2) )**2 * fx(N-1) + ( h(N-2)**2 + 2.d0*h(N-2)*h(N-1) ) * fx(N) ) &
                    / ( h(N-1)*h(N-2) * ( h(N-1) + h(N-2) ) )

end function
!-----------

end module
!*********
