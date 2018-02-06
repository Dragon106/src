module grad_pot
	use mesh
	implicit none
	
	real(8), allocatable, dimension(:,:) :: hpot, deriv_r, deriv_e 
    
contains
!*******

subroutine wrt_grad() 
    real(8) :: grady, gradz
	
	integer :: i, j
	
	allocate( hpot(Nrad,Neta), deriv_r(Nrad,Neta), deriv_e(Nrad,Neta) )

!-- read potential file

	open(1, file="POT.DAT")	
	do i = 1, Neta	
		do j = 1, Nrad
			read(1,*) rad(j), eta(i), hpot(j,i)
		enddo
		read(1,*)
	enddo	
	close(1)

!	write(*,*) 'Reading file done!'
!	write(*,*) size(hpot,1), size(hpot,2)

!-- calculating derivatives  	
	call derivative(deriv_r,deriv_e, hpot, rad,eta, Nrad,Neta)

	open(2, file = "GRADV.DAT")	
		do i = 1, Neta
			do j = 1, Nrad
					grady = sqrt(1-eta(i)**2) * ( deriv_r(j,i) - eta(i)* deriv_e(j,i) / rad(j) )
					gradz = eta(i) * deriv_r(j,i) + (1-eta(i)**2) * deriv_e(j,i) / rad(j)
					write(2,'(4(E16.8,4X))') rad(j), eta(i), grady, gradz
			enddo
			write(2,*)
		enddo
	close(2)
!	write(*,*) 'Writing GRAD.DAT file done!'	
		
	deallocate( hpot, deriv_r, deriv_e )
	
end subroutine
!-------------

!-- subroutine to derivative an array by central, forward and backward formulas
subroutine derivative(dr, de, potential, r, eta, NumR, NumE)	
	integer, intent(in) :: NumR, NumE
	real(8), intent(in) :: potential(NumR,NumE)
	real(8), intent(in) :: r(NumR), eta(NumE)
	real(8), intent(out) :: dr(NumR,NumE), de(NumR,NumE)

	integer :: i

    do i = 1, NumE
        dr(:,i) = deriv( potential(:,i), r(:) )
    enddo
    
    do i = 1, NumR 
        de(i,:) = deriv( potential(i,:), eta(:) )
    enddo
    
end subroutine
!-------------

function deriv(fx,x)
    real(8), intent(in) :: fx(:), x(:)
    real(8) :: deriv(size(fx))
    integer :: N
    real(8) :: alp(size(fx)-2), h(size(fx)-1)

    N = size(fx)

!  using forward and backward derivative formulas at the boundaries 
!    deriv(1) = ( fx(2) - fx(1) ) / ( x(2) - x(1) )
!    deriv(N) = ( fx(N) - fx(N-1) ) / ( x(N) - x(N-1) )
     
!  using central derivative formula O^1
!    deriv(2:N-1) = ( fx(3:N) - fx(1:N-2) ) / ( x(3:N) - x(1:N-2) )
    
!  using central derivative formula O^2 
!   h = x(2:N) - x(1:N-1)    
!   alp = h(2:N-1) / h(1:N-2)    
!   deriv(2:N-1) = (- alp**2 * fx(1:N-2) + fx(2:N-1) * (alp**2 - 1) + fx(3:N) ) &
!	                / ( h(1:N-2) * alp* (alp + 1) )
 
! ** using 3-point derivative using parabola interpolation 
   
   h = x(2:N) - x(1:N-1)   
   
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
