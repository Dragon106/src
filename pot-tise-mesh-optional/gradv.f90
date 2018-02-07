module grad_pot
	use grid
	implicit none
	
	real(dp), allocatable, dimension(:,:) :: hpot, deriv_r, deriv_e 

contains
!*******

subroutine wrt_gradpot()
    real(dp) :: grady, gradz 
	
	integer :: i, j

	allocate( hpot(Nrad,Neta), deriv_r(Nrad,Neta), deriv_e(Nrad,Neta) )

!-- read potential file

	open(unit = 1, file = './output/POT.DAT')	
	do i = 1, Neta	
		do j = 1, Nrad
			read(1,*) rad(j), eta(i), hpot(j,i)
		enddo
		read(1,*)
	enddo	
	close(1)

!-- calculating derivatives 
 	
	call derivative(deriv_r,deriv_e, hpot, rad, eta, Nrad, Neta)

!-- write gradv in file
	open(unit = 2, file = './output/GRADV.DAT')	
		do i = 1, Neta
			do j = 1, Nrad
					grady = sqrt(1-eta(i)**2) * ( deriv_r(j,i) - eta(i)* deriv_e(j,i) / rad(j) )
					gradz = eta(i) * deriv_r(j,i) + (1-eta(i)**2) * deriv_e(j,i) / rad(j)
					write(2,'(4(E16.8,4X))') rad(j), eta(i), grady, gradz
			enddo
			write(2,*)
		enddo
	close(2)	
		
	deallocate( hpot, deriv_r, deriv_e )

end subroutine
!-------------
 	
!-- subroutine to derivative an array by central, forward and backward formulas
subroutine derivative(dr, de, potential, r, eta, NumR,NumE)
	integer, intent(in) :: NumR, NumE
	real(dp), intent(in) :: potential(NumR,NumE)
	real(dp), intent(in) :: r(NumR), eta(NumE)
	real(dp), intent(out) :: dr(NumR,NumE), de(NumR,NumE)

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
    real(dp), intent(in) :: fx(:), x(:)
    real(dp) :: deriv(size(fx))
    integer :: N     
	real(dp) :: alp(size(fx)-2), h(size(fx)-1)

    N = size(fx)
 
!-- using 3-point derivative using parabola interpolation 
   
	h = x(2:N) - x(1:N-1)   
	   
	deriv(1) = ( -( h(2)**2 + 2.0_dp*h(1)*h(2) ) * fx(1) + ( h(2) + h(1) )**2 * fx(2) - h(1)**2 * fx(3) ) &
				  / ( h(2)*h(1) * ( h(2) + h(1) ) )

	deriv(2:N-1) = (- h(2:N-1)**2 * fx(1:N-2) + ( h(2:N-1)**2 - h(1:N-2)**2 )*fx(2:N-1) + h(1:N-2)**2 * fx(3:N) ) &
				  / ( h(2:N-1)*h(1:N-2) * ( h(2:N-1) + h(1:N-2) ) )
						
	deriv(N) = ( h(N-1)**2 * fx(N-2) - ( h(N-1) + h(N-2) )**2 * fx(N-1) + ( h(N-2)**2 + 2.0_dp*h(N-2)*h(N-1) ) * fx(N) ) &
				  / ( h(N-1)*h(N-2) * ( h(N-1) + h(N-2) ) ) 
						
end function
!-----------

end module
!*********
