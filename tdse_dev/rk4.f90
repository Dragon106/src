module rk4_m
	use type_vars
	use tdse_m
	implicit none
	
    private 

    integer,  parameter :: order = 4
    real(wp), parameter :: nodes(order-1) = [1/2._wp, 1/2._wp, 1.0_wp]
    real(wp), parameter :: weig(order)    = [1/6._wp, 1/3._wp, 1/3._wp, 1/6._wp]
    real(wp), parameter :: coeff(order-1) = [1/2._wp, 1/2._wp, 1.0_wp]

	type :: rk4_t
		integer :: neq 
		complex(wp), allocatable :: kcoef(:,:)
		complex(wp), allocatable :: ytmp(:)
		
		contains
		procedure :: pre => pre_rk4
		procedure :: update => update_rk4
		procedure :: destroy => destroy_rk4
	end type
	
	public :: rk4_t
	
contains
!*******

subroutine pre_rk4(this, neq)
	class(rk4_t), intent(inout) :: this
	integer,      intent(in)    :: neq

    call destroy_rk4(this)
 	this%neq = neq
 	   
	allocate( this%kcoef(neq,4) )
	allocate( this%ytmp(neq)    )
	
end subroutine
!-------------

subroutine destroy_rk4(this)
	class(rk4_t), intent(inout) :: this
	
	this%neq = 0
	
	if (allocated(this%kcoef)) deallocate(this%kcoef)
	
	if (allocated(this%ytmp) ) deallocate (this%ytmp)
	
end subroutine
!-------------

subroutine update_rk4(this, ode, y, t1, t2)
	class(rk4_t), intent(inout) :: this
	complex(wp),  intent(inout) :: y(this%neq)
	real(wp),     intent(in)    :: t1, t2
	type(tdse_t), intent(inout) :: ode

	real(wp) :: h
	integer  :: ik, i 

    call ode%derivfunc( this%kcoef(:,1), t1, y ) 

	h = t2-t1
	
    do ik = 2, order 
      !$OMP PARALLEL DO SIMD
      !!$OMP SIMD 
      do i = 1, this%neq
        this%ytmp(i) = y(i) + h*(coeff(ik-1) * this%kcoef(i,ik-1) ) 
      enddo
      !!$OMP END SIMD 
      !$OMP END PARALLEL DO SIMD 
      
      call ode%derivfunc( this%kcoef(:,ik), t1 + nodes(ik-1)*h, this%ytmp )
      
	enddo
	
	!$OMP PARALLEL PRIVATE(IK, I)
	do ik = 1, order
	    !$OMP DO
	    do i = 1, this%neq 
		  y(i) = y(i) + h * weig(ik) * this%kcoef(i,ik)
		enddo 
		!$OMP END DO 
	enddo
	!$OMP END PARALLEL
	
end subroutine
!-------------

end module
