module solve_tdse_m
	use type_vars
	use laser_m
	use data_m
	use tdse_m
	use output_m
	use rk4_m
	use coeffs_m
	use PtrArray_m
	implicit none
	
	private
	
	type :: solve_tdse_t
		type(data_t) :: dat
		type(tdse_t) :: tdse
		type(rk4_t)  :: integ
		type(output_t) :: output 
		type(coeffs_t) :: coeffs
		
		contains
		procedure :: pre => pre_solve
		procedure :: run => prop_tdse
		procedure :: reset => reset_solve
		procedure :: destroy => destroy_solve_tdse
	end type
	
	public :: solve_tdse_t
	
contains
!*******

subroutine pre_solve(this, fln, datdir)
	class(solve_tdse_t), intent(inout) :: this
	character(len=*),    intent(in)    :: fln, datdir
	
	integer(8) :: t1, t2, t_second
	
	call destroy_solve_tdse(this)

	write(*,'(A)')     "====================== PREPARING DATA ======================="
	write(*,'(A)')     ""
	
	call SYSTEM_CLOCK(t1)
	call this%dat%pre(fln,datdir)
	call SYSTEM_CLOCK(t2, t_second)
	
	write(*,'(A)')     repeat("=", 63) !"==============================================================="
	write(*,'(A,F,A)') 'Time preparing data consumption = ', dble(t2-t1)/t_second, ' s'
	write(*,'(A)')     ''

	call this%tdse%pre(fln,this%dat) 

    call this%integ%pre(this%dat%nbas)

    call this%coeffs%pre(fln, this%dat)
    
    call this%coeffs%setup_state()
    
	call this%output%pre(fln, this%dat, this%tdse)

end subroutine	
!-------------

subroutine reset_solve(this)
	class(solve_tdse_t), intent(in) :: this
	
	call this%coeffs%setup_state()
	
end subroutine

subroutine destroy_solve_tdse(this)
	class(solve_tdse_t), intent(inout) :: this

	call this%output%destroy()
	call this%coeffs%destroy()
	call this%integ%destroy()
	call this%tdse%destroy()
	call this%dat%destroy()
 
end subroutine
!-------------

subroutine prop_tdse(this)
	class(solve_tdse_t), intent(inout) :: this
	
	real(wp) :: TimeStep, t
	integer  :: NtMax
	
	complex(wp), allocatable :: Citmp(:)  
	complex(wp) :: tmp
	
	integer(8) :: t1, t2
	integer(8) :: t_second
	
	integer(8) :: t1_, t2_
	real(8)    :: t_elapse_solv, t_elapse_out
	
	integer :: mi, ni, nj, it
		
	Ntmax    = this%tdse%laser%Ntmax
	TimeStep = this%tdse%laser%TimeStep
	
	if( this%dat%IsMaskFunc ) then
	  allocate( Citmp(maxval(this%dat%nbasm)) )
	endif
	
	write(*,'(A)')         "====================== PROPAGATING TDSE ======================="
	write(*,'(A,F10.6,A)') 'dt = ', TimeStep, ' a.u.'

	t_elapse_solv = 0.0_8
	t_elapse_out  = 0.0_8
	
	call this%output%write( 0.0_wp, this%coeffs)

	call SYSTEM_CLOCK(t1)
		
	do it = 0, Ntmax-1
	    
	    if( mod(it+1, Ntmax/10) == 0 ) then  
 		  write(*,'(A,I0,A,I0,A, $)') char(13)//"t  =   ", it+1, "/", Ntmax, ' Ntmax' 
		endif 
		
		t = it * TimeStep
		
		call SYSTEM_CLOCK(t1_)
		call this%integ%update( this%tdse, this%coeffs%Ci, t, t+TimeStep)
		call SYSTEM_CLOCK(t2_, t_second)
		
		t_elapse_solv = t_elapse_solv + dble(t2_ - t1_) / t_second
		
		if (this%dat%isMaskFunc) then
!~ 			this%coeffs%Ci2d = matmul (this%dat%maskfunc, this%coeffs%Ci2d)
		  
		  !$OMP	PARALLEL PRIVATE (NI, NJ, TMP)
		  do mi = this%dat%mmin, this%dat%mmax  
		    		  
		    !this%coeffs%Cnm(mi)%pn = matmul(this%dat%maskfunc(mi,mi)%pn, this%coeffs%Cnm(mi)%pn)  
		    
		    !$OMP DO 
		    do ni = 1, this%dat%nbasm(mi)
		      Citmp(ni) = this%coeffs%Cnm(mi)%pn(ni) 
		    enddo 
		    !$OMP ENDDO
		    
		    !$OMP DO
		    do ni = 1, this%dat%nbasm(mi) 
		  
		      tmp = 0.0_wp  
		      do nj = 1, this%dat%nbasm(mi) 
		        tmp = tmp + this%dat%maskfunc(mi,mi)%pn(nj,ni) * Citmp(nj)
		      enddo
		    
		      this%coeffs%Cnm(mi)%pn(ni) = tmp
		    enddo
		    !$OMP END DO
		  enddo
		  !$OMP END PARALLEL
		endif
		
!~ 		  do ni = 1, this%dat%nbas 
!~ 		    Citmp(ni) = this%coeffs%Ci(ni)
!~ 		  enddo 
!~ 		  !$OMP END DO
		  
!~ 		  !$OMP DO 
!~ 		  do ni = 1, this%dat%nbas
!~ 		    tmp = 0.0_wp
!~ 			do nj = 1, this%dat%nbas 
!~ 			    tmp = tmp + this%dat%maskfunc(ni, nj) * Citmp(nj)
!~ 			enddo
!~ 			this%coeffs%Ci(ni)  = tmp 
!~ 		  enddo
!~ 		  !$OMP END DO
!~ 		  !$OMP END PARALLEL 
		
		call SYSTEM_CLOCK(t1_)
		call this%output%write(t+TimeStep, this%coeffs) 
		call SYSTEM_CLOCK(t2_)
		t_elapse_out = t_elapse_out + dble(t2_ - t1_) / t_second
		
	enddo
	call SYSTEM_CLOCK(t2, t_second)
	
	write(*,'(A)')     ""
	write(*,'(A)')     "==============================================================="
	write(*,'(A,F,A)') 'Total time consumption         = ', dble(t2-t1)/t_second, ' s'
	write(*,'(A,F,A)') 'Time solving consumption       = ', t_elapse_solv, ' s'
	write(*,'(A,F,A)') 'Time writing out consumption   = ', t_elapse_out , ' s'
	write(*,'(A)')     ""
	
!---
	if( allocated(Citmp) )  deallocate( Citmp )
	
end subroutine
!-------------

end module
