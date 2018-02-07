module calc_rate_m
	use type_vars_m
	use derivfx_m
	implicit none
	
	private
	
	type :: ioniz_t
	  character(len = 100) :: datadir, dataname, outname
	  
	  integer  :: thetamin, thetamax, thestep
	  
	  integer  :: Nt
	  real(wp) :: timestep 
	  real(wp), pointer :: time(:)
	  
	  real(wp), pointer :: prob(:)
	  real(wp), pointer :: rate(:)
	
	  contains
!~ 		procedure :: read_para     => read_para
!~ 		procedure :: read_ntmax    => read_ntmax
!~ 		procedure :: get_timestep  => calc_timestep
		procedure :: pre_calc => pre_calc
		
!~ 		procedure :: generate_vars => alloc_vars  	
!~ 		procedure :: read_data     => read_prob_data
!~ 		procedure :: calc_ion_rate => calc_ion_rate
!~ 		procedure :: write         => write_ion_rate
		procedure :: run => calc

	end type
	
	public :: ioniz_t
	
contains
!*******
subroutine pre_calc(this)
	class(ioniz_t), intent(inout) :: this
	
	call read_para(this)
	call read_ntmax(this)
	call calc_timestep(this)
	
end subroutine
!-------------

subroutine calc(this)
	class(ioniz_t), intent(inout) :: this
	
	call alloc_vars(this)
	call read_prob_data(this)
	call calc_ion_rate(this)
	call write_ion_rate(this)
	
end subroutine
!-------------

subroutine read_para(this)
	class(ioniz_t), intent(inout) :: this
	
	character(len = 100) :: datadir
	integer :: thetamin, thetamax, thestep
	integer :: idf
	
	namelist /inout_info/ datadir
	namelist /orient_info/ thetamin, thetamax, thestep
	
	open(newunit = idf, file = 'para_rate.inp', status = 'old', action = 'read')
		read(idf, nml = inout_info)
		read(idf, nml = orient_info)
	close(idf)
	
	this%datadir  = datadir
	this%thetamin = thetamin
	this%thetamax = thetamax
	this%thestep  = thestep
	
!~ 	write(*,*) this%datadir
!~ 	write(*,*) this%thetamin, this%thetamax, this%thestep

end subroutine
!-------------

subroutine read_ntmax(this)
	class(ioniz_t), intent(inout) :: this
	
	character(len=100) :: namein
	integer :: nc
	integer :: idf, ios
	
	write(namein, '(A,I0,A)') trim(this%datadir) // '/', this%thetamin, '/' // 'prob.dat'

	nc = 0
	open(newunit = idf, file = namein, status = 'old', action = 'read')
	do
		read(idf,*,iostat = ios)
		if (ios /= 0) exit
		nc = nc + 1
	enddo
	close(unit = idf)
	this%Nt = nc
	
	if( this%Nt < 2 ) then
		write(*,*) 'Data not enough, propagate more!'
		stop
	endif
!~ 	write(*,*) this%Nt

end subroutine
!-------------

subroutine calc_timestep(this)
	class(ioniz_t), intent(inout) :: this
	
	character(len=100) :: namein
	real(wp) :: t(2)
	integer :: idf
	integer :: i
	
	write(namein, '(A,I0,A)') trim(this%datadir) // '/', this%thetamin, '/' // 'prob.dat'
	
	open(newunit = idf, file = namein, status = 'old', action = 'read')
	do i = 1, 2
		read(idf,*) t(i)
	enddo
	close(idf)
	
	this%timestep = t(2) - t(1)
	write(*,*) this%timestep
	
end subroutine
!-------------

subroutine alloc_vars(this)
	class(ioniz_t), intent(inout) :: this
	
	if( allocated(this%time) ) deallocate( this%time ) 
	if( allocated(this%prob) ) deallocate( this%prob ) 
	if( allocated(this%rate) ) deallocate( this%rate )

	allocate( this%time(this%Nt) )
	allocate( this%prob(this%Nt) )
	allocate( this%rate(this%Nt) )
	
end subroutine
!-------------

subroutine read_prob_data(this)
	class(ioniz_t), intent(inout) :: this
	
	real(wp) :: dummy
	integer :: idf
	integer :: i
	
	open(newunit = idf, file = this%dataname, status = 'old', action = 'read')
	do i = 1, this%Nt
		read(idf,*) this%time(i), this%prob(i)
	enddo	
	close(unit = idf) 
	
end subroutine
!-------------

subroutine calc_ion_rate(this)
	class(ioniz_t), intent(inout) :: this
	
!~ 	call deriv( this%timestep, this%prob, this%rate )
	call deriv( this%timestep, log(1.0_wp - this%prob), this%rate )
!~ 	call deriv( this%timestep,    (1.0_wp - this%prob), this%rate )
end subroutine
!-------------

subroutine write_ion_rate(this)
	class(ioniz_t), intent(inout) :: this
	
	integer :: idf
	integer :: i
	
	open(newunit = idf, file = this%outname, status = 'unknown', action = 'write')
!~ 		write(idf,'(2E16.8)') ( this%time(i), this%rate(i), i = 1, this%Nt )
		write(idf,'(2E16.8)') ( this%time(i), -this%rate(i), i = 1, this%Nt )
	close(unit = idf)
	
end subroutine
!-------------

end module
!*********
