module coeffs_m
	use type_vars
	use data_m
	use PtrArray_m	
	implicit none
	
    private 	

	type :: coeffs_t
	    type(data_t), pointer :: dat =>null() 
	
		complex(wp), pointer :: Ci(:)
		complex(wp), pointer :: Ci2d(:,:)
	    type(CPtrArray1d), pointer :: Cnm(:)
	    
	    integer :: mini, nini

	    contains
	    procedure :: pre => pre_coeffs
	    procedure :: setup_state => setup_ini_state
	    procedure :: destroy => destroy_coeffs
	end type
	
	public :: coeffs_t 
	
contains
!*******

subroutine pre_coeffs(this, fln, dat)
	class(coeffs_t),  intent(inout) :: this
	character(len=*), intent(in)    :: fln
	type(data_t), target, intent(inout) :: dat 

	call destroy_coeffs(this) 

    this%dat => dat  

    allocate ( this%Ci(this%dat%nbas) )
	this%Ci2d(1:this%dat%nbas,1:1) => this%Ci(:)
	
	allocate( this%Cnm(this%dat%mmin : this%dat%mmax) )
	call PtrArray_associate (this%Cnm, this%Ci, this%dat%nbasm)

!---
	
	call read_gs_para(this, fln)
    call prt_gs_para(this)
    	
end subroutine
!-------------

subroutine setup_ini_state(this)
	class(coeffs_t), intent(in) :: this
	
    this%Ci = 0.0_wp

    if( this%dat%mmin > this%mini &
   .or. this%dat%mmax < this%mini &
   .or. this%dat%nbasm(this%mini) < this%nini ) &
		stop 'set_init_population Error: mini & nini out of range [mmin, mmax]!' 
    
    this%Cnm(this%mini)%pn(this%nini) = 1.0_wp

end subroutine
!-------------

subroutine read_gs_para(this, fln)
    class(coeffs_t),  intent(inout) :: this
    character(len=*), intent(in)    :: fln
    
	integer :: mini, nini
	integer :: idf
	
	mini = this%dat%mmin
	nini = 1
	
	namelist /gs_info/ mini, nini
	open(newunit = idf, file = fln, status = 'old', action = 'read')
		read(unit = idf, nml = gs_info)
	close(unit = idf)

	this%mini = mini
	this%nini = nini
	
end subroutine
!-------------

subroutine prt_gs_para(this)
	class(coeffs_t), intent(in) :: this
	
	write(*,'(A)')    "=== GROUND STATE INFO ==="
	write(*,'(A)')    ''
	write(*,'(A,I0)') 'nini = ', this%nini
	write(*,'(A,I0)') 'mini = ', this%mini
	write(*,'(A)')    ''
	
end subroutine
!-------------

subroutine destroy_coeffs(this)
	class(coeffs_t), intent(inout) :: this
	
	if( allocated(this%Ci) ) deallocate( this%Ci )
	
	this%Ci2d => null()
	
	if( allocated(this%Cnm) ) then 
	    call destroy_PtrArray(this%Cnm) 
		deallocate (this%Cnm)
	endif 

	this%dat => null()

end subroutine
!-------------
	
end module
