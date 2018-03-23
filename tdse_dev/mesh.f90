module mesh_m
	use type_vars
	implicit none
	
    private 

    type :: mesh_t 
      integer :: ikind
      integer :: Nrad, Neta
      real(wp) :: Rmax
    
      real(wp), allocatable :: rad(:), wrad(:)
      real(wp), allocatable :: eta(:), sneta(:), weta(:)
      real(wp), allocatable :: wre(:,:)
    
    contains
      procedure :: read_info => read_mesh_info
      procedure :: print => prt_mesh
      procedure :: create => create_mesh 
      procedure :: write => wrt_mesh
      procedure :: read => read_mesh
      procedure :: destroy => destroy_mesh
      procedure :: setup => setup_mesh 
    end type 

   public :: mesh_t 

contains
!*******

subroutine create_mesh(this)
  class(mesh_t), intent(inout) :: this

  call destroy_mesh(this)
  
  allocate( this%rad(this%nrad) )
  allocate( this%eta(this%neta), this%sneta(this%Neta) )
  allocate( this%wrad(this%Nrad), this%weta(this%neta), this%wre(this%Nrad,this%Neta) )

end subroutine 
!-------------

subroutine destroy_mesh(this)
  class(mesh_t), intent(inout) :: this
  
  if( allocated (this%rad) )   deallocate (this%rad)
  if( allocated (this%eta) )   deallocate (this%eta)
  if( allocated (this%sneta) ) deallocate (this%sneta)  
  if( allocated (this%wrad) )  deallocate (this%wrad)
  if( allocated (this%weta) )  deallocate (this%weta)
  if( allocated (this%wre) )   deallocate (this%wre)  

end subroutine 
!

subroutine read_mesh_info(this, inp_fln)
    class(mesh_t), intent(inout) :: this 
	character(len=*), intent(in) :: inp_fln
	
	integer  :: ikind, Nrad, Neta
	real(wp) :: Rmax
	
	integer :: idf
	
	namelist /mesh_info/ ikind, Nrad, Neta, Rmax
	
	open(newunit = idf, file = inp_fln, status = 'old', action = 'read')
		read(idf, nml = mesh_info)
	close(unit = idf)

    this%ikind = ikind
    this%Nrad  = Nrad
    this%Neta  = Neta
    this%Rmax  = Rmax

end subroutine
!-------------

subroutine prt_mesh(this)
	class(mesh_t), intent(in) :: this 

	write(*,'(a)')    "=== MESH INFO ==="
    write(*,'(a,g0)') 'ikind = ', this%ikind
    write(*,'(a,g0)') 'NRAD  = ', this%Nrad
    write(*,'(a,g0)') 'NETA  = ', this%Neta
    write(*,'(a,g0)') 'RMAX  = ', this%Rmax    
    write(*,*)
    
end subroutine
!------------- 

subroutine setup_mesh(this) 
	class(mesh_t), intent(inout) :: this

	integer :: ir, ie

	if (this%ikind == 1) then
		call cgqf ( this%Nrad, this%ikind, 0.0_wp, 0.0_wp, 0.0_wp, this%Rmax, this%rad, this%wrad )
	    call cgqf ( this%Neta, this%ikind, 0.0_wp, 0.0_wp, -1.0_wp, 1.0_wp, this%eta, this%weta )
	elseif (this%ikind == 2) then
		write(*,*) 'Experimental feature (under developed), please choose 1!'
	else
		write(*,*) 'Please choose ikind = 1 or 2!'
	endif
	
!-- calc weight	
	forall (ir = 1:this%Nrad, ie = 1:this%Neta)
		this%wre(ir,ie) = this%rad(ir)**2* this%wrad(ir) * this%weta(ie)
	end forall
	
end subroutine
!-------------

subroutine wrt_mesh(this,datdir)
    class(mesh_t), intent(inout) :: this
    character(len=*), intent(in) :: datdir
    
    integer :: idf
	integer :: ir, ie
		
	open(newunit = idf, file = trim(datdir)// '/' // 'Rinfo.in', &
						status = 'unknown', action = 'write')
		write(idf, '(2E16.8)') (this%rad(ir), this%wrad(ir), ir = 1, this%Nrad)
	close(unit = idf)

	open(newunit = idf, file = trim(datdir) // '/' // 'Einfo.in', &
						status = 'unknown', action = 'write')
		write(idf, '(3E16.8)') (this%eta(ie), sqrt(1.0_wp - this%eta(ie)**2), &
								this%weta(ie), ie = 1, this%Neta)
	close(unit = idf)
	
end subroutine
!-------------

subroutine read_mesh(this, datdir)
    class(mesh_t),    intent(inout) :: this
    character(len=*), intent(in)    :: datdir
	
	integer :: idf
	integer :: ir, ie
	
	open(newunit = idf, file = trim(datdir) // '/' // 'Rinfo.in', status = 'old', action = 'read')
		read(idf,*) (this%rad(ir), this%wrad(ir), ir = 1, this%Nrad)
	close(unit = idf)

	open(newunit = idf, file = trim(datdir) // '/' // 'Einfo.in', status = 'old', action = 'read')
		read(idf,*) (this%eta(ie), this%sneta(ie), this%weta(ie), ie = 1, this%Neta)
	close(unit = idf)
	
	forall (ir = 1:this%Nrad, ie = 1:this%Neta)
 		this%wre(ir,ie) = this%rad(ir)**2 * this%wrad(ir) * this%weta(ie) !!! for newest TISE 14 Oct
!~ 		this%wre(ir,ie) = this%wrad(ir) * this%weta(ie)
	end forall
	
end subroutine
!-------------

end module
