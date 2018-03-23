module data_m
	use type_vars 
	use consts
	use mesh_m
	use PtrArray_m
	USE OMP_LIB
	implicit none
	
	private
	
	type :: data_t 
	  	integer, allocatable :: nbasm(:)
	  	integer  :: mmin, mmax
	  	real(wp) :: emax
		
		logical :: isRmOrbital
		integer :: NrmOrb
		integer, pointer :: rmOrbIdx(:,:)
		
		integer :: nbas
		real(wp), pointer :: Ei(:)
		type(RPtrArray1d), pointer :: Enm(:) 
	    
		type(CPtrArray2d), pointer :: dipy(:,:), dipz(:,:)
		
!		complex(wp), pointer :: dipy(:,:), dipz(:,:)
		
		logical  :: isDynCore = .false.
		real(wp) :: alphaxy, alphazz
		
		type(CPtrArray2d), pointer :: dzVp(:,:), dyVp(:,:)
		
!		complex(wp), pointer :: dzVp(:,:), dyVp(:,:)
		complex(wp), pointer :: Vabs(:,:)
		
		type(mesh_t) :: mesh
		
		logical  :: isMaskFunc
		real(wp) :: z0
		
		type(CPtrArray2d), pointer :: maskfunc(:,:)
!		real(wp), pointer :: maskfunc(:,:)
		
		contains
		
		procedure :: pre     => pre_data
		procedure :: destroy => destroy_data
	end type 
	
	public :: data_t

contains 
!*******

subroutine pre_data(this, fln, datdir) 
	class(data_t),    intent(inout) :: this
	character(len=*), intent(in)    :: fln, datdir

    call destroy_data(this) 

!---
    call read_inp(this, fln) 
    
	call get_nbasm(this, datdir)
	
	call prt_basis_info(this)
	
	call read_energy(this, datdir)
	
	call read_mesh(this, datdir)
	
	call rm_orbital(this, fln)

    call calc_dipole(this, fln, datdir)

end subroutine 
!-------------

subroutine destroy_data(this)
	class(data_t), intent(inout) :: this
	
	this%emax = 0.0d0
	this%mmin = 0
	this%mmax = 0
	this%nbas = 0 
	
	if( allocated(this%nbasm) ) deallocate( this%nbasm )
	if( allocated(this%Ei) )    deallocate( this%Ei )

!	if( allocated(this%dipy) )  deallocate( this%dipy )	
!	if( allocated(this%dipz) )  deallocate( this%dipz )	
	
	call this%mesh%destroy()
	
end subroutine
!-------------

subroutine read_inp(this, inp_fln)
	class(data_t),    intent(inout) :: this
	character(len=*), intent(in)    :: inp_fln
	
	integer  :: mmin, mmax
	real(wp) :: emax
	
	integer  :: idf
	
	namelist /basis_info/ mmin, mmax, emax
	
	open (newunit = idf, file = inp_fln, status = 'old', action = 'read')
		read(unit = idf, nml = basis_info)
	close (unit = idf)

	this%mmin = mmin
	this%mmax = mmax
	this%emax = emax

end subroutine
!-------------

subroutine prt_basis_info(this)
	class(data_t), intent(in) :: this
	
	write(*,'(A)')         "=== BASIS INFO ==="
	write(*,'(A, F8.3,A)') 'Emax = ', this%emax, ' a.u. '
	write(*,'(A, I0)')     'NBAS = ', this%nbas
	write(*,'(A, I0)')     'mmin = ', this%mmin
	write(*,'(A, I0)')     'mmax = ', this%mmax
	write(*,'(A)')         ""
		
end subroutine
!-------------

subroutine rm_orbital(this, fln)
	class(data_t),    intent(inout) :: this
	character(len=*), intent(in)    :: fln
	
	logical :: isRmOrbital
	character(len=100) :: diagram_file
	integer :: idf, ios, nrow, NrmOrb
	integer :: i
	
	isRmOrbital = .false.
	
	namelist /RmOrb_info/ isRmOrbital, diagram_file
	
	open (newunit = idf, file = fln, status = 'old', action = 'read')
		read(unit = idf, nml = RmOrb_info)
	close(unit = idf)

	this%isRmOrbital = isRmOrbital
	
	if ( .not. this%isRmOrbital ) then 
	
	  write(*,'(A)') "=== REMOVE ORBITALS INFO ==="
	  write(*,'(A)') 'Propagating full orbitals.'
	  write(*,'(A)') ""
	  return 
	
	endif 
!--	
	open(newunit = idf, file = diagram_file, status = 'old', action = 'read')
		nrow = 0
		do
			read(idf, *, iostat = ios) i, i
			if (ios /= 0) exit
			nrow = nrow + 1
		enddo
		
	close(unit = idf)
!--	
	this%NrmOrb = nrow
	allocate( this%rmOrbIdx(this%NrmOrb,2) )
	
	open(newunit = idf, file = diagram_file, status = 'old', action = 'read')
		read(idf, *) (this%rmOrbIdx(i,:), i = 1, this%NrmOrb)
	close(unit = idf) 
!--	
	NrmOrb = this%NrmOrb!!!
	
	do i = 1, this%NrmOrb 
		if( this%rmOrbIdx(i,1) < this%mmin & 
	   .or. this%rmOrbIdx(i,1) > this%mmax ) then
			
!~ 			stop 'RmOrbital: mIndex out of range.'
			write(*,'(2(A, I0))') 'RmOrbital: mIndex out of range m = ', this%rmOrbIdx(i,1), ', n = ', this%rmOrbIdx(i,2)
			NrmOrb = NrmOrb - 1 !!!
			
	    !--- if (m > n or n < 1 ) stop program
	    elseif( this%rmOrbIdx(i,2) > this%nbasm(this%rmOrbIdx(i,1)) &
	       .or. this%rmOrbIdx(i,2) < 1 ) then 
	        
	        stop 'RmOrbital: nIndex out of range.'
	        
		endif 	
	enddo

!--	
	write(*,'(A)')            "=== REMOVE ORBITALS INFO ==="
	write(*,'(A,I0)')         'Number of orbitals removed = ', NrmOrb !this%NrmOrb !!!
!~ 	write(*,'(2I5,3x,F10.6)') (this%rmOrbIdx(i,:), &
!~ 							  this%Enm(this%rmOrbIdx(i,1)) % pn( this%rmOrbIdx(i,2) ), &
!~ 							  i = 1, this%NrmOrb) 
	do i = 1, this%NrmOrb
		if (this%rmOrbIdx(i,1) >= this%mmin & 
	  .and. this%rmOrbIdx(i,1) <= this%mmax ) then 
	    write(*,'(2I5,3x,F10.6)') this%rmOrbIdx(i,:), &
 							    this%Enm(this%rmOrbIdx(i,1)) % pn( this%rmOrbIdx(i,2) )
 		endif					    
	enddo
	write(*,'(A)')            "" 
	
end subroutine
!-------------

subroutine read_mesh(this, datdir)
	class(data_t),    intent(inout) :: this
	character(len=*), intent(in)    :: datdir
	
    call this%mesh%read_info( trim(datdir) // '/' // 'mesh.in' )
    call this%mesh%print()
    call this%mesh%create() 
    call this%mesh%read(datdir) 

end subroutine
!-------------

subroutine get_nbasm(this, datdir)
    class(data_t),    intent(inout) :: this
	character(len=*), intent(in)    :: datdir

	character(len = 100) :: namein
	real(wp) :: En
	integer :: ncount	
	integer :: ios, idf, idum
	integer :: m
	
	allocate (this%nbasm (this%mmin : this%mmax) )
	
	do m = this%mmin, this%mmax
		ncount = 0

		write(namein, '(A,I0,A)') trim(datdir) // "/" // 'EM-', abs(M), '.in'
		open(newunit = idf, file = namein, status = 'old', action = 'read')
		
		do 
			read(idf, *, iostat=ios) idum, En
			
			if (ios /= 0) exit
			
			if (En <= this%Emax) then
				ncount = ncount + 1
			else
				exit
			endif
			
		enddo
		
		this%nbasm(m) = ncount
		
		close (unit = idf)
		
	enddo

	this%nbas = sum( this%nbasm )
	
end subroutine
!-------------

subroutine read_energy(this, datdir)
	class(data_t),    intent(inout) :: this
	character(len=*), intent(in)    :: datdir
	
	character(len=100) :: namein
	integer :: offset
	integer :: idf, idum
	integer :: m, n
	
	allocate( this%Ei(this%nbas) )
	
	do m = this%mmin, this%mmax
		write(namein, '(A,I0,A)') trim(datdir) // "/" // 'EM-', abs(M), '.in'
		open(newunit = idf, file = namein, status = 'old', action = 'read')
		
		offset = sum(this%nbasm(this%mmin:m-1))
		do n = 1, this%nbasm(m)
			read(idf, *) idum, this%Ei( offset + n )
		enddo
		
		close(unit = idf)

	enddo
	
	allocate( this%Enm(this%mmin:this%mmax) ) 
	call PtrArray_associate(this%Enm, this%Ei, this%nbasm )
	
end subroutine
!-------------

subroutine read_wfm(this, datdir, wfm, m)
    class(data_t),    intent(in)                 :: this
	character(len=*), intent(in)                 :: datdir
	integer,          intent(in)                 :: m
	real(wp),         intent(inout), allocatable :: wfm(:,:,:)

	character(len = 100) :: namein
	integer :: idf, ios
	integer :: i
	
	if( allocated(wfm) ) deallocate (wfm)		
	allocate (wfm(this%mesh%Nrad, this%mesh%Neta, this%nbasm(m) ) )
	
	write(namein, '(A,I0,A)') trim(datdir) // '/' // 'WF-', abs(m), '.in'

!-- reading	binary data file
	open(newunit = idf, file = namein, form ='unformatted', status = 'old', action = 'read')
		read(idf, iostat = ios) wfm 
	close(unit = idf)

	if ( this%isRmOrbital ) then
	
	  do i = 1, this%NrmOrb
		if( m == this%rmOrbIdx(i,1) ) then 
		  wfm(:,:,this%rmOrbIdx(i,2)) = 0.0_wp
		endif
	  enddo 
	
	endif
	
end subroutine
!-------------	

subroutine get_maskfunc_para(this, fln)
	class(data_t),    intent(inout) :: this
	character(len=*), intent(in)    :: fln
	
	logical  :: isMaskFunc
	real(wp) :: z0
	integer  :: idf
	
	isMaskFunc = .false.
	
	namelist /maskfunc_info/ isMaskFunc, z0

	open(newunit = idf, file = fln, status = 'old', action = 'read')
		read(unit = idf, nml = maskfunc_info)
	close(unit = idf)
	
	this%isMaskFunc = isMaskFunc
	this%z0 = z0

!--	
	if (z0 >= this%mesh%Rmax) stop 'Mask function error: z0 out of box range!'
!--	
	write(*,'(A)')            "=== MASK FUNCTION INFO ==="
	if (this%isMaskFunc) then
		write(*,'(A,F6.2,A)') 'Mask function starts at = ', this%z0, ' a.u..'
	else
		write(*,'(A)')        'Mask function is not used.'
	endif	
	write(*,'(A)')            ""
	
end subroutine
!-------------

function mask_function_sin(z, z0, a)
	real(wp), intent(in) :: z(:), z0, a
	real(wp) :: mask_function_sin(size(z))
	
	real(wp) :: zmax
	
	zmax = z0 + a
	
	!$OMP PARALLEL WORKSHARE 
	where ( z >= z0 .and. z <= zmax ) 
		mask_function_sin = sin( (zmax - z) * pi / (2.0_wp*a) )**(1/4.0_wp)
	else where
		mask_function_sin = 1.0_wp
	end where
	!$OMP END PARALLEL WORKSHARE 
		
end function
!-----------

function mask_function(z, z0, a)
	real(wp), intent(in) :: z(:), z0, a
	real(wp) :: mask_function(size(z))
	
	real(wp) :: zmax
	
	zmax = z0 + a
	
	!$OMP PARALLEL WORKSHARE 
	where ( z > z0 .and. z <= zmax ) 
		mask_function = cos( (z - z0) * pi / (2.0_wp*a) )**(1/8.0_wp)
	else where
		mask_function = 1.0_wp
	end where
	!$OMP END PARALLEL WORKSHARE 
		
end function
!-----------

subroutine calc_dipole(this, fln, datdir)
	use PtrArray_m
	
	class(data_t),     intent(inout) :: this
	character(len=*),  intent(in)    :: fln, datdir
	type(CPtrArray2d), allocatable   :: dipym(:,:), dipzm(:,:)
	type(RPtrArray2d), allocatable   :: maskfuncm(:,:)
  
    real(wp), allocatable :: wfmi(:,:,:), wfmj(:,:,:) 
    real(wp), allocatable :: rz(:,:), ry(:,:), maskfunc_re(:,:)
	
	real(wp), allocatable :: rzVp(:,:), ryVp(:,:)
	type(CPtrArray2d), allocatable :: dyVpm(:,:), dzVpm(:,:)
	
	real(wp) :: rxy, rzz

  	integer :: ir, ie, m, mi, mj, ni, nj

!---
  	write(*,'(A)') "=== DIPOLE ELEMENTS CALCULATING ==="
  	write(*,'(A)') 'Please take a cup of tea, it takes few seconds to hours...'
  	write(*,'(A)') ""

!-----------	
  	call get_maskfunc_para(this, fln)
  	
  	if (this%isMaskFunc) then
		
		allocate( maskfunc_re(this%mesh%nrad,this%mesh%neta) )
		
		!$OMP PARALLEL DO
		do ie = 1, this%mesh%neta
			maskfunc_re(:,ie) = mask_function(this%mesh%rad, this%z0, this%mesh%Rmax - this%z0) &
								* this%mesh%wre(:,ie)
		enddo
		!$OMP END PARALLEL DO		
		
		if (allocated (this%maskfunc) ) then 
			call destroy_PtrArray ( this%maskfunc )
			deallocate (this%maskfunc)
		endif 
!		allocate ( this%maskfunc(this%nbas, this%nbas) )
		
		allocate ( this%maskfunc(this%mmin : this%mmax, this%mmin : this%mmax) )
		
		do m = this%mmin, this%mmax
			allocate( this%maskfunc(m,m) %pn(1: this%nbasm(m), 1:this%nbasm(m) ) )
			this%maskfunc(m,m)%pn = 0.0_wp
		enddo
		

!		allocate ( maskfuncm(this%mmin: this%mmax, this%mmin : this%mmax) )		
!		call PtrArray_associate(maskfuncm, this%maskfunc, this%nbasm)
		
!		this%maskfunc = 0.0_wp
		
	endif
	
!-----------	
	
	call get_dyncore_para(this, fln)
	
	if (this%isDynCore) then
		allocate( rzVp(this%mesh%Nrad, this%mesh%Neta), ryVp(this%mesh%Nrad, this%mesh%Neta) )

		rxy = (this%alphaxy)**(1/3.0_wp)
		rzz = (this%alphazz)**(1/3.0_wp)
		
		!$OMP PARALLEL DO PRIVATE(IE, IR) COLLAPSE (2)
		do ie = 1, this%mesh%neta
		
			do ir = 1, this%mesh%nrad 
		
			   rzVp(ir,ie) = this%mesh%rad(ir) * this%mesh%eta(ie) &
							 * (1.0_wp - this%alphazz / this%mesh%rad(ir)**3) &
							 * this%mesh%wre(ir,ie)
			   ryVp(ir,ie) = this%mesh%rad(ir) * this%mesh%sneta(ie) &
							 * (1.0_wp - this%alphaxy / this%mesh%rad(ir)**3) &
							 * this%mesh%wre(ir,ie)
			   if ( ( (this%mesh%rad(ir) * this%mesh%sneta(ie) / rxy)**2 &
					+ (this%mesh%rad(ir) * this%mesh%eta(ie)   / rzz)**2) <= 1.0_wp) then
					rzVp(ir,ie) = 0.0_wp
					ryVp(ir,ie) = 0.0_wp
			   endif
			enddo
		enddo		
		!$OMP END PARALLEL DO
		
!		if (allocated (this%dyVp)) deallocate (this%dyVp)
!		if (allocated (this%dzVp)) deallocate (this%dzVp)

		if (allocated (this%dyVp)) then
			call destroy_PtrArray (this%dyVp)
			deallocate (this%dyVp)
		endif
		
		if (allocated (this%dzVp)) then
			call destroy_PtrArray (this%dzVp)
			deallocate (this%dzVp)
		endif
		
		allocate (this%dyVp(this%mmin:this%mmax,this%mmin:this%mmax))
		allocate (this%dzVp(this%mmin:this%mmax,this%mmin:this%mmax))
		
		do m = this%mmin,this%mmax
			allocate ( this%dzVp(m,m)% pn(1:this%nbasm(m), 1:this%nbasm(m) ) )
			this%dzVp(m,m)%pn = 0.0_wp
		enddo
		
		do m = this%mmin, this%mmax - 1
			allocate ( this%dyVp(m,m+1) % pn(1:this%nbasm(m+1), 1:this%nbasm(m) ) )
			this%dyVp(m,m+1)%pn = 0.0_wp
		enddo
		
!		allocate( this%dyVp(this%nbas, this%nbas), this%dzVp(this%nbas, this%nbas) )
!		allocate( dyVpm(this%mmin : this%mmax, this%mmin : this%mmax), &
!				  dzVpm(this%mmin : this%mmax, this%mmin : this%mmax) )

!		call PtrArray_associate (dyVpm, this%dyVp, this%nbasm)
!		call PtrArray_associate (dzVpm, this%dzVp, this%nbasm)
  
!		this%dyVp = 0.0_wp
!		this%dzVp = 0.0_wp
		
	endif 
	
!-----------------------------------------------------------------------	
    allocate( rz(this%mesh%Nrad, this%mesh%Neta), ry(this%mesh%Nrad, this%mesh%Neta) )

    !$OMP PARALLEL WORKSHARE
    forall (ir=1:this%mesh%Nrad, ie = 1:this%mesh%Neta)
      rz(ir, ie) = this%mesh%rad(ir) * this%mesh%eta(ie)   * this%mesh%wre(ir,ie)
      ry(ir, ie) = this%mesh%rad(ir) * this%mesh%sneta(ie) * this%mesh%wre(ir,ie)
    end forall
    !$OMP END PARALLEL WORKSHARE
!------------  
    
!    if (allocated (this%dipy)) deallocate (this%dipy)
!    if (allocated (this%dipz)) deallocate (this%dipz)
!    allocate( this%dipy(this%nbas, this%nbas), this%dipz(this%nbas, this%nbas) )
!    allocate( dipym(this%mmin : this%mmax, this%mmin : this%mmax), &
!    		  dipzm(this%mmin : this%mmax, this%mmin : this%mmax) )
  

!    call PtrArray_associate (dipym, this%dipy, this%nbasm)
!    call PtrArray_associate (dipzm, this%dipz, this%nbasm)
    
!    this%dipz = 0.0_wp 
!    this%dipy = 0.0_wp

	if ( allocated (this%dipy) ) then
		call destroy_PtrArray( this%dipy )
		deallocate ( this%dipy )
	endif
	if ( allocated (this%dipz) ) then
		call destroy_PtrArray( this%dipz )
		deallocate ( this%dipz )
	endif
	
	allocate (this%dipy(this%mmin:this%mmax, this%mmin:this%mmax))
	allocate (this%dipz(this%mmin:this%mmax, this%mmin:this%mmax))
	
	do m = this%mmin, this%mmax
		allocate ( this%dipz(m,m) % pn(1:this%nbasm(m), 1:this%nbasm(m) ) )
		this%dipz(m,m)%pn = 0.0_wp
	enddo
	
	do m = this%mmin, this%mmax - 1
		allocate ( this%dipy(m,m+1) % pn(1:this%nbasm(m+1), 1:this%nbasm(m) ) )
		this%dipy(m,m+1)%pn = 0.0_wp
	enddo
!----------------------  
!~ goto 1988
    do mi = this%mmin, this%mmax
    
!~       write(*,*) mi
    
	  call read_wfm(this, datdir, wfmi, mi) 
	  !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(NJ,NI) 
	  do nJ = 1, this%nbasm(mi); do nI = 1, this%nbasm(mI)    
!		dipzm(mi,mi)%pn(ni,nj) = sum(wfmi(:,:,ni)* rz(:,:) * wfmi(:,:,nj) ) 
		this%dipz(mi,mi)%pn(ni,nj) = sum( wfmi(:,:,ni) * rz(:,:) * wfmi(:,:,nj) ) 
	  enddo; enddo
	  !$OMP END PARALLEL DO 
	  
	  if (this%isDynCore) then
	  !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(NJ,NI) 
	    do nJ = 1, this%nbasm(mi); do ni = 1, this%nbasm(mi)
!		  dzVpm(mi,mi)%pn(ni,nj) = sum(wfmi(:,:,ni)* rzVp(:,:) * wfmi(:,:,nj) ) 	
		  this%dzVp(mi,mi)%pn(ni,nj) = sum(wfmi(:,:,ni)* rzVp(:,:) * wfmi(:,:,nj) ) 	
	    enddo; enddo
	  !$OMP END PARALLEL DO 
	  endif
	  
	  if( this%isMaskFunc ) then
	  !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(NJ,NI) 
		do nJ = 1, this%nbasm(mi); do nI = 1, this%nbasm(mi)
!			maskfuncm(mi,mi)%pn(ni,nj) = sum(wfmi(:,:,ni)* maskfunc_re(:,:) * wfmi(:,:,nj) ) 
			this%maskfunc(mi,mi)%pn(ni,nj) = sum(wfmi(:,:,ni)* maskfunc_re(:,:) * wfmi(:,:,nj) ) 
		enddo; enddo
	  !$OMP END PARALLEL DO 
	  endif
	  
	  do mj = mi, this%mmax
	  
	  if (mj - mi == -1) then
!~ 	    write(*,*) 'mi = ', mi, 'mj = ', mj 
	
		call read_wfm(this, datdir, wfmj, mj) 
	    !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(NJ,NI) 
		do nJ = 1, this%nbasm(mJ); do nI = 1, this%nbasm(mI)
		  this%dipy(mi,mj)%pn(nj,ni) = -0.5_wp*im* sum(wfmi(:,:,ni) * ry(:,:) * wfmj(:,:,nj) )
		enddo; enddo
	    !$OMP END PARALLEL DO 
		
		if (this%isDynCore) then
	      !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(NJ,NI) 
		  do nJ = 1, this%nbasm(mJ); do nI = 1, this%nbasm(mI)
		    this%dyVp(mi,mj)%pn(nj,ni) = -0.5_wp*im* sum(wfmi(:,:,ni) * ryVp(:,:) * wfmj(:,:,nj) )
		  enddo; enddo
		  !$OMP END PARALLEL DO 
		endif
	
	  elseif(mj - mi == +1) then
	
		call read_wfm(this, datdir, wfmj, mj) 
	
		!$OMP PARALLEL DO COLLAPSE(2) PRIVATE(NJ,NI) 
		do nJ = 1, this%nbasm(mJ); do nI = 1, this%nbasm(mI)
		  this%dipy(mi,mj)%pn(nj,ni) = +0.5_wp*im* sum(wfmi(:,:,ni) * ry(:,:) * wfmj(:,:,nj) )
		enddo; enddo
		!$OMP END PARALLEL DO 
		
		if (this%isDynCore) then
		  !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(NJ,NI) 
		  do nJ = 1, this%nbasm(mJ); do nI = 1, this%nbasm(mI)
		    this%dyVp(mi,mj)%pn(nj,ni) = +0.5_wp*im* sum(wfmi(:,:,ni) * ryVp(:,:) * wfmj(:,:,nj) )
		  enddo; enddo
		  !$OMP END PARALLEL DO 
		endif 
		
      endif 

	  enddo
    enddo
!~ 1988 continue 	
!	!$OMP PARALLEL DO PRIVATE(MI,MJ,NI,NJ) 
!	do mi = this%mmin, this%mmax - 1
!		do mj = mi+1, this%mmax
!		    do nj = 1, this%nbasm(mj) ; do ni = 1, this%nbasm(mi)
!			  dipym(mi,mj)%pn(ni,nj) = conjg(dipym(mj,mi)%pn(nj,ni))
!			enddo; enddo 
!		enddo
!	enddo
!	!$OMP END PARALLEL DO
!~  read(*,*) 
!	if (this%isDynCore) then
!		!$OMP PARALLEL DO PRIVATE(MI,MJ)
!	  	do mi = this%mmin, this%mmax - 1
!		  do mj = mi+1, this%mmax
!			do nj = 1, this%nbasm(mj); do ni = 1, this%nbasm(mi)
!			  dyVpm(mi,mj)%pn(ni,nj) = conjg(dyVpm(mj,mi)%pn(nj,ni)) 
!			enddo; enddo  
!		  enddo
!		enddo
!		!$OMP END PARALLEL DO
!	endif

	
	if (allocated(wfmi) ) deallocate (wfmi)	
	if (allocated(wfmj) ) deallocate (wfmj)
	
	deallocate( rz, ry )
  
end subroutine
!-------------

subroutine get_dyncore_para(this, fln)
	class(data_t),    intent(inout) :: this
	character(len=*), intent(in)    :: fln
	
	logical  :: isDynCore
	real(wp) :: alphaxy, alphazz

	integer :: idf
	
	namelist /dynamic_core/ isDynCore, alphaxy, alphazz
	
	isDynCore = .false.
	
	open(newunit = idf, file = fln, status = 'old', action = 'read')
		read(unit = idf, nml = dynamic_core)
	close(unit = idf)
	
	this%isDynCore = isDynCore
	this%alphaxy = alphaxy
	this%alphazz = alphazz

	write(*,'(A)')           "=== DYNAMIC CORE INFO ==="
	if (this%isDynCore) then
		write(*,'(A,F10.6)') 'ALPHAXY = ', this%alphaxy
		write(*,'(A,F10.6)') 'ALPHAZZ = ', this%alphazz
	else
		write(*,'(A)')       "Without including dynamic core polarization."
	end if
	write(*,'(A)')           ""
	
end subroutine
!-------------

end module

