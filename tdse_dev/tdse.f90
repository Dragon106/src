module tdse_m
	use type_vars
	use consts
	use laser_m
	use data_m
	use PtrArray_m
	implicit none
	
	private 
	
	type :: tdse_t
	    type(elaser_t) :: laser 
        type(data_t), pointer :: dat => null()   
        
        complex(wp), pointer :: Er(:,:)
        type(CPtrArray2d), pointer :: dipy(:,:), dipz(:,:)
        
!~         complex(wp), pointer :: dipy(:,:), dipz(:,:)
        real(wp), pointer  :: Ei(:,:) 
        
        type(CPtrArray1d), pointer :: dCidtm(:), Cim(:)
        type(RPtrArray1d), pointer :: Eim(:)
        
        complex(wp), pointer :: zdipy(:,:), zdipz(:,:) 
        complex(wp), pointer :: zEi(:,:) 
                
        contains
        procedure :: pre       => pre_tdse
        procedure :: derivfunc => derivfunc
        procedure :: destroy   => destroy_tdse
	end type
	
	public :: tdse_t
	
contains
!*******

subroutine pre_tdse(this, fln, dat)
	class(tdse_t),    intent(inout) :: this
	character(len=*), intent(in)    :: fln 
	type(data_t),     target        :: dat
	
	call destroy_tdse(this) 
    call this%laser%read_params(fln)
    
    this%dat => dat 
	
!~ 	allocate( this%Er(this%dat%nbas, this%dat%nbas) )
	allocate (this%dCidtm (this%dat%mmin : this%dat%mmax) )

	allocate (this%Cim (this%dat%mmin : this%dat%mmax) )

	allocate ( this%Eim(this%dat%mmin:this%dat%mmax) )
	
	this%Ei(1:this%dat%nbas,1:1) => this%dat%Ei(1:this%dat%nbas)
	call PtrArray_associate (this%Eim, this%dat%Ei, this%dat%nbasm )

	if (this%dat%isDynCore) then
	this%dipy => this%dat%dyVp
	this%dipz => this%dat%dzVp
	else
	this%dipy => this%dat%dipy
	this%dipz => this%dat%dipz
	endif

!~ 	allocate( this%zEi(this%dat%nbas,1) )
!~ 	allocate( this%zdipy(this%dat%nbas, this%dat%nbas) )
!~ 	allocate( this%zdipz(this%dat%nbas, this%dat%nbas) )
	
!~ 	this%zEi = -im*this%Ei
!~ 	this%zdipy = -im*this%dipy
!~ 	this%zdipz = -im*this%dipz
	
 return 

	open(1, file = 'energy.dat')
		write(1,'(E16.8)') this%Ei(:,:)
	close(1)

!~ 	write(*,*) maxval( abs(this%dipz - conjg(transpose(this%dipz)) ) ), & 
!~ 			   maxval( abs(this%dipy - conjg(transpose(this%dipy)) ) ) 
stop
end subroutine
!-------------

subroutine destroy_tdse(this)
	class(tdse_t), intent(inout) :: this

	call this%laser%destroy() 
	
	if (allocated (this%Er) ) deallocate( this%Er )
	this%Ei => null()
	this%dipy => null()
	this%dipz => null()

	if ( allocated (this%dCidtm) ) then
		call destroy_PtrArray ( this%dCidtm )
		deallocate ( this%dCidtm )
	endif	
	
	if ( allocated (this%Cim) ) then
		call destroy_PtrArray ( this%Cim )
		deallocate ( this%Cim )
	endif	

	if ( allocated (this%Eim) ) then
		call destroy_PtrArray( this%Eim )
		deallocate (this%Eim)
	endif

	this%dat => null() 
	
end subroutine 
!-------------

subroutine derivfunc(this, dCidt, t, Ci) 
	class(tdse_t), intent(inout) :: this
	real(wp),      intent(in)    :: t
	complex(wp),   intent(inout) :: Ci(this%dat%nbas)
	complex(wp),   intent(out)   :: dCidt(this%dat%nbas)
	
    real(wp) :: Ey, Ez 	
	integer  :: mi, i, j
	complex(wp) :: tmpy, tmpz
	
	call this%laser%elaser(t, Ey, Ez) 
    
    call PtrArray_associate(this%dCidtm, dCidt, this%dat%nbasm)
    call PtrArray_associate(this%Cim, Ci, this%dat%nbasm)
    
  !$OMP PARALLEL PRIVATE(MI,I,J, TMPY,TMPZ) 
!~     !$OMP DO COLLAPSE(2)
!~     do j = 1, this%dat%nbas
!~     do i = 1, this%dat%nbas 
!~ 		this%Er(J,I) = Ey*this%dipy(i,j) + Ez*this%dipz(i,j)
!~ 	enddo
!~ 	enddo
!~     !$OMP END DO 
    
  
  do mi = this%dat%mmin, this%dat%mmax
  
!	dCidtm(mi)%pn = dCidtm(mi)%pn + matmul( this%dipz(mi,mi)%pn, Cim(mi)%pn ) 
	
	!$OMP DO 
	do i = 1, this%dat%nbasm(mi)

		tmpz = 0.0_wp 
		do j = 1, this%dat%nbasm(mi)
			tmpz = tmpz + ( this%dipz(mi,mi)%pn(j,i) ) * this%Cim(mi)%pn(j) 
		enddo		
		
	    tmpy = 0.0_wp
	    if( mi+1 <= this%dat%mmax ) then  
	    do j = 1, this%dat%nbasm(mi+1)
			tmpy = tmpy + this%dipy(mi,mi+1)%pn(j,i) * this%Cim(mi+1)%pn(j)
		enddo	
		endif 
		if (mi-1 >= this%dat%mmin ) then
		do j = 1, this%dat%nbasm(mi-1)
!~ 			tmpy = tmpy + this%dipy(mi,mi-1)%pn(i,j) * Cim(mi-1)%pn(j)
			tmpy = tmpy + conjg(this%dipy(mi-1,mi)%pn(i,j)) * this%Cim(mi-1)%pn(j)
		enddo
		endif
		
!		dCidtm(mi)%pn(i) = dCidtm(mi)%pn(i) + (tmpy_p1 + tmpy_m1) * Ey

		this%dCidtm(mi)%pn(i) = -im * ( this%Eim(mi)%pn(i)*this%Cim(mi)%pn(i) + Ez*tmpz + Ey*tmpy ) 
	enddo
    !$OMP END DO NOWAIT
  enddo
  !$OMP END PARALLEL 

  call destroy_PtrArray ( this%dCidtm )
  call destroy_PtrArray ( this%Cim )

!	do i = 1, this%dat%nbas
!	    tmp = 0.0_wp
!	    do j = 1, this%dat%nbas
!~ 	      tmp = tmp + this%Er(J,I)*Ci(j)
!	      tmp = tmp + ( Ey*this%dipy(i,j) + Ez*this%dipz(i,j) ) * Ci(j)
!	    enddo

!		dCidt(i) = -im * ( this%dat%Ei(i)*Ci(i) +  tmp )
!	enddo
  
!return 

!~ 	write(*,*) maxval( abs(this%Er - conjg(transpose(this%Er)) ) )  
!	write(*,*) maxval( abs( this%Er ))  

!    dCidt = this%zEi*Ci + matmul(this%Er,Ci)
 
end subroutine
!-------------

!~ subroutine derivfunc_(this, dCidt, t, Ci) 
!~ 	class(tdse_t), intent(inout) :: this
!~ 	real(wp), intent(in) :: t
!~ 	complex(wp), intent(in) :: Ci(this%dat%nbas,1)
!~ 	complex(wp), intent(out) :: dCidt(this%dat%nbas,1)
	
!~     real(wp) :: Ey, Ez 	
	
!~ 	call this%laser%elaser(t, Ey, Ez) 
    
!~ 	this%Er = Ey*this%dipy + Ez*this%dipz 

!~     dCidt = -im * ( this%Ei*Ci + matmul(this%Er,Ci) )

!~ return 

!~ 	write(*,*) maxval( abs(this%Er - conjg(transpose(this%Er)) ) )  
!	write(*,*) maxval( abs( this%Er ))  

!    dCidt = this%zEi*Ci + matmul(this%Er,Ci)

!~ end subroutine
!-------------

end module
 

