module output_m
	use type_vars
	use data_m
	use tdse_m
	use coeffs_m
	implicit none
	
	private
	
	type :: output_t
	    type(data_t), pointer :: dat  => null()
	    type(tdse_t), pointer :: tdse => null() 
	    
		character(len = 100)  :: outdir = './output-tdse/'
		
		logical :: isOpened = .false.
		logical :: isout_laser, isout_dip, isout_norm, isout_prob, isout_energy, isout_coeff
		
		integer :: idf_laser, idf_dip, idf_norm, idf_prob, idf_energy, idf_coeff
		
		contains
		procedure :: pre     => pre_output
		procedure :: setdir  => set_outputdir
		procedure :: write   => write_output
		procedure :: destroy => destroy_output
	end type
	
	public :: output_t
	
contains 
!*******

subroutine pre_output(this, fln, dat, tdse) 
	class(output_t),  intent(inout) :: this
	character(len=*), intent(in)    :: fln
	type(data_t), target, intent(inout) :: dat 
	type(tdse_t), target, intent(inout) :: tdse
	
	call destroy_output(this) 
	
	this%dat  => dat 
	this%tdse => tdse 
	
	call read_para(this, fln)
	
end subroutine
!-------------

subroutine read_para(this, fln)
	class(output_t),  intent(inout) :: this
	character(len=*), intent(in)    :: fln
	integer :: idf
	logical :: isout_laser, isout_dip, isout_norm, isout_prob, isout_energy, isout_coeff
	
	isout_laser  = .true.
	isout_dip    = .true.
	isout_norm   = .false.
	isout_prob   = .false.
	isout_energy = .false.
	isout_coeff  = .false.
	
	namelist /output_info/ isout_laser, isout_dip, isout_norm, isout_prob, &
						   isout_energy, isout_coeff
	open (newunit = idf, file = fln, status = 'old', action = 'read')
		read(unit = idf, nml  = output_info)
	close(unit = idf)

	this%isout_laser  = isout_laser
	this%isout_dip    = isout_dip
	this%isout_norm   = isout_norm
	this%isout_prob   = isout_prob
	this%isout_energy = isout_energy
	this%isout_coeff  = isout_coeff
	
end subroutine
!-------------

subroutine set_outputdir(this, outdir) 
	class(output_t),  intent(inout) :: this
	character(len=*), intent(in)    :: outdir
	
	this%outdir = trim(outdir) // "/"
		
	call make_dir(outdir) 
	
	call open_files(this) 
	
	call this%tdse%laser%write(outdir)
	
end subroutine 	
!-------------

subroutine open_files(this) 
	class(output_t), intent(inout) :: this

	call close_files(this)
  
    this%isOpened = .true.
  
	if( this%isout_laser ) &
		open(newunit = this%idf_laser, file=trim(this%outdir)//'elaser.dat', &
		     status='unknown', action = 'write')
	
	if( this%isout_dip ) &
		open(newunit = this%idf_dip, file=trim(this%outdir)//'dip.dat', &
		     status='unknown', action = 'write')
	
	if( this%isout_norm ) &
		open(newunit = this%idf_norm, file=trim(this%outdir)//'norm.dat', &
		     status='unknown', action = 'write')
		     
	if( this%isout_prob ) &
		open(newunit = this%idf_prob, file=trim(this%outdir)//'prob.dat', &
		     status='unknown', action = 'write')
		     
	if( this%isout_energy ) &
		open(newunit = this%idf_energy, file=trim(this%outdir)//'energy.dat', &
		     status='unknown', action = 'write')
		     
	if( this%isout_coeff ) &
		open(newunit = this%idf_coeff, file=trim(this%outdir)//'coeff.dat', &
		     status='unknown', action = 'write')
	
end subroutine 	
!-------------

subroutine close_files(this)
	class(output_t), intent(inout) :: this	
	
	if( this%isOpened ) then

		if( this%isout_laser  ) close(unit=this%idf_laser )
	
		if( this%isout_dip    ) close(unit=this%idf_dip   )
	
		if( this%isout_norm   ) close(unit=this%idf_norm  )
		
		if( this%isout_prob   ) close(unit=this%idf_prob  )
		
		if( this%isout_energy ) close(unit=this%idf_energy)
		
		if( this%isout_coeff ) close(unit=this%idf_coeff)
		
		this%isOpened = .false.
		
	endif
	
end subroutine	
!-------------

subroutine destroy_output(this)
	class(output_t), intent(inout) :: this

	call close_files(this) 
	
	this%dat  => null()
	this%tdse => null() 
	
end subroutine
!-------------

subroutine write_output(this, t, coeffs)
	class(output_t), intent(inout) :: this
	real(wp),        intent(in)    :: t
	type(coeffs_t),  intent(inout) :: coeffs
	
	real(wp) :: Ey, Ez
		
	if( .not. this%isOpened ) stop 'Write_output: File is opened.'

!----- write laser 
	if( this%isout_laser ) then 

	block
		call this%tdse%laser%elaser(t, Ey, Ez)
		write(this%idf_laser, '(3E16.8)') t, Ez, Ey
	end block
		
	endif 

!----- write normalization of wave function
	if( this%isout_norm ) then
	block
		real(wp) :: norm
		norm = sum( abs(coeffs%Ci)**2 )
		write(this%idf_norm, '(2E16.8)') t, norm
	end block
	endif

!----- write dipole
	if( this%isout_dip ) then 
	
	block
!~ 		real(wp) :: dy(1,1), dz(1,1)
		complex(wp) :: dy, dz
		complex(wp) :: tmpy, tmpz
		real(wp)    :: dippar, dipper
		
		integer     :: mi, mj, ni, nj
		
!~ 		dy = matmul( matmul(conjg(transpose(coeffs%Ci2d)), this%dat%dipy), coeffs%Ci2d ) + Ey
!~ 		dz = matmul( matmul(conjg(transpose(coeffs%Ci2d)), this%dat%dipz), coeffs%Ci2d ) + Ez
		dy = 0.0_wp
		dz = 0.0_wp
		!$OMP PARALLEL PRIVATE(MI,MJ, NI, NJ, TMPY, TMPZ) REDUCTION(+: DY, DZ)
		do mi = this%dat%mmin, this%dat%mmax
			!$OMP DO 
			do ni = 1, this%dat%nbasm(mi)
				tmpz = 0.0_wp 
				do nj = 1, this%dat%nbasm(mi)
				    tmpz =  tmpz + this%dat%dipz(mi,mi) % pn(nj,ni) * coeffs%Cnm(mi) % pn(nj)
				enddo
				!dzcitmp(ni) = tmp
			
				tmpy = 0.0_wp
				if (mi+1 <= this%dat%mmax) then
					do nj = 1, this%dat%nbasm(mi+1)
						tmpy = tmpy + this%dat%dipy(mi,mi+1) % pn(nj,ni) * coeffs%Cnm(mi+1) % pn(nj)
					enddo
				endif
				if (mi-1 >= this%dat%mmin) then
					do nj = 1, this%dat%nbasm(mi-1)
!~ 						tmpy = tmpy + this%dipy(mi,mi-1) % pn(nj,ni) * this%coeffs%Cnm(mi-1) % pn(nj)
						tmpy = tmpy + conjg(this%dat%dipy(mi-1,mi) % pn(ni,nj)) * coeffs%Cnm(mi-1) % pn(nj)
					enddo
				endif
				dy = dy + conjg(coeffs%Cnm(mi)%pn(ni) ) * tmpy
				dz = dz + conjg(coeffs%Cnm(mi)%pn(ni) ) * tmpz 
			enddo
			!$OMP END DO NOWAIT 
		enddo
		!$OMP END PARALLEL
		
!~ 		!$OMP PARALLEL DO PRIVATE(I, J, DYTMP, DZTMP) REDUCTION (+:DY, DZ)
!~ 		do i = 1, this%dat%nbas
!~ 		  dytmp = 0.0_wp
!~ 		  dztmp = 0.0_wp
!~ 		  do j = 1, this%dat%nbas
!~ 			dytmp = dytmp + coeffs%Ci(j) * this%dat%dipy(i,j)
!~ 			dztmp = dztmp + coeffs%Ci(j) * this%dat%dipz(i,j)
!~ 		  enddo

!~ 			dy = dy + dytmp * conjg(coeffs%Ci(i))
!~ 			dz = dz + dztmp * conjg(coeffs%Ci(i))
!~ 		enddo
!~ 		!$OMP END PARALLEL DO
		
!~ 		dy = dy + Ey 
!~ 		dz = dz + Ez
!~ 		dy = dy + Ey * sum(abs(coeffs%Ci)**2 )
!~ 		dz = dz + Ez * sum(abs(coeffs%Ci)**2 )
		
		dippar = dz * this%tdse%laser%cstheta + dy * this%tdse%laser%sntheta
		dipper = dy * this%tdse%laser%cstheta - dz * this%tdse%laser%sntheta
		
		write(this%idf_dip, '(7E16.8)') t, dippar, dipper, &
										real(dz), imag(dz), real(dy), imag(dy)
	end block	
	endif 

!----- write probability ionization
	if( this%isout_prob ) then
	block
		real(wp) :: prob, prob_gs
		prob    = 1.0_wp - sum( abs(coeffs%Ci)**2, mask = this%dat%Ei < 0.0_wp ) 
		prob_gs = 1.0_wp - abs( coeffs%Cnm(coeffs%mini) % pn(coeffs%nini) )**2 
		write(this%idf_prob, '(3E16.8)') t, prob, prob_gs
	end block
	endif

!----- write energy under laser-dressed
	if( this%isout_energy) then
!~ 		write(this%idf_energy, '(2E16.8)') t, this%dat%Enm(coeffs%mini) % pn(coeffs%nini) 
!~ 		write(this%idf_energy, '(2E16.8)') t, abs(coeffs%Cnm(coeffs%mini) %pn(coeffs%nini))**2 &
!~ 													 * this%dat%Enm(coeffs%mini) %pn(coeffs%nini) 
	  call calc_energy_1()

	endif

!----- write coeffs of bound states under laser-dressed
	if( this%isout_coeff ) then
	
	block
	  integer :: mi, ni
		do mi = this%dat%mmin, this%dat%mmax
		  
		  do ni = 1, this%dat%nbasm(mi)
			if( this%dat%Enm(mi) % pn(ni)  < 0.0_wp ) then
			  write(*,*) mi, ni
		      write(this%idf_coeff, '(*(E16.8))') t, coeffs%Cnm(mi) % pn(ni)
		    endif
		  enddo
		  write(this%idf_coeff,*)
		  
		enddo
	end block
	
	endif

!------ ... 

contains
!-------

	subroutine calc_energy_1()
	  real(wp)    :: H0exval, Vexval
	  complex(wp) :: tmpy, tmpz
	  
	  integer     :: mi, mj, ni, nj
	  
	  H0exval = sum(abs(coeffs%Ci)**2 * this%dat%Ei, mask = this%dat%Ei <=0.0_wp )
	  
	  Vexval = 0.0_wp
	  
	  !$OMP PARALLEL PRIVATE(MI,MJ, NI, NJ, TMPY, TMPZ) REDUCTION(+: Vexval)
	  do mi = this%dat%mmin, this%dat%mmax
		!$OMP DO 
		do ni = 1, this%dat%nbasm(mi)
			if (this%dat%Enm(mi) % pn(ni) > 0.0_wp) cycle
			tmpz = 0.0_wp 
			do nj = 1, this%dat%nbasm(mi)
			  if (this%dat%Enm(mi) % pn(nj) <= 0.0_wp) &
				tmpz =  tmpz + this%tdse%dipz(mi,mi) % pn(nj,ni) * coeffs%Cnm(mi) % pn(nj)
			enddo
			!dzcitmp(ni) = tmp
		
			tmpy = 0.0_wp
			if (mi+1 <= this%dat%mmax) then
				do nj = 1, this%dat%nbasm(mi+1)
				  if (this%dat%Enm(mi+1) % pn(nj) <= 0.0_wp) &
					tmpy = tmpy + this%tdse%dipy(mi,mi+1) % pn(nj,ni) * coeffs%Cnm(mi+1) % pn(nj)
				enddo
			endif
			if (mi-1 >= this%dat%mmin) then
				do nj = 1, this%dat%nbasm(mi-1)
				   if (this%dat%Enm(mi-1) % pn(nj) <= 0.0_wp) &
	!~ 						tmpy = tmpy + this%dipy(mi,mi-1) % pn(nj,ni) * this%coeffs%Cnm(mi-1) % pn(nj)
					 tmpy = tmpy + conjg(this%tdse%dipy(mi-1,mi) % pn(ni,nj)) * coeffs%Cnm(mi-1) % pn(nj)
				enddo
			endif
			
			Vexval = Vexval + conjg(coeffs%Cnm(mi)%pn(ni) ) * (tmpy*Ey + tmpz*Ez)
		
		 enddo
		 !$OMP END DO NOWAIT 
	  enddo
	  !$OMP END PARALLEL
	  write(this%idf_energy, '(3E16.8)') t, H0exval, H0exval + Vexval
	  
	end subroutine
	!-------------
	
	subroutine calc_energy_1p()
	  real(wp)    :: H0exval, Vexval
	  complex(wp) :: tmpy, tmpz
	  integer     :: mi, mj, ni, nj
	  
	  H0exval = abs(coeffs%Cnm(coeffs%mini) % pn(coeffs%nini))**2 &
			  * this%dat%Enm(coeffs%mini)   % pn(coeffs%nini) 
	  
	  Vexval = 0.0_wp
	  
	  !$OMP PARALLEL PRIVATE(MI,MJ, NI, NJ, TMPY, TMPZ) REDUCTION(+: Vexval)
	  do mi = this%dat%mmin, this%dat%mmax
		!$OMP DO 
		do ni = 1, this%dat%nbasm(mi)
			if ( mi /= coeffs%mini .or. ni /= coeffs%nini) cycle
			tmpz = 0.0_wp 
			do nj = 1, this%dat%nbasm(mi)
			  if (nj == coeffs%nini) &
				tmpz =  tmpz + this%tdse%dipz(mi,mi) % pn(nj,ni) * coeffs%Cnm(mi) % pn(nj)
			enddo
		
!~ 			tmpy = 0.0_wp
!~ 			if (mi+1 <= this%dat%mmax) then
!~ 				do nj = 1, this%dat%nbasm(mi+1)
!~ 				  if (this%dat%Enm(mi+1) % pn(nj) <= 0.0_wp) &
!~ 					tmpy = tmpy + this%tdse%dipy(mi,mi+1) % pn(nj,ni) * coeffs%Cnm(mi+1) % pn(nj)
!~ 				enddo
!~ 			endif
!~ 			if (mi-1 >= this%dat%mmin) then
!~ 				do nj = 1, this%dat%nbasm(mi-1)
!~ 				   if (this%dat%Enm(mi-1) % pn(nj) <= 0.0_wp) &
!~ 					tmpy = tmpy + conjg(this%tdse%dipy(mi-1,mi) % pn(ni,nj)) * coeffs%Cnm(mi-1) % pn(nj)
!~ 				enddo
!~ 			endif
			
			Vexval = Vexval + conjg(coeffs%Cnm(mi)%pn(ni) ) * (tmpy*Ey + tmpz*Ez)
		
		 enddo
		 !$OMP END DO NOWAIT 
	  enddo
	  !$OMP END PARALLEL
	  write(this%idf_energy, '(3E16.8)') t, H0exval, H0exval + Vexval
	  
	end subroutine
	!-------------

	subroutine calc_energy_2()
	  real(wp)    :: H0exval, Vexval
	  complex(wp) :: tmpy, tmpz
	  integer     :: mi, mj, ni, nj
	  
	  H0exval = sum(abs(coeffs%Ci)**2 * this%dat%Ei )
	  
	  Vexval = 0.0_wp
	  
	  !$OMP PARALLEL PRIVATE(MI, MJ, NI, NJ, TMPY, TMPZ) REDUCTION(+: Vexval)
	  do mi = this%dat%mmin, this%dat%mmax
		!$OMP DO 
		do ni = 1, this%dat%nbasm(mi)
			tmpz = 0.0_wp 
			do nj = 1, this%dat%nbasm(mi)
				tmpz =  tmpz + this%tdse%dipz(mi,mi) % pn(nj,ni) * coeffs%Cnm(mi) % pn(nj)
			enddo
		
			tmpy = 0.0_wp
			if (mi+1 <= this%dat%mmax) then
				do nj = 1, this%dat%nbasm(mi+1)
					tmpy = tmpy + this%tdse%dipy(mi,mi+1) % pn(nj,ni) * coeffs%Cnm(mi+1) % pn(nj)
				enddo
			endif
			if (mi-1 >= this%dat%mmin) then
				do nj = 1, this%dat%nbasm(mi-1)
					tmpy = tmpy + conjg(this%tdse%dipy(mi-1,mi) % pn(ni,nj)) * coeffs%Cnm(mi-1) % pn(nj)
				enddo
			endif
			
			Vexval = Vexval + conjg(coeffs%Cnm(mi)%pn(ni) ) * (tmpy*Ey + tmpz*Ez)
		
		 enddo
		 !$OMP END DO NOWAIT 
	  enddo
	  !$OMP END PARALLEL
	  write(this%idf_energy, '(3E16.8)') t, H0exval, H0exval + Vexval

	end subroutine
	!-------------
	
end subroutine
!-------------

subroutine make_dir( dir )
	character(len=*), intent(in) :: dir
	
	call execute_command_line( 'mkdir -p ' // dir )
	
end subroutine
!-------------

end module
!*********
