program interp_potential
	use type_vars
	use consts
	use read_pars
	use grid
	use linear_interp
	use cubic_spline
	use grad_pot
	implicit none
	
	integer :: nr1, ne1
	real(dp), allocatable :: r1(:), e1(:), pot0(:,:), poti(:,:)
	real(dp), allocatable :: gradv_y0(:,:), gradv_z0(:,:)
	
	real(dp), allocatable :: gradv_yi(:,:), gradv_zi(:,:)
	
!---

	call read_input()
	call read_pot0()
	
	call create_grids()
	call prt_grids()
!stop	
	call cre_pot()
	call pot_interp()
	call wrt_interp_pot()
	
	call gradv_interp()
	call wrt_interp_gradv()
	
	call wrt_gradpot()
	
	call clr_vars()
	call clr_grids()
	
contains	
!*******

subroutine read_pot0()
	character*100 :: sdum
	integer :: i, j, nc

    nr1 = 0
	open(unit = 1, file = './original/POT.DAT')
	do
		read(1,'(a)') sdum
		if (sdum == "") exit
		nr1 = nr1 + 1
	enddo
	close(1)
!--
	nc = 0
	open(unit = 1,file = './original/POT.DAT')
	do
		read(1,*, end = 999)
		nc = nc + 1
		ne1 = nc / (nr1+1)
	enddo
999		rewind 1
	close(1)

!--
	allocate ( r1(nr1), e1(ne1),  pot0(nr1,ne1) )
	allocate ( gradv_y0(nr1,ne1), gradv_z0(nr1,ne1) )
	
	open(unit = 1, file = './original/POT.DAT')	
	do i = 1, Ne1	
		do j = 1, Nr1
			read(1,*) r1(j), e1(i), pot0(j,i)
		enddo
		read(1,*)
	enddo	
	close(1)


	open(unit = 1, file = './original/GRADV.DAT')	
	do i = 1, Ne1	
		do j = 1, Nr1
			read(1,*) r1(j), e1(i), gradv_y0(j,i), gradv_z0(j,i)
		enddo
		read(1,*)
	enddo	
	close(1)
	
!-- print to screen parameters of initial mesh

	write (*, '(a)') '=== FORMER GRIDS ==='
    write (*, '("NR1  = ", i11)') nr1
    write (*, '("NE1  = ", i11)') ne1 
    write (*, '("Rmax = ", f11.6 )') r1(nr1)
    write (*, '(a)')

end subroutine
!-------------

subroutine cre_pot()
	integer :: ie, ir

	allocate ( poti(nrad,neta) )
	allocate ( gradv_yi(nrad,neta), gradv_zi(nrad,neta) )
	
!~ 	open (unit = 1, file = 'POTtest.DAT')
!~ 	do ie = 1, neta
!~ 		do ir = 1, nrad
!~ 			poti(ir,ie) = -1.0_dp / rad(ir)
!~ 			write(1,*) rad(ir), eta(ie), poti(ir,ie)
!~ 		enddo
!~ 		write(1,*)
!~ 	enddo
!~ 	close(1)
	
end subroutine
!-------------

subroutine pot_interp()
	real(dp), allocatable :: pott(:,:)
	
	integer :: i, ie
	
	allocate ( pott(nrad,ne1) )
	
	do ie = 1, ne1
		pott(:,ie) = -1.0_dp / rad(:)
	enddo
!-- 
	do i = 1, ne1	
		call poly(r1, pot0(:,i), rad, pott(:,i) )
!~ 		call cspl(r1, pot0(:,i), rad, pott(:,i) )
	enddo
	
	do i = 1, nrad
		call poly(e1, pott(i,:), eta, poti(i,:))
!~ 		call cspl(e1, pott(i,:), eta, poti(i,:))
	enddo
	
	deallocate(pott)

end subroutine
!-------------

subroutine gradv_interp()
	real(dp), allocatable :: gradv_yt(:,:), gradv_zt(:,:)
	
	integer :: ir, ie
	
	allocate ( gradv_yt(nrad,ne1), gradv_zt(nrad,ne1) )

	do ie = 1, ne1
		gradv_yt(:,ie) = 1.0_dp / rad**2 * sqrt( 1.0_dp - eta(ie)**2 )
		gradv_zt(:,ie) = 1.0_dp / rad**2 * eta(ie)
	enddo	

	do ie = 1, ne1
		call poly(r1, gradv_y0(:,ie), rad, gradv_yt(:,ie) )
		call poly(r1, gradv_z0(:,ie), rad, gradv_zt(:,ie) )
	enddo
	
	do ir = 1, nrad
		call poly(e1, gradv_yt(ir,:), eta, gradv_yi(ir,:) )
		call poly(e1, gradv_zt(ir,:), eta, gradv_zi(ir,:) )
	enddo
	
	deallocate (gradv_yt, gradv_zt)

end subroutine
!-------------

subroutine wrt_interp_pot()
	integer :: ir, ie
	
	open(unit = 1, file = './output/POT.DAT', status = 'unknown')
	do ie = 1, neta
		do ir = 1, nrad
			write(1,'(3(E16.8,x))') rad(ir), eta(ie), poti(ir, ie)
		enddo
		write(1,'')
	enddo
	close(1)
	
	open(unit = 2, file = './output/RPOT.DAT', status = 'unknown')
		write(2, '(4(E16.8,X))') (rad(ir), rad(ir) * poti(ir,1), &
								  rad(ir) * poti(ir,neta/2), rad(ir) * poti(ir,neta), ir=1,nrad )
	close(2)
	
end subroutine
!-------------

subroutine wrt_interp_gradv()
	integer :: ir, ie
	
	open (unit = 1, file = './output/GRADV.DAT')
	do ie = 1, neta
		do ir = 1, nrad
			write(1, '(*(E16.8,X))') rad(ir), eta(ie), gradv_yi(ir,ie), gradv_zi(ir,ie)
		enddo
		write(1,'()')
	enddo
	
	close(1)
	
end subroutine
!-------------

subroutine clr_vars()

	deallocate (r1, e1, pot0, poti)
	deallocate (gradv_y0, gradv_z0, gradv_yi, gradv_zi)

end subroutine
!-------------

end program
!##############################
!### CODE nhiệt tình quá ha ###
!##############################
