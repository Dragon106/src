!***********************************************************************
! MODULE CONTROLS INFORMATION ABOUT ATOMS, ELECTRON AND ELECTRON DENSITY 
!***********************************************************************

module eledens_info
	use type_vars
	use consts
	use grid
	use utils
	implicit none

	integer :: Natom, Units, Nele, Norb
	integer,  allocatable :: Nzatom(:)
	real(dp), allocatable :: xatom(:), yatom(:), zatom(:)

	real(dp), allocatable, dimension (:,:,:) :: wf, den1e
	real(dp), allocatable, dimension (:,:)   :: eden, seden, gr_seden, ge_seden
	real(dp), allocatable, dimension (:,:,:) :: drwf, dewf

contains
!*******

subroutine read_Rinfo()
	character :: symbol

	integer   :: i 
	
	open(unit=1, file = 'Rinfo.inp')

		read(1,*)
		read(1,*) NAtom, Units
		read(1,*)
		read(1,*)
		
		allocate( Nzatom(NAtom), xatom(NAtom), yatom(NAtom), zatom(NAtom) )
		
		Nele = 0
		
		do i = 1, NAtom
			read(1,*) symbol, Nzatom(i), xatom(i), yatom(i), zatom(i)
		enddo
			Nele = sum( Nzatom(:) )
	close(1)
	
	if (Units .eq. 2) then
		xatom(:) = xatom(:) / aulength
		yatom(:) = yatom(:) / aulength
		zatom(:) = zatom(:) / aulength
 	endif
	
	Norb = ceiling( dble(Nele)/2 )
	
end subroutine
!-------------

subroutine write_edens()
!	character, intent(in) :: wdir*(*) 
	integer :: ir, ie 

!	open(unit = 100, file = getfln(wdir,'edens.dat'))
	open(unit = 100, file = './output/edens.dat')

	do ie = 1, neta
		do ir = 1, nrad
			write(100, '(4(e,3x))') eden(ir,ie), seden(ir,ie), & 
									gr_seden(ir,ie), ge_seden(ir,ie)  
		enddo
		write(100, '')
	enddo

	close(100)

end subroutine	
!-------------

subroutine read_edens()
!	character, intent(in) :: fln*(*)
	integer :: ir, ie 

!	open(unit = 100, file = fln, status = 'old') 
	open(unit = 100, file = './output/edens.dat', status = 'old') 

		read(100, *) ( ( eden(ir,ie), seden(ir,ie), &
						gr_seden(ir,ie), ge_seden(ir,ie), ir = 1, nrad ), ie = 1, neta ) 

	close(100)

end subroutine
!-------------

subroutine pre_calc_dens()
	
	allocate ( den1e(Nrad, Neta, Norb), eden(Nrad, Neta), seden(Nrad, Neta), &
			   gr_seden(Nrad, Neta), ge_seden(Nrad, Neta) )
			   
end subroutine
!-------------

subroutine clr_dens()

	deallocate ( den1e, eden, seden, gr_seden, ge_seden)
	
end subroutine
!-------------

subroutine clr_atm_info() 

	deallocate ( Nzatom, xatom, yatom, zatom )

end subroutine
!-------------

end module
!*********
