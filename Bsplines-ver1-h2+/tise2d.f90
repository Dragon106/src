program tise2d
    use consts
    use mesh
    use utils 
    use read_pars
    use basis
    use molpot
    use hamil 
    use wfs
    use grad_pot !!!
    
    implicit none

    call prog_init()
!stop
    call cre_basis()
    call cre_overlap()
    call prt_basis()
    call wrt_overlap()
    call wrt_basis()

    call cre_ham()
    call ham_diag(coeffs,ovlp,engs)
    call Nmax_calc()
    call prt_eigenv()
 
    call calc_wf()
!    call check_norm()
    call wrt_eigenve()
!    call wrt_20wf()
!    call calc_dip()
        
    call prog_end()

contains
!*******

subroutine prog_init() 
    
    character :: fln*100

    call getarg(1,fln)
    
    if( fln == '' ) &
        call stop_error('Enter input file name!!!')

    call read_input(trim(fln))

    call cre_mesh()
    call prt_mesh()

    call cre_pot()

!--	read file POT.DAT of the molecule
!    call read_pot(trim(potfln))

!-- for hydrogen atom: calc potential
!~     call calc_pot()
!~     call wrt_grad()

!-- for H2+
	call calc_hydro_mol()
!~ 	call wrt_grad() !!! numerical derivatives
    
end subroutine 
!-------------

subroutine prog_end()

    call clr_ham()
    call clr_pot()
    call clr_overlap()
    call clr_mesh() 

end subroutine 
!-------------

subroutine calc_pot()

    integer :: i, j

!   forall (i = 1:nrad, j = 1:neta)
    forall (i = 1:nrad)

!-- calc pot of hydrogen atom
!       pot(i,j) = -1.0d0 / sqrt( ( rad(i)**2 * (1.0d0 - eta(j)**2) ) + (rad(i)*eta(j))**2 ) 
!       pot(i,j) = -1.0d0 / rad(i)  

        pot(i,:) = -1.0d0 / rad(i)

!-- calc pot of harmonic oscillator
!       pot(i,j) = 0.5d0 * rad(i)*rad(i)
    end forall

!-- write potential in file for TDSE --
    
    open (unit=1,file='POT.DAT')
		do j = 1, nrad
			do i = 1, neta
				write(1, '(3(E24.16,X))') rad(j), eta(i), pot(j,i)
			enddo
			write(1,*)
		enddo
    close(1)
 
end subroutine
!-------------

subroutine calc_hydro_mol()
	real(8), allocatable :: v1(:,:), v2(:,:)
	real(8), allocatable :: dr_v(:,:), de_v(:,:)

	real(8) :: z1, z2 !!! position of the nucleii

	integer :: ir, ie
	
	z1 =  R/2.0d0
	z2 = -R/2.0d0
	
	allocate ( v1(nrad,neta), v2(nrad,neta) )
	allocate ( dr_v(nrad,neta), de_v(nrad,neta) )
		
	forall (ir = 1:nrad, ie = 1:neta)
	
		v1(ir,ie) = -1.0d0 / sqrt ( rad(ir)**2 - 2* z1 * rad(ir)*eta(ie) + z1**2 )
		v2(ir,ie) = -1.0d0 / sqrt ( rad(ir)**2 - 2* z2 * rad(ir)*eta(ie) + z2**2 )
		pot(ir,ie)  = v1(ir,ie) + v2(ir,ie)
		
		dr_v(ir,ie) = - (rad(ir) - z1 * eta(ie)) * v1(ir,ie)**3 - (rad(ir) - z2 * eta(ie)) * v2(ir,ie)**3
		de_v(ir,ie) =  rad(ir) * ( z1 * v1(ir,ie)**3 + z2 * v2(ir,ie)**3 )
		
		grady_v(ir,ie) = sqrt( 1 - eta(ie)**2 ) * ( dr_v(ir,ie) - eta(ie) * de_v(ir,ie) / rad(ir) )
		gradz_v(ir,ie) = eta(ie) * dr_v(ir,ie) + ( 1-eta(ie)**2 ) * de_v(ir,ie) / rad(ir)
		
	end forall 	
	
	open (unit=1,file='POT.DAT')
	open (unit=2,file='GRADV.DAT')
		do ie = 1, neta
			do ir = 1, nrad
				write(1, '(3(E24.16,X))') rad(ir), eta(ie), pot(ir,ie)
				write(2, '(4(E24.16,X))') rad(ir), eta(ie), grady_v(ir,ie), gradz_v(ir,ie)
			enddo
			write(1,*)
			write(2,*)
		enddo
    close(1)
    close(2)

	deallocate (v1, v2)
	deallocate (dr_v, de_v)
	
end subroutine
!-------------

end program 
!**********
