!***********************************************************************
!	CALCULATING ELECTRON AND SPIN DENSITY FROM GAUSSIAN WAVE FUNCTIONS !
!   INPUT: input-io.wf from Gaussian (io = 1:Norb)					   !	
!   OUTPUT: edens.dat 												   !
!***********************************************************************

module edens_gau
	use type_vars
	use consts
	use grid
	use eledens_info
	use basis_gauss
	use utils
	implicit none
	
	integer :: Nbasg
	
	integer,  dimension(:,:), allocatable :: iatom, nxg, nyg, nzg
	real(dp), dimension(:,:), allocatable :: coeffg, zetag
		
contains
!*******

subroutine calc_dens_gau()

	call read_gauss_coeffs()	
	call dens_gau()
	call checknorm_gau()
	call write_edens()
	
end subroutine
!-------------

subroutine read_gauss_coeffs()
	character :: filename*100
	integer   :: ibs, io, ios, idum
	
	write (*,'(a)')'=== READ WF FROM GAUSSIAN ==='

	do io = 1, Norb
		Nbasg = 0
		call make_filename(filename, 'input-gauss/input-',io,'.wf')
		open(unit = 2, file = filename)
		do
			read(2, *, iostat = ios)
			if ( ios /= 0 ) exit
			Nbasg = Nbasg + 1
		enddo
			
		write(*,'(A, I3, A, I4)') 'Nbas of orbital ', io, ': ', Nbasg
		
		close(2)		
	enddo
	
	write (*,'') 
!--

	allocate( iatom(nbasg,norb), nxg(nbasg,norb), nyg(nbasg,norb), & 
              nzg(nbasg,norb), coeffg(nbasg,norb), zetag(nbasg,norb) )
              
	do io = 1, Norb
		call make_filename(filename, 'input-gauss/input-',io,'.wf')
		open(unit=2, file = filename, status = 'old')
		do ibs = 1, Nbasg
			read(2,*) idum, iatom(ibs,io), nxg(ibs,io), nyg(ibs,io), nzg(ibs,io), &
					  coeffg(ibs,io), zetag(ibs,io)
		enddo
		close(2)		
	enddo
		
end subroutine	
!-------------

subroutine dens_gau()
	real(dp), dimension(:,:), allocatable :: rp, ep
	real(dp), dimension(:,:), allocatable :: dr_rp, de_rp, dr_ep, de_ep

	real(dp), dimension(:),   allocatable :: gaur, gaue, dgaur, dgaue 
	real(dp), dimension(:,:), allocatable :: gaure, drgaure, degaure
		
	real(dp), dimension(:),   allocatable :: den, grse_tmp, gese_tmp
	real(dp), dimension(:,:,:), allocatable :: grse_io, gese_io

	real(dp) :: keta(nrad, neta), int_phi
	
	real(dp) :: za, zeta
	integer  :: nx, ny, nz, n2
	integer  :: i, ir, ie, ia, io, ibs, jbs
 
	n2 = nrad * neta
	
	allocate ( rp(n2,natom), ep(n2,natom) )
	allocate ( dr_rp(n2,natom), de_rp(n2,natom), dr_ep(n2,natom), de_ep(n2,natom) )
	
	allocate ( gaur(n2), gaue(n2), dgaur(n2), dgaue(n2) )
	allocate ( gaure(n2,nbasg), drgaure(n2,nbasg), degaure(n2,nbasg) )
	
	allocate ( den(n2), grse_tmp(n2), gese_tmp(n2) )
	allocate ( grse_io(nrad,neta,norb), gese_io(nrad,neta,norb) )

!-----------

    do ia = 1, natom
        za   = zatom(ia)
        
		do ie = 1, neta
		do ir = 1, nrad
			i = ir + (ie-1) * nrad
			rp(i,ia) = sqrt( rad(ir) * rad(ir) - 2.0_dp * rad(ir) * za * eta(ie) + za*za )
			ep(i,ia) = (rad(ir) * eta(ie) - za ) / rp(i,ia) 
			
			dr_rp(i,ia) = ( rad(ir) - za*eta(ie) ) / rp(i,ia)
			de_rp(i,ia) = - rad(ir)*za / rp(i,ia)
			
			dr_ep(i,ia) = eta(ie) / rp(i,ia) - 1.0_dp / ( rp(i,ia)*rp(i,ia) ) &
					      * ( rad(ir)*eta(ie) - za ) * dr_rp(i,ia)
			de_ep(i,ia) = rad(ir) / rp(i,ia) - 1.0_dp / ( rp(i,ia)*rp(i,ia) ) &
					      * ( rad(ir)*eta(ie) - za ) * de_rp(i,ia)
					                         
		enddo
		enddo
    enddo
    
!------------
	
	do io = 1, Norb
!---
		do ibs = 1, Nbasg

		    ia   = iatom(ibs,io)
			zeta = zetag(ibs,io)
			nx   =   nxg(ibs,io) 
			ny   =   nyg(ibs,io)
			nz   =   nzg(ibs,io)
			
			call  gaurad( gaur, rp(:,ia), zeta, nx, ny, nz)
			call dgaurad(dgaur, rp(:,ia), zeta, nx, ny, nz)

			call  gaueta( gaue, ep(:,ia), zeta, nx, ny, nz) 
			call dgaueta(dgaue, ep(:,ia), zeta, nx, ny, nz)

			gaure(:,ibs) = gaur(:) * gaue(:) 

			drgaure(:,ibs) = gaue * dgaur * dr_rp(:,ia) + gaur * dgaue * dr_ep(:,ia)
			degaure(:,ibs) = gaue * dgaur * de_rp(:,ia) + gaur * dgaue * de_ep(:,ia) 	
									
		enddo
!---

		den = 0.0d0
		grse_tmp = 0.0d0
		gese_tmp = 0.0d0
		
		do ibs = 1, Nbasg
			do jbs = 1, Nbasg
			
                int_phi =  ditgsc( nyg(ibs,io) + nyg(jbs,io), nxg(ibs,io) + nxg(jbs,io) )

				den(:)   = den(:) +  coeffg(ibs,io) * gaure(:,ibs) &
						  * coeffg(jbs,io) * gaure(:,jbs) &
						  * int_phi
						  
				grse_tmp = grse_tmp + 2*coeffg(ibs,io) * coeffg(jbs,io) & 
						  * gaure(:,jbs) * drgaure(:,ibs)  &
						  * int_phi
						  
 				gese_tmp = gese_tmp + 2*coeffg(ibs,io) * coeffg(jbs,io) &
 						  * gaure(:,jbs) * degaure(:,ibs) &
 						  * int_phi

			enddo
		enddo
!--
		
		den1e(:,:,io)   = reshape (den, [nrad, neta])	
		grse_io(:,:,io) = reshape (grse_tmp, [nrad, neta])
 		gese_io(:,:,io) = reshape (gese_tmp, [nrad, neta]) 
		
	enddo
	
!------------ 

	eden  = 2 * sum ( den1e(:,:,1:Norb-1), 3) + ( 2-mod(Nele,2) ) * den1e(:,:,Norb)	
	seden = sum (den1e, 3)
!--	
	forall (ir = 1:nrad, ie = 1:neta)
		keta(ir,ie) = - sqrt (1-eta(ie)**2) / rad(ir)
	endforall

	gr_seden = sum (grse_io, 3)
 	ge_seden = sum (gese_io, 3) * keta

end subroutine
!-------------

subroutine checknorm_gau() 	

    write(*,'(a,f6.2)') 'Number of electrons = ', sum(w12*eden) *2*pi    

end subroutine
!-------------

function ditgsc(n1,n2)

!   ditgsc = integrate ( sin(phi)**n1 * cos(phi)**n2 * d(phi), 0, 2.d0*pi ) / (2*pi) )

	integer, intent(in) :: n1, n2
	real(dp) :: ditgsc  

	if( mod(n1,2) /= 0 .or. mod(n2,2) /= 0 ) then
		ditgsc = 0.0_dp 
	else
		ditgsc = exp( log_gamma(n1 + 1.0_dp) + log_gamma(n2 + 1.0_dp) - log_gamma(n1/2+n2/2 + 1.0_dp) & 
					 -log_gamma(n1/2 + 1.0_dp) - log_gamma(n2/2 + 1.0_dp) - (n1+n2)*log(2.0_dp) )
	endif

end function
!----------- 

end module
!*********
