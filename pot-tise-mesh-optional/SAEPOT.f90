! **********************************************************************
!  PROGRAM FOR CALCULATING SAE POTENTIAL FOR A LINEAR MOLECULE		   !
!  INPUT : WF FROM TISE (WF.in), or WF FROM GAUSSIAN (input-io.wf)	   !			
!  OUTPUT: SAE POTENTIAL: POT.DAT, RPOT.DAT and edens.dat			   !
! **********************************************************************

program saepot_calc
	use type_vars
	use consts
	use read_pars
	use grid
	use eledens_info
	use edens_tise
	use edens_gau
	use grad_pot
		
	implicit none
	
	real(dp), allocatable :: pot(:,:)
	real(dp) :: r, e, rp1, rp2, rp3 

!---

	call read_Rinfo
	call read_input()
	call create_grids()	
	call prt_grids()
	call prt_pars()
	
	call pre_calc_dens()
	
!-- calc electron density
!~  	call calc_dens_gau()  ! from Gaussian
 	call calc_dens_tise() ! from TISE wf

	call read_edens()     ! read eden, seden, gr_seden, ge_eden   
	
	call cre_pot()
	call calc_pot()
	call wrt_pot()
	call wrt_gradpot() !!!
	call pot_fin()

!***********************************************************************
	
contains
!-------

subroutine calc_pot()
!-- calc SAE potential
	integer :: ir, ie
	
	do ie = 1, Neta
		e = eta(ie)		
		do ir = 1, Nrad
			r = rad(ir)
			if (r .LE. rasym) then
				pot(ir,ie) =  saepot(ir,ie)
			else
				pot(ir,ie) = -1.d0 / r
			endif
		enddo
	enddo
	
end subroutine
!-------------

subroutine wrt_pot()
	integer :: ir, ie

!-- write SAE potential in file	
	open (unit = 1, file = './output/POT.DAT')
	do ie = 1, neta
		do ir = 1, nrad
			write(1, '(3(E16.8,X))') rad(ir), eta(ie), pot(ir,ie)
		enddo	
		write(1,'')
	enddo
	close(1)

!-- write effective charge in file	
	open(unit = 2, file = './output/RPOT.DAT')
	do ir = 1, Nrad
		r  = rad(ir)
		e = eta(1)
		rp1 = r * pot(ir,1)
		e = eta( Neta/2 )
		rp2 = r * pot(ir, neta/2)
		e = eta(neta)
		rp3 = r * pot(ir, neta)
		write(2,'(4(E16.8,X))') r, rp1, rp2, rp3
	enddo
	close(2)
	
end subroutine
!-------------

! -- FUNCTION TO GET SAE POTENTIAL
!    V(r) = V_en(r) + V_el(r) + V_exc(r) + V_cor(r)

function saepot(i,j)
	integer, intent(in) :: i, j
	real(dp) :: saepot
	
	saepot = enpot() + harpot() + excpot(i,j) + iVc*corpot(i,j) ! Notice ()! () is important in case of FUNCTION!
	
end function
!-----------

function enpot
	real(dp) :: enpot
	
	enpot = - sum( Nzatom(:) / sqrt(R**2 + Zatom(:)**2 - 2.0_dp* R * Zatom(:) * E) )
		
end function
!-----------

function harpot
	real(dp) :: harpot
	
	real(dp) :: temp, Rp, Ep, Angle, Ree, wre ! p = prime
	
	integer :: i, j, k

	temp = 0.0_dp

	do i = 1, Nrad
		Rp = rad(i) ! r'
		do j = 1, Neta
			Ep = eta(j) ! eta'
			do k = 1, NPGAU
				Angle = sqrt(1 - E**2) * sqrt(1 - Ep**2) * cos(PGAU(k)) + E*Ep! Angle between r and r'
				Ree  = sqrt(R**2 + Rp**2 - 2.0_dp * R * Rp * Angle) ! |r - r'|
				
!~ 				if (Ree .ne. 0.0_dp) then
				if (Ree .GT. 1.0d-10) then
				wre = w12(i,j)
				temp = temp + wre * WPGAU(k) * EDen(i,j) / Ree
				else
					stop
				endif
							
			enddo
		enddo
	enddo

	harpot = temp	

end function
!-----------

function excpot(i,j)	
	integer, intent(in) :: i, j ! i = 1, NDVR2; j = 1, NDVR1
	real(dp) :: excpot
	real(dp) :: chi, VLDA, VGC, sp43, gra
	
    real(dp) :: eps, sp
    
    eps = tiny(0.0d0) !epsilon(0.d0)
    
    sp = seden(i,j)
    
    if (sp .LE. eps) then 
        sp = eps
        gra = 0.0d0 
    else
        gra  = sqrt( gr_seden(i,j)**2 + ge_seden(i,j)**2 )
    endif             
        
    sp43 = sp**(4.d0/3.d0)
    
!    if( sp43 .le. eps ) sp43 = eps
!    if( gra  .le. eps ) gra  = eps
    
	chi    =  gra / sp43  

	VLDA   = - ALPHA * ( 6.d0/pi * sp )**(1.d0/3.d0)
	VGC    = - BETA * chi**2 * sp**(1.d0/3.d0) / ( 1 + 3.d0 * BETA * chi * ASINH(chi) )

	excpot = VLDA + VGC
	
end function
!-----------

function corpot(i,j)
    integer, intent(in) :: i, j
    real(dp) :: corpot

!-- Using ferromagnetic parameters	
    real(dp), parameter :: a = 0.0621814/2
    real(dp), parameter :: b = 7.06042
    real(dp), parameter :: c = 18.057
    real(dp), parameter :: x0 = -0.32500
    
    real(dp) :: q
    
    real(dp) :: rs
    real(dp) :: ec, x, xx, xx0
    
    rs = (3.0d0 / 4.0d0 / pi / seden(i,j) )**(1.0d0/3.0d0)
    
    x = sqrt(rs)
    Xx = rs + b*x + c
    Xx0 = x0**2 + b*x0 + c
    Q = sqrt(4*c - b**2)

!-- In Hartree unit, the formula (in Rydberg unit) must be added the factor 1/2
    ec = 0.5 * a * ( log(x**2 / Xx) + 2*b/Q * atan(Q / (2*x + b) ) &
                   - b*x0 / Xx0 * (log( (x-x0)**2 / Xx ) + 2*(b + 2*x0) /Q *atan( Q/(2*x + b) ) ) )

    corpot = ec - a/6 * ( (x - x0)*c - b*x*x0 ) / ( (x - x0) * Xx )
   
end function
!-----------


subroutine cre_pot()

	allocate ( pot(nrad,neta) )

end subroutine
!-------------

subroutine pot_fin()

	deallocate ( pot )
	
	call clr_grids()
	call clr_dens()
	call clr_atm_info()
	
end subroutine
!-------------

end program
!**********

