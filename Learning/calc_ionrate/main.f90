!***********************************************************************
! CALCULATION OF IONIZATION RATE FROM IONIZATION PROBABILITY BY USING  !
! THE FORMULA OF TELNOV [PRA2014]									   !
!**********************************************************************!
program main
	use type_vars_m
	use derivfx_m
	use calc_rate_m
	implicit none
	
	type(ioniz_t) :: ionrate
	
	integer :: i
	
!~ 	call ionrate%read_para()
!~ 	call ionrate%read_ntmax()
!~ 	call ionrate%get_timestep()
	call ionrate%pre_calc()
	
	do i = ionrate%thetamin, ionrate%thetamax, ionrate%thestep
	
		write(ionrate%dataname, '(A,I0,A)') trim(ionrate%datadir) &
											// '/', i, '/' // 'prob.dat'
		write(ionrate%outname,  '(A,I0,A)') trim(ionrate%datadir) &
											// '/', i, '/' // 'rate.dat'
		
!~ 		call ionrate%generate_vars()
!~ 		call ionrate%read_data() 
!~ 		call ionrate%calc_ion_rate()
!~ 		call ionrate%write()
		call ionrate%run()
	
	enddo
	
end program
!**********
