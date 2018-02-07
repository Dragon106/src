!*** CALCULATION OF HHG FROM DIPOLE: z, dz, d2z, and ACC-DIP FILES ***

program calculation_hhg
	use pars
	use vars
	use read_inputs
	use deriv_dip
	use hhg_mod
	use utils
	use splines !!!
	implicit none
	
	real(dp), allocatable :: time(:)

	character :: FileName*50, inp_fln*50
	
	integer   :: itheta, i, ios

!--
	call get_command_argument(1, inp_fln)
!~ 		if (inp_fln == '') inp_fln = 'input'
		if (inp_fln == '') inp_fln = 'para_hhg.inp'
		
	call read_pars( inp_fln )
	
	write(*,'(a, f6.1, a)') '   Smooth HHG with sigma       = ', sigma, '*omega_laser'
	if ( isSpline3der ) then
		write(*,'(3x,A)') 'Derivatives use Splines3Der.'
	else
		write(*,'(3x,A)') 'Derivatives use finite differences.'
	endif
	write(*,*) '===================================================='

	call cre_vars()

!--------------------------------------------------------------------------

!-- write output file: hhg.dat and phase.dat in format: 
!   order, hhg_x, hhgs_x, hhg_y, hhgs_y
!   order, phase_x, phase_y
!--
 	
	do itheta = theta_min, theta_max, theta_step
	  theta =  itheta * Pi / 180.0_dp
	  write(*,'(A,I4,A)') '   Alignment angle             = ', itheta, ' degree(s)'

!				  
!------ calc hhg from dipole: z, dz, d2z or gradV			  
!	

	  allocate( time(0:Ntmax) )
		
	  if ( isout_dip .or. isout_vec .or. isout_acc ) then
	    
		call alloc_arr_xy( dip, 0, Ntmax )

!-- read dipole file 
!~ 		write(FileName, '(a,i0,a)') trim(dat_dir) // 'dip-',itheta,'.dat' !!! for dip-0.dat, dat-5.dat,...
		write(FileName, '(a,i0,a)') trim(dat_dir) // '/', itheta, '/' // 'dip.dat' !!! for './dat_dir/0/dip.dat',... 
		!--- accompanied with writing dz, d2z and turning on mp_wrt_hhg and mp_wrt_pha in hhg_mod.f90
		
		open(unit = 11, file = FileName, status = 'old', action = 'read')
!--	
		do i = 0, Ntmax
			read(11,*,iostat = ios) time(i), dip%x(i), dip%y(i)
			if (ios /= 0) then
			    write(*,*)
				write(*,'(3x,3(A,I0))') 'WARNING! File dip.dat ', itheta, ', Nt only has ', i-1, '/', Ntmax
			    write(*,*)
				Ntmax = i - 1
				exit
			endif
		enddo
		
		close(11)
		
!~ 		call alloc_arr_xy( dip, 0, Ntmax ) !!!
		call resize_arr_xy( dip, 0, Ntmax ) 
		
!~ 		write(*,*) size(dip%x) 
	  endif
!		
!-- calc hhg from z
!

      if( isout_dip ) then

!-- calc hhg signal (no w**4)
		call calc_hhg( order, dip, hhg, hhgodd, hhgs, phase )

!-- write hhg and phase
        call wrt_hhg( trim(dat_dir), 'hhg-z', 'hhgo-z', order, hhg, hhgodd, hhgs, itheta )
        call wrt_pha( trim(dat_dir), 'phase-z', order, phase, itheta )
        
      endif
    
!
!-- calc dz and hhg from dz	
!

      if (isout_vec) then 

!-- calc dz
        call alloc_arr_xy( dz, 0, Ntmax )
        
        if ( isSpline3der ) then
			call spline3ders(time, dip%x, time, dynew=dz%x)
			call spline3ders(time, dip%x, time, dynew=dz%x)
		else
			call deriv( TimeStep, dip%x, dz%x )
			call deriv( TimeStep, dip%y, dz%y )
		endif

!--	write dz	
!~ 		write(FileName,'(a,i0,a)') trim(dat_dir) // 'dz-',itheta,'.dat'
		write(FileName,'(a,i0,a)') trim(dat_dir) // '/', itheta, '/' // 'dz.dat' !!! 
		open(unit = 1, file = FileName)
			write(1,'(3(E16.8,x))') ( time(i), dz%x(i), dz%y(i), i = 0, Ntmax )
		close(1)	

!-- calc hhg signal (no w**2) from dz		
		call calc_hhg( order, dz, hhg, hhgodd, hhgs, phase )

!-- write hhg and phase		
        call wrt_hhg( trim(dat_dir), 'hhg-dz', 'hhgo-dz', order, hhg, hhgodd, hhgs, itheta )
        call wrt_pha( trim(dat_dir), 'phase-dz', order, phase, itheta )
        
        call dealloc_arr_xy( dz )
        
      endif

!    
!-- calc d2z and hhg from d2z
!
    
      if ( isout_acc ) then
        call alloc_arr_xy( d2z, 0, Ntmax )
        
!-- calc d2z
!		write(*,*) size(dip%x), size(d2z%x)
		if (isSpline3der) then
			call spline3ders(time, dip%x, time, d2ynew=d2z%x)
			call spline3ders(time, dip%y, time, d2ynew=d2z%y)
		else
			call second_deriv( TimeStep, dip%x, d2z%x )
			call second_deriv( TimeStep, dip%y, d2z%y )
		endif
		
!--	write d2z	
!~ 		write(FileName, '(a,i0,a)') trim(dat_dir) // 'd2z-',itheta,'.dat'
		write(FileName, '(a,i0,a)') trim(dat_dir) // '/', itheta, '/' // 'd2z.dat'
		open(unit = 2, file = FileName)
			write(2,'(3(E16.8,x))') ( time(i), d2z%x(i), d2z%y(i), i = 0, Ntmax )
		close(2)
			
!--	calc hhg from d2z		
		call calc_hhg(order, d2z, hhg, hhgodd, hhgs, phase)

!-- write hhg and phase from d2z
        call wrt_hhg( trim(dat_dir), 'hhg-d2z', 'hhgo-d2z', order, hhg, hhgodd, hhgs, itheta )
        call wrt_pha( trim(dat_dir), 'phase-d2z', order, phase, itheta )

        call dealloc_arr_xy( d2z )

      endif
    
!
!------ calc hhg from acceleration dipole
!
    
      if ( isout_grd ) then		
        call alloc_arr_xy( acc_dip, 0, Ntmax )
		
!-- read accleration file
		write(FileName,'(a,i0,a)') trim(dat_dir) // 'acc-dip-',itheta,'.dat'
		open(unit = 111, file = FileName)
!--		
		do i = 0, Ntmax
			read(111,*,iostat = ios) time(i), acc_dip%x(i), acc_dip%y(i)
			if (ios .ne. 0) then
				write(*,'(3x,3(A,I0))') 'WARNING! File acc-dip.dat ', itheta, ', Nt only has ', i-1, '/', Ntmax
				write(*,*)
				Ntmax = i - 1
				exit
			endif
		enddo
		
		close(111)
		
!-- calc hhg from acc-dip (Ehrenfest's theorem)		
		call calc_hhg( order, acc_dip, hhg, hhgodd, hhgs, phase )

!-- write hhg and phase		
        call wrt_hhg( trim(dat_dir), 'hhg-acc', 'hhgo-acc', order, hhg, hhgodd, hhgs, itheta )
        call wrt_pha( trim(dat_dir), 'phase-acc', order, phase, itheta )

	  endif
!--

	  deallocate( time )
			
    enddo

!--------------------------------------------------------------------------

	call clr_vars()

end program
!**********
