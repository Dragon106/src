program main
	use type_vars
	use solve_tdse_m
	implicit none
	
	character(len = 100) :: input_fln, datdir, outdir

	call read_arguments(input_fln, datdir, outdir)
	call run(input_fln, datdir, outdir)

contains
!*******

subroutine read_arguments(input_fln, datdir, outdir)
	character(len = *), intent(out) :: input_fln, datdir, outdir
	
	call get_command_argument(1, input_fln)
		if (input_fln == '') input_fln = 'para_data.inp'
		
	call get_command_argument(2, datdir)
		if (datdir == '')    datdir    = './input-tise/'
	
	call get_command_argument(3, outdir)
	    if (outdir == '')    outdir    = './output-tdse/'
	
end subroutine
!-------------

subroutine run(input_fln, datdir, outdir) 
	character(len=*), intent(in) :: input_fln, datdir, outdir

	character(len = 100) :: dirname
	integer :: thetamin, thetamax, theStep
	
	integer :: itheta

	type(solve_tdse_t) :: solve_tdse

	call solve_tdse%pre(input_fln, datdir) 

	call read_orient(input_fln, thetamin, thetamax, theStep)
	
	do itheta = thetamin, thetamax, theStep
		
		write(*,'(A)')     repeat("*",63)
		write(*,'(A,I0)') 'THETA = ', itheta
		
		call solve_tdse%reset()
		
		call solve_tdse%tdse%laser%setup( real(itheta,wp) ) 
		
		write(dirname, '(A,I0)') trim(outdir) // "/", itheta 
		
		call solve_tdse%output%setdir( dirname ) 
		
		call solve_tdse%run()
		
	enddo
	
	call solve_tdse%destroy()
	
end subroutine
!-------------

subroutine read_orient(fln, thetamin, thetamax, theStep)
	character(len=*), intent(in) :: fln
	
	integer :: thetamin, thetamax, theStep
	
	integer :: idf
	
	namelist /orient_info/ thetamin, thetamax, theStep
	
	open (newunit = idf, file = fln, status = 'old', action = 'read')
		read(unit = idf, nml = orient_info)
	close(unit = idf)
	
end subroutine
!-------------
	
end program
