module read_inputs
    use pars
    use vars
    
    implicit none

contains
!*******

subroutine read_pars(fln)
    character(len=*), intent(in) :: fln
	
	character(len=100) :: namein
    integer :: idfln

    namelist /orient_info/ theta_min, theta_max, theta_step
    namelist /mol_info/ PotI!, nini, mini
    namelist /laser_info/ Ncycle, PeakInt, omega0, phacep, add, Nt
    namelist /hhg_info/ sigma, domega, dat_dir, &
						isout_acc, isout_dip, isout_vec, isout_grd, isSpline3der
!-- 
	isout_dip = .false.
	isout_vec = .false.
	isout_acc = .true.
	isout_grd = .false.
	isSpline3der = .true.
	
	open(newunit = idfln, file = fln, status = 'old', action = 'read')
	    read(unit = idfln, nml = orient_info) 
	    rewind(unit = idfln)
	
	    read(unit = idfln, nml = mol_info)
	    rewind(unit = idfln)
	    
!~ 	    read(unit = idfln, nml = laser_info)
!~ 	    rewind(unit = idfln)

	    read(unit = idfln, nml = hhg_info)
	close(unit = idfln)

!-- read laser_para in 'dat_dir' folder !!! Aug 10 2017
	write(namein,'(A,I0,A)') trim(dat_dir) // '/', theta_min, '/' // 'para_laser.in'
	open(newunit = idfln, file = namein, status = 'old', action = 'read')
		read(unit = idfln, nml = laser_info)
	close(unit = idfln)
	
!-- convert to atomic unit
	
	Evec0 = dsqrt( PeakInt/auI )
!   phacep = phacep * Pi / 180.0_dp
    omega0 = omega0 / auE
    TimeFWHM = Ncycle * 2.0_dp * Pi / omega0 ! in a.u. time
    Ntmax = Ncycle * Nt
    TimeStep = TimeFWHM / dble(NtMax)
    PotI = PotI / auE
    
!-- calc NhhgMax = Ncutoff + 50, Npts_hhg = NhhgMax / Omega_step 
    
    Up = 0.25_dp * Evec0**2 / omega0**2
    
    NhhgMax  = ceiling( (3.17_dp*Up + PotI) / omega0 + 50 )

!-- print to screen
    
    write(*,'')
	write(*,*) '================= CALC HHG SPECTRA ================='
    write(*,'(A,I5)')'   Maximum harmonic order      = ', NhhgMax
    
!    write(*, '(A)', advance = 'no') '   Enter step of order: domega   = '
!    read (*,*) domega

    write(*,'(A,F6.2)') '   step of order               =  ', domega
    write(*,*) '===================================================='
    
    Npts_hhg = ceiling( NhhgMax / domega )
    
end subroutine
!-------------

end module
!*********
