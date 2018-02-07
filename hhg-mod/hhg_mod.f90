module hhg_mod
    use pars
    use vars
    use integ_weigs_m
    implicit none

    interface calc_hhg
        module procedure calc_hhg
        module procedure calc_hhg_arr_xy
    end interface

contains  
!*******

subroutine calc_hhg_arr_xy(ho, dip, hhg, hhgodd, hhgs, phase) !!!
    real(dp),     intent(in)    :: ho(:)
    type(arr_xy), intent(in)    :: dip
    type(arr_xy), intent(inout) :: hhg, hhgodd, hhgs, phase
  
    call calc_hhg(ho, dip%x, hhg%x, hhgodd%x, hhgs%x, phase%x)
    call calc_hhg(ho, dip%y, hhg%y, hhgodd%y, hhgs%y, phase%y)

end subroutine 
!-------------

subroutine calc_hhg(ho, dip, hhg, hhgodd, hhgs, phase)
	real(dp), intent(in)  :: ho(:), dip(0:)
	real(dp), intent(out) :: hhgodd(:)
	real(dp), intent(out) :: hhg(size(ho)), hhgs(size(hhg)), phase(size(hhg))
	
	real(dp) :: OmegaMax, omega!, omega_t
	real(dp) :: signal_real, signal_imag
	
	integer :: i
	
!~ 	write(*,*) lbound(dip), ubound(dip)
	
	OmegaMax = Pi/TimeStep
	
	if ( NhhgMax*omega0 >= OmegaMax ) then
       write(*,*) 'Please increase the value of Nt!'
       Stop
    endif

block
    integer :: ndip
    real(8) :: time(size(dip)), wt(size(dip))
    
    ndip = size(dip)
    time = [0:ndip-1]*TimeStep
    
    call integ_weigs(time, wt)
    
    do i = 1, size(hhg)
		omega = ho(i) * omega0
		
!		signal_real = 0.0_dp
!		signal_imag = 0.0_dp
 		
!		do j = 0, size(dip)
!--			omega_t = order(i)*omega0 * dble(j)*TimeStep
!			omega_t = omega * dble(j)*TimeStep			
!			signal_real = signal_real + dip(j) * cos(omega_t)
!			signal_imag = signal_imag + dip(j) * sin(omega_t)
!		enddo
		
		signal_real = sum( wt(:) * dip(:) * cos( time * omega ) )
		signal_imag = sum( wt(:) * dip(:) * sin( time * omega ) )
		
!		hhg(i) = TimeStep*TimeStep * abs(sum( dip(0:ndip) * exp( time * omega ) ) )
!		hhg(i) = TimeStep**2 * abs(signal)
		
!~ 		hhg(i) = TimeStep**2 * (signal_real**2 + signal_imag**2)
		hhg(i) = signal_real**2 + signal_imag**2
		if (hhg(i) == 0.0_dp) hhg(i) = tiny(0.0_dp)
		
		phase(i) = atan2(signal_imag, signal_real) * 180.0_dp/Pi
    
    enddo

end block

    call calc_odd_hhg(ho, hhg, hhgodd)
    call calc_hhg_smooth(hhg, hhgs)

end subroutine
!-------------

!** Smooth spectra procedure follow M. Lein article **
! S(omega') = \int { [S(omega) * exp { - (omega - omega')^2} / sigma^2] d(omega) } 
! sigma = n * omega0, n = 1, 2, 3...

subroutine calc_hhg_smooth(hhg, hhg_smooth)
	real(dp), intent(in)  :: hhg(:)
	real(dp), intent(out) :: hhg_smooth(size(hhg))

	real(dp)  :: omega, omega_p, alpha, omega_step, phase_temp
	
	integer :: i, j
	
    alpha = sigma * omega0
	omega_step = dble(NhhgMax) / Npts_hhg
	
	do i = 1, size(hhg)
		hhg_smooth(i) = 0.0_dp
		omega = dble(i) * omega_step * omega0
		do j = 1, size(hhg)
			omega_p = dble(j) * omega_step * omega0
			phase_temp = - 1.0_dp * (omega_p - omega)**2 / alpha**2
			hhg_smooth(i) = hhg_smooth(i) + hhg(j) * exp(phase_temp)
		enddo	
		hhg_smooth(i) = hhg_smooth(i) * omega_step * omega0
		if (hhg_smooth(i) == 0.0_dp) hhg_smooth(i) = 1.0d-6
	enddo

end subroutine
!-------------

subroutine calc_odd_hhg(ho, hhg, hhgo)
    real(8), intent(in)  :: ho(:), hhg(:)
    real(8), intent(out) :: hhgo(:)
    
    real(8) :: delta
    integer :: Nhhg
    integer :: n
        
    Nhhg = nint( ho(size(ho)) )
!    write(*,*) NhhgMax
    
    delta = ho(2) - ho(1)

block
    logical :: mask(1:size(ho))
      
    do n = 1, Nhhg, 2
    
        mask = ho >= dble(n-1) .and. ho < dble(n+1)
        hhgo(n) = sum( hhg, mask = mask )
       
    enddo
    hhgo = hhgo * delta

end block
    
end subroutine
!-------------

subroutine wrt_hhg(outdir, fln1, fln2, ho, hhg, hhgo, hhgs, iangle)
    character(len=*), intent(in) :: outdir, fln1, fln2
    integer,   intent(in) :: iangle
    real(dp),  intent(in) :: ho(:)
    type(arr_xy), intent(in) :: hhg, hhgo, hhgs
    
    character :: nameout*100
    integer   :: idfln
    integer   :: k
    
!~      write(nameout, '(a,i0,a)') trim(outdir) // '/' // trim(fln1) // '-', iangle, '.dat'
    write(nameout, '(a,i0,a)')  trim(outdir) // '/', iangle, '/' // trim(fln1) // '.dat' !!! 
    
!    print*, nameout, size(ho), size(hhg%x), size(hhg%y), size(hhgs%x), size(hhgs%y) 
    
    open(newunit = idfln, file = nameout)
    	write(idfln, '(5(E16.8,X))') ( ho(k), hhg%x(k), hhgs%x(k), & 
									          hhg%y(k), hhgs%y(k), k = 1, size(ho) )

    close(unit = idfln)

!--    
!~     write(nameout, '(a,i0,a)') trim(outdir) // '/' // trim(fln2) // '-', iangle, '.dat'
    write(nameout, '(a,i0,a)') trim(outdir) // '/', iangle, '/' // trim(fln2) // '.dat' !!! 
    
    open(newunit = idfln, file = nameout)
    	write(idfln, '(I, 4E16.8)') ( k, hhgo%x(k), log10(hhgo%x(k)), & 
									     hhgo%y(k), log10(hhgo%y(k)), k = 1, size(hhgo%x), 2 )

    close(unit = idfln)
    
end subroutine
!-------------

subroutine wrt_pha(outdir, fln, ho, phase, iangle)
    character(len=*), intent(in) :: outdir, fln
    integer,   intent(in) :: iangle
    real(dp),  intent(in) :: ho(:)
    type(arr_xy), intent(in) :: phase
    
    character :: nameout*100
    integer   :: idfln
    integer   :: k
    
!~     write(nameout, '(a,i0,a)') trim(outdir) // '/' // trim(fln)// '-', iangle, '.dat'
    write(nameout, '(a,i0,a)') trim(outdir) // '/', iangle, '/' // trim(fln) // '.dat' !!! for './iangle/fln.dat'

    open(newunit = idfln, file = nameout)
    	write(idfln, '(3(E16.8,X))') ( ho(k), phase%x(k), phase%y(k), k = 1, size(ho) )
    close(unit = idfln)
    
end subroutine
!-------------

end module
