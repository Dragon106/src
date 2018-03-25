module laser_m
	use type_vars
	use consts
	implicit none
	
	private

    type :: elaser_t
		integer  :: itype_env

 		integer  :: Ncycle
		integer  :: Ncycle_turns, Ncycle_const
		integer  :: Nt, Ntmax

 		real(wp) :: PeakInt, omega0, phacep, Add
 		real(wp) :: Evec0, TimeFWHM, TimeStep

 		real(wp) :: theta, cstheta, sntheta
 		
 		contains
 		
 		procedure :: read_params => read_laser_para
 		procedure :: setup => setup_laser
 		procedure :: elaser => elaser
 		procedure :: write => wrt_laser_info
 		procedure :: destroy => destroy_laser	
    end type
    
    public :: elaser_t
	
contains 
!*******

subroutine read_laser_para(this, inp_file)
	class(elaser_t),  intent(inout) :: this
	character(len=*), intent(in)    :: inp_file

	integer  :: itype_env
	integer  :: Ncycle, Ncycle_turns, Ncycle_const
	integer  :: Nt
	real(wp) :: PeakInt, omega0, phacep, Add
	
	integer :: idf

	namelist /laser_info/ itype_env, Ncycle, Ncycle_turns, Ncycle_const, &
			              PeakInt, omega0, phacep, Add, Nt

	open(newunit = idf, file = inp_file, status = 'old', action = 'read')
		read(unit = idf, nml = laser_info)
	close(unit = idf)
	
	call destroy_laser(this)

	this%itype_env    = itype_env
	this%Ncycle       = Ncycle
	this%Ncycle_turns = Ncycle_turns
	this%Ncycle_const = Ncycle_const
	this%Nt           = Nt
	this%PeakInt      = PeakInt
	this%omega0       = omega0
	this%phacep       = phacep
	this%Add	      = Add

!-- convert to atomic unit
	
	this%Evec0    = sqrt( this%PeakInt / auI )
	this%phacep   = this%phacep * Pi / 180.0_wp
    this%omega0   = this%omega0 / auE
    this%TimeFWHM = this%Ncycle * 2.0_wp * Pi / this%omega0 ! in a.u. time
    this%Ntmax    = this%Ncycle * this%Nt
    this%TimeStep = this%TimeFWHM / real(this%NtMax,wp)
	
	call prt_laser_info(this)
	
end subroutine 
!-------------

subroutine prt_laser_info(this)
	class(elaser_t), intent(in) :: this
	
	write(*,'(A)')        		 "=== LASER PARAMETERS ==="
	write(*,'(A,I0,A)')   		 'Number of cycles  =   ', this%Ncycle, ' cycles'
	write(*,'(A,F6.2,A,F8.3,A)') 'Time FWHM         = ', this%Ncycle *2.0_wp*pi/this%omega0 * auT, ' fs             = ', &
														 this%Ncycle *2.0_wp*pi/this%omega0 , ' a.u.T'
	write(*,'(A,F6.2,A,F8.3,A)') 'Peak of intensity = ', this%PeakInt, ' x 10^14 W/cm^2 = ', this%Evec0, ' a.u.I'
	write(*,'(A,F6.2,A,F8.3,A)') 'Energy            = ', this%omega0*auE, ' eV             = ', this%omega0, ' a.u.E'
	write(*,'(A,F6.2,A)') 		 'Ponderomotive     = ', 0.25_wp * this%Evec0**2 / this%omega0**2, ' a.u.E'
	write(*,'(A)')        		 ""
	write(*,*) this%itype_env, this%Ncycle_turns, this%Ncycle_const

end subroutine
!-------------

subroutine wrt_laser_info(this,outdir)
	class(elaser_t),  intent(in) :: this
	character(len=*), intent(in) :: outdir
	
	character (len=100) :: nameout
	
	integer :: idf
	
	nameout = trim(outdir)// '/' // 'para_laser.in'
	
	open(newunit = idf, file = nameout, status = 'unknown', action = 'write')
		write(idf,'(A)')         '&LASER_INFO'
		write(idf,'(4x,A,I0)')   'Ncycle  = ', this%Ncycle
		write(idf,'(4x,A,F8.3)') 'PeakInt = ', this%PeakInt
		write(idf,'(4x,A,F8.3)') 'omega0  = ', this%omega0 * auE
		write(idf,'(4x,A,F8.3)') 'phacep  = ', this%phacep * rad2deg
		write(idf,'(4x,A,F8.3)') 'Add     = ', this%Add
		write(idf,'(4x,A,I0)')   'Nt      = ', this%Nt
		write(idf,'(A)')         '/'
	close(unit = idf)
	
end subroutine
!-------------

subroutine setup_laser(this, theta) 
  	class(elaser_t), intent(inout) :: this
	real(wp),        intent(in)    :: theta
	
    this%theta   = theta*deg2rad 
    this%cstheta = cos(this%theta)
    this%sntheta = sin(this%theta)

end subroutine 
!-------------

subroutine write_laser(this, fln)
    class(elaser_t),  intent(in) :: this 
	character(len=*), intent(in) :: fln
	
	real(wp) :: t, Ey, Ez
	integer :: idf
	integer :: i
	
	open(newunit = idf, file = fln, status = 'unknown', action = 'write')

	do i = 0, this%Ntmax
		t = i * this%TimeStep
		call this%elaser(t, Ey, Ez)
		write(idf,'(3E16.8)') t, Ez, Ey
	enddo
	
	close(idf)

end subroutine 
!------------- 

subroutine elaser(this, t, Ey, Ez)
	class(elaser_t), intent(in)  :: this
	real(wp),        intent(in)  :: t
	real(wp),        intent(out) :: Ey, Ez
	
	real(wp) :: Et

!~ 	if (this%itype_env == 1) call elaser_sinsquared()
!~ 	if (this%itype_env == 2) call elaser_trapezoid()

	select case (this%itype_env)
	case(1)
		call elaser_sinsquared()
	case(2)
		call elaser_trapezoid()
	case default
		write(*,'(A)') 'Please choose 1 (sin-squared) or 2 (trapezoid)!'
		stop
	end select 

!~ 	select case (this%itype_env)
!~ 	case("SinSquared") 
!~ 		call elaser_sinsquared()
!~ 	case("Trapezoid")
!~ 		call elaser_trapezoid()
!~ 	end select 

contains
!-------

	subroutine elaser_sinsquared()
		
		Et = this%Evec0 * sin(pi*t/this%TimeFWHM)**2 * sin(this%omega0*t + this%phacep)
		
		Ey = Et * this%sntheta
		Ez = Et * this%cstheta
		
	end subroutine
	!-------------

	subroutine elaser_trapezoid()

		real(wp) :: t_on, t_off

		t_on  = this%Ncycle_turns * 2*pi / this%omega0
		t_off = (this%Ncycle_const + this%Ncycle_turns) * 2*pi / this%omega0

		if (t <= t_on) then
			Et = t / t_on * cos(this%omega0 * t + this%phacep) &
			   + 1.0_wp / t_on * sin(this%omega0 * t + this%phacep) / this%omega0
		elseif( t > t_on .and. t < t_off ) then 
			Et = cos(this%omega0 * t + this%phacep )
		else 
			Et = (this%TimeFWHM - t) / (this%TimeFWHM - t_off) * cos(this%omega0 * t + this%phacep) &
			   - 1.0_wp / (this%TimeFWHM - t_off) * sin(this%omega0 * t + this%phacep) / this%omega0 
		endif
			Et = - this%Evec0 * Et
	
		Ey = Et * this%sntheta
		Ez = Et * this%cstheta

	end subroutine
	!-------------

end subroutine
!-------------

subroutine destroy_laser(this)
	class(elaser_t), intent(inout) :: this
	
	this%Ncycle  = 0
	this%Nt      = 0
	this%PeakInt = 0.0_wp
	this%omega0  = 0.0_wp
	this%phacep  = 0.0_wp
	this%add     = 0.0_wp
	
end subroutine
!-------------

end module
