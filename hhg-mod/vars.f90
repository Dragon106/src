module vars
	use pars
	
	implicit none
	
	integer :: Ncycle, Nt, Ntmax
	real(dp) :: add
	
	integer :: NhhgMax, Npts_hhg
	integer :: theta_min, theta_max, theta_step
	real(dp) :: theta
	
	integer :: nini, mini
	real(dp) :: PotI
	
	real(dp) :: TimeFWHM, TimeStep, PeakInt, Evec0, omega0, phacep
	real(dp) :: Up
	
	real(dp), allocatable :: order(:)
	real(dp)  :: sigma, domega
	
	logical :: isout_dip, isout_vec, isout_acc, isout_grd
	logical :: isSpline3der
	character :: dat_dir*100
	
	type :: arr_xy !my_type
	  integer 			    :: lb = 0, &
							   ub = 0 
	  real(dp), allocatable :: x(:), y(:)
	end type
	
	type(arr_xy) :: acc_dip, dip
	type(arr_xy) :: dz, d2z
	type(arr_xy) :: hhg, hhgs, hhgodd !!!!! 11Apr
	type(arr_xy) :: phase

!### lession 1 : pattern reconige  ###

contains
!*******

subroutine alloc_arr_xy(arr, b1, b2)
  type(arr_xy), intent(inout) :: arr
  integer,   intent(in)  :: b1, b2
  optional :: b2
  integer  :: lb, ub
  
  if ( present(b2) ) then
    lb = b1 
    ub = b2
  else
    lb = 1
    ub = b1
  endif
  
  arr%lb = lb 
  arr%ub = ub 

  if( allocated(arr%x) ) deallocate( arr%x )
  allocate( arr%x(lb:ub) )
  
  if( allocated(arr%y) ) deallocate( arr%y )
  allocate( arr%y(lb:ub) )
  
end subroutine
!-------------

subroutine resize_arr_xy( arr, b1, b2 ) 
  type(arr_xy), intent(inout) :: arr
  integer,   intent(in)  :: b1, b2
  optional :: b2
  integer  :: lb, ub
  
  type(arr_xy) :: arr_tmp

  if ( present(b2) ) then
    lb = b1 
    ub = b2
  else
    lb = 1
    ub = b1
  endif

  call alloc_arr_xy( arr_tmp, arr%lb, arr%ub )
  
  arr_tmp%x = arr%x
  arr_tmp%y = arr%y
  
  call dealloc_arr_xy( arr ) 
  call alloc_arr_xy( arr, b1, b2 ) 
 !!! not general 
 
  lb = max(arr%lb, arr_tmp%lb); ub = min(arr%ub, arr_tmp%ub)
  
  arr%x(lb:ub) = arr_tmp%x(lb:ub) 
  arr%y(lb:ub) = arr_tmp%y(lb:ub) 

  call dealloc_arr_xy( arr_tmp )

end subroutine 

subroutine dealloc_arr_xy(arr)
    type(arr_xy), intent(inout) :: arr

	arr%lb = 0 
	arr%ub = 0 

    if ( allocated(arr%x) ) deallocate( arr%x )
    if ( allocated(arr%y) ) deallocate( arr%y )

end subroutine
!-------------

subroutine cre_vars()

    allocate( order(Npts_hhg) )
        order(:) = dble((/1:Npts_hhg/)) * NhhgMax/Npts_hhg !!!

    call alloc_arr_xy( hhg,    Npts_hhg )
    call alloc_arr_xy( hhgs,   Npts_hhg )
    call alloc_arr_xy( phase,  Npts_hhg )
    call alloc_arr_xy( hhgodd, NhhgMax )
    
end subroutine
!-------------

subroutine clr_vars()

	if(allocated(order)) deallocate( order )
	
	call dealloc_arr_xy( acc_dip )
	call dealloc_arr_xy( dip )
	call dealloc_arr_xy( hhg )
	call dealloc_arr_xy( hhgs )
	call dealloc_arr_xy( phase )
    call dealloc_arr_xy( hhgodd )

end subroutine
!-------------

end module
!*********
