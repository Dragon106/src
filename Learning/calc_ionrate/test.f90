program test
	use type_vars_m
	use calc_rate_m

	implicit none
	
	type(ioniz_t) :: ionrate
	
	call ionrate%read_para()
	call ionrate%read_ntmax()
	call ionrate%generate_vars()
	call ionrate%get_timestep()

end program
