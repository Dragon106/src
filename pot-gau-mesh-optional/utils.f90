module utils

	implicit none 

contains 
!*******

subroutine make_filename( a1na2, a1, n, a2 ) 

	integer, intent(in) :: n
	character(len=*), intent(in)  :: a1,a2
	character(len=*), intent(out) :: a1na2
	character(len=5) :: ID

    write(ID,'(I5)') n ! write n into ID variable and change it from character type to integer type
    a1na2= trim(a1)//trim(adjustl(ID))//trim(a2)
 
end subroutine
!-------------

subroutine stop_error( msg )
	
	character*(*), intent(in) :: msg

	write (*,'(a)') msg
	stop 1

end subroutine 
!-------------

subroutine getwdir( wdir ) 

	character, intent(out) :: wdir*(*)

	call getarg(1, wdir)

	if( wdir == '' ) &
	    call stop_error("Enter working directory!")

end subroutine
!-------------

function getfln( workdir, filename )
	character*(*), intent(in) :: workdir, filename
	character :: getfln*( len_trim(workdir) + len_trim(filename) + 1 )

	getfln = trim(workdir) // '/' // trim(filename)

end function
!-----------

end module 
!*********
