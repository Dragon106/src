module utils

	implicit none

contains
!*******

subroutine make_filename(a1na2,a1,n,a2)
	integer, intent(in) :: n
	character(len=*), intent(in)  :: a1,a2
	character(len=*), intent(out) :: a1na2
	character(len=5) :: ID

    write(ID,'(I5)') n 	! write n into ID var and change it from character type to integer one   
    a1na2= trim(a1)//trim(adjustl(ID))//trim(a2)
    
end subroutine
!-------------

subroutine stop_error( msg )
	character*(*), intent(in) :: msg

	write (*,'(a)') msg
	stop 1

end subroutine 
!-------------

end module
!*********
