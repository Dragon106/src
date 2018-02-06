module utils
    
    implicit none 

contains
!*******
subroutine stop_error(msg)
    character, intent(in) :: msg*(*)

    write(*,'(a)') msg
    stop 1 

end subroutine 
!-------------

function getfln(workdir, filename)
character*(*), intent(in) :: workdir, filename
character :: getfln*( len_trim(workdir) + len_trim(filename) + 1 )

getfln = trim(workdir) // '/' // trim(filename)

end function
!-----------

subroutine makedir( dir )
character, intent(in) :: dir*(*) 

call system( "mkdir -p " // dir )

end subroutine 
!-------------

end module 
!*********
