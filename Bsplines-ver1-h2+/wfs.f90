module wfs 

    use mesh
    use basis
    use hamil
    
    implicit none

    integer :: Nmax    

    real(8) :: Emax

    real(8), allocatable :: en(:)
    real(8), allocatable :: wf(:,:,:)

contains 
!*******

subroutine Nmax_calc()
    integer :: i 
    
    do i = 1, NBAS
        if ( engs(i) > Emax ) exit
    enddo
    Nmax = i - 1
    
    allocate( en(Nmax) ) 
    
    en(:) = engs(1:Nmax)   
    
end subroutine
!-------------

subroutine calc_wf()
    
    integer :: ir, ie, n
    real(8), allocatable :: basrT(:,:), cn(:,:)
    
    allocate ( cn(Nbr,Nbe) )    
    allocate ( wf(Nrad, Neta, Nmax) )
    allocate ( basrT(Nbr, Nrad) )
    
    basrT = transpose(basr)

!$omp parallel do private(n, cn, ir, ie)
    do n = 1, Nmax 
    
        cn = reshape(coeffs(:,n), [Nbr,Nbe]) 
    
        forall(ir = 1:Nrad, ie = 1:Neta) 
            wf(ir,ie,n) = sum( cn * matmul(basrT(:,ir:ir), base(ie:ie,:) ) ) / rad(ir)
        end forall
    
    enddo 
!$omp end parallel do

write(*,'(a)') "WF DONE!"

end subroutine
!-------------

subroutine clr_wf()

    deallocate ( wf )

end subroutine
!-------------

subroutine calc_dip()
    
    real(8), allocatable :: dip(:), dipx(:), dipy(:), dipz(:)

    integer :: n, ir, ie
    
    real(8), allocatable :: rw12(:,:), rew12(:,:)
     
    allocate ( dip(Nmax), dipz(Nmax) )
    allocate ( rw12(Nrad, Neta), rew12(Nrad, Neta) )
    
    forall (ir = 1: Nrad, ie = 1: Neta)
        rw12 (ir,ie) = rad(ir) * w12(ir,ie)
        rew12(ir,ie) = rad(ir) * eta(ie) * w12(ir,ie)
    end forall
    
    forall (n = 1:Nmax)
        dip(n)  = sum( wf(:,:,n)**2 * rw12(:,:) ) 
        dipz(n) = sum (wf(:,:,n)**2 * rew12(:,:) ) 
    end forall
    
    deallocate (rw12, rew12)
!--
        
    open(unit=1,file='dipole.dat')
        write(1,'(i,2F)') (n, dip(n), dipz(n), n = 1, Nmax)
    close(1)
    
    write(*,'(A)') "=== DIPOLE INFO ==="
    write(*,*)     sum(dip), sum(dipz) 
    write(*,*) 
    
    deallocate (dip, dipz)
    
end subroutine
!-------------

subroutine wrt_20wf()

    integer :: i, j
    
    open(unit = 1, file = '20wf.in')

        do i = 1, Nrad        
            do j = 1, Neta
                write(1,'(*(E24.16))') rad(i), eta(j), wf(i,j,1:20)
            enddo
            write(1,*)
        enddo
        
    close(1)
    
end subroutine 
!-------------

subroutine prt_eigenv()

    integer :: i

    write (*, '(a)') '=== ENERGY INFO ==='

    write (*, '(a, f11.6)') 'Emax = ', Emax  
    write (*, '(a, i11  )') 'Nmax = ', Nmax
    write (*, '')
    
    if( engs(Nbas) < Emax ) then 
        write(*, '(a)') 'o_O You do not reach to Emax!'
        write(*, '')
    endif    

    write (*, '(a)') 'FIRST TEN ENERGY LEVELS'
    write (*, '(i2,3x,f11.6)') ( i, engs(i), i = 1, 10 )

    write (*, '')

end subroutine 
!-------------

subroutine wrt_eigenve()
    integer :: i, j, ni
    
    open(unit = 1, file = 'EM.in')
        write(1,'(I,E)') (ni, en(ni), ni = 1, Nmax)
    close(1)
    
!     open(unit = 1, file = 'WF.in')
    open(unit = 1, file = 'WF.in', form = 'unformatted') ! binary data
!     do ni = 1, Nmax
!         do i = 1, Nrad
!             do j = 1, Neta
!                write(1,'(E24.16)') wf(i,j,ni)
                write(1) wf
!             enddo
!         enddo
!         write(1,*)
!     enddo
    close(1)

end subroutine
!-------------

subroutine read_wf()
    integer :: i, j, n

!~     open(unit = 106, file = 'WF.in')
    open(unit = 106, file = 'WF.in', form = 'unformatted')
    
!~     do n = 1, Nmax
!~         do i = 1, Nrad
!~             do j = 1, Neta
!                read(106,*) wf(i,j,n)
                read(106) wf
!~             enddo
!~         enddo
!~         read(106,*)
!~     enddo
        
    close(106) 

end subroutine
!-------------

subroutine check_norm()
    
    integer :: i, j
    real(8) :: S, eps
        
    eps = 1.0d-10

open (unit = 1, file = 'norm.dat')
    
    do j = 1, Nmax
    do i = 1, j
        S = sum( wf(:,:,i) * wf(:,:,j) * w12 )
        
        if( i == j ) then
            if( abs(S-1.0d0) > eps ) then
                write(1,'(2I,X,2F)') i, j, S, abs(S-1.0d0)
            endif 
!                call stop_error('Wave function is not normalized') 
        else 
            if( abs(S) > eps ) then
                write(1,'(2I,X,2F)') i, j, S, abs(S)
            endif
!                call stop_error('Wave function is not orthogonal')
        endif
        
    enddo
    enddo

close(1)
    
end subroutine
!-------------

end module 
!*********
