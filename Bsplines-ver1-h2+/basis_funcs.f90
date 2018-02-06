
!** Associated Legendre polynomials in recurrence form **
!   (l-m) * P(l,m) = x*(2l-1) * P(l-1,m) + (l+m-1) * P(l-2,m) => P(l,m)
!   P(m,m)   = (-1)**m * (2m-1)!! * (1 - x**2) ** (m/2)
!   P(m+1,m) = x * (2m+1) * P(m,m)

subroutine cre_nplg(ql, qm, x, plg)

    integer, intent(in) :: ql, qm
    real(8), intent(in) :: x(:)
    real(8), intent(out) :: plg(size(x))
    
    real(8) :: pmm(size(x)), pmmp1(size(x)), pll(size(x))
    
    integer :: i, ll, mi
    
!    if (m .lt. 0 .or. m .gt. l .or. abs(xi) .gt. 1._8) pause 'Bad arguments!'
!    if (mi .lt. 0 .or. mi .gt. li) stop 'Bad arguments!'

    mi = abs(qm) 

    if( mi > ql ) then
        plg = 0.0d0
        return
    endif
    
    

!-
    pmm = 1.d0 ! P00 = 1
    
!- General    

    if (mi .gt. 0) then
        pmm = (-1.d0)**mi * fact2(2*mi-1) * (1.d0 - x**2)**(mi/2.d0)
    endif
    
!-- if (l = m) then plg = P(m,m)    
    if (ql .eq. mi) then; plg = pmm 
    
!-- if (l != m)
    else
        pmmp1 = x * (2*mi+1) * pmm
        
        if (ql .eq. (mi+1)) then
            plg = pmmp1
        else
            do ll = mi+2, ql
                pll = ( x * (2*ll-1) * pmmp1 - (ll+mi-1) * pmm ) / (ll - mi) 
                pmm = pmmp1                                                
                pmmp1 = pll                                                
            enddo
            plg = pll
        endif        
    endif

    if( qm < 0 ) then
        plg = (-1)**mi * plg 
    endif

    plg = plg * sqrt( (2*ql + 1.d0)/2 * fact(ql-mi) / fact(ql+mi) )     
    
end subroutine
!--

!function deriv_plg(li,mi,xi)
!    integer, intent(in) :: li, mi
!    real(8), intent(in) :: xi
    
!    real(8) :: deriv_plg
    
!    deriv_plg = (li * xi * plg(li,mi,xi) - (li+mi) * plg(li-1,mi,xi) )&
!                / sqrt( 1.d0 - xi**2 )

!end function
!--
function fact(n)
! calculate factorial
    integer, intent(in) :: n
    real(8) :: fact 
    
    fact = gamma(n + 1.0d0)

end function
!-----------

function fact2(n) 
! calculate double factorial  
    integer, intent(in) :: n
    real(8) :: fact2
    
    integer :: i
    
    fact2 = 1.d0
    do i = 1, n, 2
        fact2 = fact2*i
    enddo
    
end function 
!-----------

subroutine cre_bspline(bps, deg, pts, bspl, dbspl)
    integer, intent(in) :: deg
    real(8), intent(in) :: bps(:), pts(:) 
    real(8), allocatable, intent(out) :: bspl(:,:), dbspl(:,:)
    
    real(8), allocatable :: bt0(:,:), bt1(:,:)
    real(8), allocatable :: knots(:)
    
    integer :: nbps, npts, nbs_max, nbs, nknots
    real(8) :: temp, eps
    integer :: i, ix, k
        
    nbps = size(bps)
    npts = size(pts) 
    nbs_max = nbps + deg - 2
    nbs = nbs_max
    nknots = nbps + 2*deg - 2

    if( allocated(bspl) .or. allocated (dbspl) ) & 
        call stop_error('Bspline basis are allocated') 
    
!    allocate( bspl(npts,nbs), dbspl(npts,nbs) ) 
    allocate( bspl(npts,nbs-1), dbspl(npts,nbs-1) ) !!!!! Tutu
    allocate( knots(nknots) )
    allocate( bt0(npts,nknots-1), bt1(npts,nknots-1) )
    
    knots(1 : deg) = bps(1)
    knots(deg+1 : nbps+deg-2) = bps(2 : nbps-1)
    knots(nbps+deg-1 : nknots) = bps(nbps)

    eps = epsilon(0.d0)

!-- calc: k=1
    do i = 1, nknots-1 
    
        where ( pts(:) .GE. knots(i) .AND. pts(:) .LT. knots(i+1) )
            bt1(:,i) = 1.d0
        else where
            bt1(:,i) = 0.d0
        end where
            
    enddo
        
!-- calc: Bspline of order k=2 -> deg
! in order to get dspline of order k-th, need bspline order k-1
    do k = 2, deg 
    
        bt0(:, 1 : nknots-k+1) = bt1(:, 1 : nknots-k+1)
                    
        do i = 1, nknots-k
            temp = knots(i+k-1) - knots(i)  
            if ( temp .eq. 0.d0 ) temp = eps                     
            bt1(:,i) = ( pts(:) - knots(i) ) / temp * bt0(:,i) 

            temp = knots(i+k) - knots(i+1) 
            if ( temp .eq. 0.d0 ) temp = eps 
            bt1(:,i) = bt1(:,i) + ( knots(i+k) - pts(:) ) / temp * bt0(:,i+1) 
        enddo
        
    enddo

!    bspl(:,:) = bt1(:,1:nbs)
!~     bspl(:,:) = bt1(:,1:nbs-1)
    bspl(:,:) = bt1(:,2:nbs-1)
        
!-- calc: derivative of Bspline of order k=deg   
    do i = 1, nbs
        temp = knots(i+deg-1) - knots(i)
        if (temp .eq. 0.d0) temp = eps
        bt1(:,i) = ( bt0(:,i) / temp )
        
        temp = knots(i+deg) - knots(i+1)
        if (temp .eq. 0.d0) temp = eps
        bt1(:,i) = (deg-1) * ( bt1(:,i) - bt0(:,i+1) / temp )            
    enddo

!    dbspl(:,:) = bt1(:,1:nbs)
!~     dbspl(:,:) = bt1(:,1:nbs-1)
    dbspl(:,:) = bt1(:,2:nbs-1)
        
    deallocate( knots )
    deallocate( bt0, bt1 )
 
end subroutine

