module PtrArray_m
  use type_vars
  implicit none

  type :: CPtrArray1d
    complex(wp), pointer :: pn(:) 
  end type 

  type :: RPtrArray1d
    real(wp), pointer :: pn(:) 
  end type 

  type :: CPtrArray2d 
    complex(wp), pointer :: pn(:,:)
  end type
  
  type :: RPtrArray2d
    real(wp), pointer :: pn(:,:) 
  end type 

  interface PtrArray_associate
    module procedure RPtrArray1d_associate
    module procedure RPtrArray2d_associate 
	module procedure RPtrArray2d_associate_1
    module procedure CPtrArray1d_associate 
    module procedure CPtrArray2d_associate 
    module procedure CPtrArray2d_associate_1
  end interface 
  
  interface destroy_PtrArray
    module procedure destroy_RPtrArr1d
    module procedure destroy_RPtrArr2d
    module procedure destroy_CPtrArr1d
    module procedure destroy_CPtrArr2d
  end interface 
  
contains 
!*******

elemental subroutine destroy_RPtrArr1d(Ptr)
  type(RPtrArray1d), intent(inout) :: Ptr 

  Ptr%pn => null()
    
end subroutine 
!-------------

elemental subroutine destroy_CPtrArr1d(Ptr)
  type(CPtrArray1d), intent(inout) :: Ptr 

  Ptr%pn => null()
    
end subroutine 
!-------------

elemental subroutine destroy_RPtrArr2d(Ptr)
  type(RPtrArray2d), intent(inout) :: Ptr

  Ptr%pn => null()
  
end subroutine 
!-------------

elemental subroutine destroy_CPtrArr2d(Ptr)
  type(CPtrArray2d), intent(inout) :: Ptr

  Ptr%pn => null()
  
end subroutine 
!-------------

subroutine CPtrArray1d_associate(ptr, arr, blkSize)  
  complex(wp), intent(inout), target :: arr(:)
  integer, intent(in) :: blkSize(:) 
  type(CPtrArray1d), intent(inout) :: ptr(size(blkSize))
  
  integer :: nblock, offset
  integer :: ib
  
  nblock = size(blkSize)
  
  if( sum(blkSize) /= size(arr) ) then 
    stop 'CPtrArray1d_associate Error: Size mismatched!'
  endif 
  
  do ib = 1, nblock
    offset = sum( blkSize(1 : ib-1) ) 
  	ptr(ib)%pn => arr( offset + 1 : offset + blkSize(ib) ) 
  enddo

end subroutine 
!-------------

subroutine RPtrArray1d_associate(ptr, arr, blkSize)  
  real(wp), intent(inout), target :: arr(:)
  integer, intent(in) :: blkSize(:) 
  type(RPtrArray1d), intent(inout) :: ptr(size(blkSize))
  
  integer :: nblock, offset
  integer :: ib
  
  nblock = size(blkSize)
  
  if( sum(blkSize) /= size(arr) ) then 
    stop 'RPtrArray1d_associate Error: Size mismatched!'
  endif 
  
  do ib = 1, nblock
    offset = sum( blkSize(1 : ib-1) ) 
  	ptr(ib)%pn => arr( offset + 1 : offset + blkSize(ib) ) 
  enddo

end subroutine 
!-------------

subroutine CPtrArray2d_associate(ptr, arr, blkSize1, blkSize2) 
  complex(wp), intent(inout), target :: arr(:,:)
  integer, intent(in) :: blkSize1(:), blkSize2(:)
  type(CPtrArray2d), intent(inout) :: ptr(size(blkSize1), size(blkSize2))
  
  integer :: nblock1, nblock2
  integer :: offset1, offset2
  integer :: ib1, ib2

  if ( sum(blkSize1) /= size(arr,1) .or. & 
       sum(blkSize2) /= size(arr,2) ) then 
    stop 'CPtrArray2d_associate Error: Size mismatched!'
  endif 

  nblock1 = size(blkSize1)
  nblock2 = size(blkSize2)
  
  do ib1 = 1, nblock1
    offset1 = sum( blkSize1(1 : ib1-1) ) 
    
	do ib2 = 1, nblock2
	  offset2 = sum( blkSize2(1 : ib2-1) )
	  
	  Ptr(ib1, ib2)%pn => arr(offset1 + 1 : offset1 + blkSize1(ib1), &
	  						  offset2 + 1 : offset2 + blkSize2(ib2))
	enddo
  enddo
	
end subroutine 
!-------------

subroutine RPtrArray2d_associate(ptr, arr, blkSize1, blkSize2) 
  real(wp), intent(inout), target :: arr(:,:)
  integer, intent(in) :: blkSize1(:), blkSize2(:)
  type(RPtrArray2d), intent(inout) :: ptr(size(blkSize1), size(blkSize2))
  
  integer :: nblock1, nblock2
  integer :: offset1, offset2
  integer :: ib1, ib2

  if ( sum(blkSize1) /= size(arr,1) .or. & 
       sum(blkSize2) /= size(arr,2) ) then 
    stop 'CPtrArray2d_associate Error: Size mismatched!'
  endif 

  nblock1 = size(blkSize1)
  nblock2 = size(blkSize2)
  
  do ib1 = 1, nblock1
    offset1 = sum( blkSize1(1 : ib1-1) ) 
    
	do ib2 = 1, nblock2
	  offset2 = sum( blkSize2(1 : ib2-1) )
	  
	  Ptr(ib1, ib2)%pn => arr(offset1 + 1 : offset1 + blkSize1(ib1), &
	  						  offset2 + 1 : offset2 + blkSize2(ib2))
	enddo
  enddo
	
end subroutine 
!-------------

subroutine CPtrArray2d_associate_1(ptr, arr, blkSize)
  complex(wp), intent(inout), target :: arr(:,:)
  integer, intent(in) :: blkSize(:)
  type(CPtrArray2d), intent(inout) :: ptr(size(blkSize), size(blkSize) )
  
  call CPtrArray2d_associate(ptr, arr, blkSize, blkSize) 

end subroutine  
!-------------

subroutine RPtrArray2d_associate_1(ptr, arr, blkSize)
  real(wp), intent(inout), target :: arr(:,:)
  integer, intent(in) :: blkSize(:)
  type(RPtrArray2d), intent(inout) :: ptr(size(blkSize), size(blkSize) )
  
  call RPtrArray2d_associate(ptr, arr, blkSize, blkSize) 

end subroutine 
!-------------

end module 
