  integer function nn_site(i,j)
    use info
    implicit none
!... i - current site global index
!    j - given direction to the nn (from 1 to 4 in the anticlockwise order)
! nn_site - gives the global number of the neighbour corresponding to i=>j
    integer i,j,i1
    real*8 r(2)
    r(:)=x_lat(i,:)+Xnn(j,:) 
    call find_site(r,i1)
    nn_site=i1
  end function nn_site

  subroutine find_site(r,i)
    use info
    implicit none
!........ on the whole lattice
!   i - gives the global index
    real*8 r(2),x(2)
    integer i1,i
    logical ask_equiv
    do i1=1,LL3
       i=i1  ; x(:)=x_lat(i1,:)-r(:)
       if(ask_equiv(x,BL)) return
    end do
    i=0
  end subroutine find_site

  logical function ask_equiv(r,XJ)
!........................................................................
!    It checks whether vector r() is equivalent to lattice translations
!  reciprocal to vectors XJ (up to 2*pi):
!
!  if XJ - reciprocal lattice vectors, then r() must be in the direct space;
!  if XJ - direct lattice vectors, then r() must be in the reciprocal space;
!........................................................................
!  .true. - if r() is equivalent to latice translations
!  .false.- if not.
!........................................................................
    use param
    implicit none
    real*8 r(2),XJ(2,2),x
    integer i
    ask_equiv=.false.
!
!____ multiply on all direct lattice vectors: the result must be integer
! if r() is an arbitrary sum of reciprocal lattice vectors (note, r() is
! defined up to 2*pi because of the choice of XJ used to generate it)
!
    do i=1,2
       x=r(1)*XJ(i,1)+r(2)*XJ(i,2)
       if( abs(x-nint(x)) .gt. tiny ) return
    end do
    ask_equiv=.true.
  end function ask_equiv

 subroutine posit(char,len,len00)
! the space in char is filled from the right; this finds the first
! non-empty character from the left, lev00
   character*(*)  char
   integer len,i,j,len00
   do i=len,1,-1
      j=i
      if(char(i:i).eq.' ') then
         len00=j+1
         return
      end if
   end do
   len00=1
 end subroutine posit

 logical function if_in_latt(i1,i2)
   use info
   implicit none
   integer i1,i2,i
   do i=1,LL3
      if(n1(i).eq.i1.and.n2(i).eq.i2) then
         if_in_latt=.true.
         return
      end if
   end do
   if_in_latt=.false.
 end function if_in_latt
