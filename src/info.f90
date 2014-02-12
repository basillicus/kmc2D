module info
  use param
  use lattice
  implicit none
!... lattice sites info
! i - counts all lattice sites (global)
! x_lat(i),n1(i),n2(i) - give access to the sites and their indices
!                        via the global index i
! n1(i),n2(i) - indeces in expanding the position i in lattice vectors
! Xnn(j,2) - vectors to 8 neighbours (j=1,..,8) counted anticlockwise

  real*8 Xnn(8,2)
  real*8,dimension(:,:),allocatable :: x_lat,x_lat_ext
  integer,dimension(:),allocatable :: n1,n2
  integer,dimension(:,:),allocatable :: neighb_M
!... direct and reciprocal space lattice vectors (for PBC)
  real*8 AL(2,2),BL(2,2)
!....... occ for printing the map of sites
  integer,dimension(:),allocatable :: occ1

  CONTAINS

!=====================================================================================
!================================== T H E   P B C    C A S E =========================
!=====================================================================================

 subroutine build_latt_pbc()
    use lattice
    integer i,k,i1,i2,j
    integer nn_site
!
!............. for periodic boundary conditions (PBC): 
!              (i) construct direct super-lattice vectors AL(#,xy)
!              (ii) construct reciprocal super-lattice vectors BL(#,xy)
!
    AL(:,:)=L*AJ(:,:)  ;  BL(:,:)=BJ(:,:)/L
    if(iPrnt.gt.4) then
       write(9,*)'...> Direct super-lattice vectors: '
       write(*,*)'...> Direct super-lattice vectors: '
       do i=1,2
          write(9,'(a,i1,a,2(x,f10.5))') '[',i,'] ',AL(i,:)
          write(*,'(a,i1,a,2(x,f10.5))') '[',i,'] ',AL(i,:)
       end do
       write(9,*)'...> Reciprocal super-lattice vectors: '
       do i=1,2
          write(9,'(a,i1,a,2(x,f10.5))') '[',i,'] ',BL(i,:)
       end do
    end if
!
!............... vectors to nearest neighbours out of any site 
!
 !    Xnn(1,1)= 1.0   ;  Xnn(1,2)= 0.0 
 !    Xnn(2,1)= 1.0   ;  Xnn(2,2)= 1.0  
 !    Xnn(3,1)= 0.0   ;  Xnn(3,2)= 1.0
 !    Xnn(4,1)=-1.0   ;  Xnn(4,2)= 1.0
 !    Xnn(5,1)=-1.0   ;  Xnn(5,2)= 0.0
 !    Xnn(6,1)=-1.0   ;  Xnn(6,2)=-1.0
 !    Xnn(7,1)= 0.0   ;  Xnn(7,2)=-1.0
 !    Xnn(8,1)= 1.0   ;  Xnn(8,2)=-1.0
     Xnn(1,1)= AJ(1,1)        ;  Xnn(1,2)= 0.0 
     Xnn(2,1)= AJ(1,1)+AJ(2,1) ;  Xnn(2,2)=AJ(2,2)   
     Xnn(3,1)= AJ(2,1)        ;  Xnn(3,2)= AJ(2,2)
     Xnn(4,1)=-AJ(1,1)+AJ(2,1) ;  Xnn(4,2)= AJ(2,2)
     Xnn(5,1)=-AJ(1,1)        ;  Xnn(5,2)= 0.0
     Xnn(6,1)=-AJ(1,1)-AJ(2,1) ;  Xnn(6,2)=-AJ(2,2)
     Xnn(7,1)=-AJ(2,1)       ;  Xnn(7,2)=-AJ(2,2)
     Xnn(8,1)= AJ(1,1)-AJ(2,1) ;  Xnn(8,2)=-AJ(2,2)

!............. build an extended lattice of sites for the drawing
!
    k=0
    do i1=-1,L+1
       do i2=-1,L+1
          k=k+1  ; x_lat_ext(k,:)=i1*AJ(1,:)+i2*AJ(2,:)
          if(iPrnt.eq.5) write(9,*) k,i1,i2,x_lat_ext(k,:)
       end do
    end do
!
!............. build the lattice of grid points 
!
    i=0 
    do i1=0,L-1
       do i2=0,L-1
          i=i+1 ; x_lat(i,:)=i1*AJ(1,:)+i2*AJ(2,:) ; n1(i)=i1 ; n2(i)=i2  
       end do
    end do
!______________________ print
    if(iPrnt.eq.4) then
       do i=1,LL3
          write(9,'(3(i5,x),2(f10.5,x))') i,n1(i),n2(i),x_lat(i,:)
       end do
    end if
!
!.............. give a map of the lattice with site numbers
!
!    allocate(occ1(0:LL3)) ; call draw(0,occ1,.true.,0.0d0) ;  deallocate(occ1)
!
!.............. info arrays for the nearest neighbours sites
!
    do i=1,LL3
       do j=1,8
          neighb_M(i,j)=nn_site(i,j)
       end do
       If(iPrnt.ge.4) write(9,'(a,i5,a,4(x,i5),a)')  &
          'M site ',i,' has as its neighbours sites: ',(neighb_M(i,j),j=1,8) &
          ,'  along dir = 1-8'
    end do
   
    write(*,*)' ... finished with neighb_M ...'
!
!........... write all info into the file, unit=11
    write(11) AL,BL,x_lat_ext,x_lat,n1,n2,neighb_M

  end subroutine build_latt_pbc

end module info
