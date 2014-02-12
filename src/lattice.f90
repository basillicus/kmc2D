module lattice
  implicit none
!... direct and reciprocal lattice vectors: the format is this ->  (#,xy)
  real*8 :: AJ(2,2),BJ(2,2)
  real*8 :: pi=3.1415926535898

CONTAINS

  subroutine setup()
!................................................................................
!___ lattice vectors 
    AJ(1,1)=1.0            ;      AJ(1,2)=0.0
    AJ(2,1)=cos(pi/3.0)    ;      AJ(2,2)=sin(pi/3.0)

!___ Construct the reciprocal lattice vectors
    call get_reciprocal_vectors (AJ,BJ)

    write(9,*)'...> Direct lattice vectors: '
    write(9,'(a,2(x,f10.5))') '[1] ',AJ(1,:)
    write(9,'(a,2(x,f10.5))') '[2] ',AJ(2,:)
    write(9,*)'...> Reciprocal lattice vectors: '
    write(9,'(a,2(x,f10.5))') '[1] ',BJ(1,:)
    write(9,'(a,2(x,f10.5))') '[2] ',BJ(2,:)

  end subroutine setup

subroutine get_reciprocal_vectors (A,B)
    implicit none
    real*8 :: A(2,2)
    real*8 :: B(2,2)
    real*8     :: denom
    ! A_i*B_j = delta_(ij)
    ! Using Cramer's Rule for solving the system
    denom = A(1,1)*A(2,2) - A(1,2)*A(2,1)
    B(1,1) =  A(2,2)/denom
    B(1,2) = -A(2,1)/denom
    B(2,1) = -A(1,2)/denom
    B(2,2) =  A(1,1)/denom
end subroutine get_reciprocal_vectors

end module lattice
