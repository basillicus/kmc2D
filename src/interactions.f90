module interactions
! Definition of the matrix of interactions
integer :: M_int(4,4,4)

! ---- Used for Difussion ---------
!      i  : Current position of the molecule
!      j  : j-value of the position of the other molecule the current one may interact with
!           (is parallel to the current one) 
! j=direction(i) : it points to the position of another molecule in 4 possible 
!                   j-directions for each of the 3 cases of i:
!    i=1-4  -> "j-vectors" to check for the horizontal case
!    i=5-8  -> "j-vectors" to check for the vertical-right case
!    i=9-12 -> "j-vectors" to check for the vertical-left case
! Depending on the initial state of the molecule, you have to check different
! sites (i => direction j) around the site i for other molecules which it can interact with. 
integer  :: direction(12)=(/3,4,7,8,4,5,8,1,5,7,1,3/)

! ---- Used for Isomerization ---------
! isomer_M(i,kind) : Final state (occ) of the molecule (after isomerization)
!     i: Initial state of the molecule (i=1...4)
!     kind: Which end of the molecule is to be isomerized:
!          kind = 1 -> isomerize the end "o"
!          kind = 2 -> isomerize the end "x"
!     Other states are converted by symetry
!               5 => 1
!               6 => 2
!               . => .
!               . => .
!              12 => 4
integer  :: isomer_M(4,2)=reshape((/4,3,2,1,3,4,1,2/),(/4,2/))

! dir_isomer(i,kind) : "j" vector to check whether there is any interaction
!                      with the molecule of the same group in this direction 
!      i: the initial state of A (i=1...12)
!      kind: the kind of the isomerization (kind=1,2)
integer  :: dir_isomer(12,2) = reshape &
((/7,4,7,4,8,5,8,5,1,7,1,7,8,3,3,8,1,4,4,1,3,5,5,3/),(/12,2/))

CONTAINS

subroutine init_interactions(ignore_single_HB)
implicit none

logical ignore_single_HB
! Fill the values of the matrix of interactions
!  Scheeme of the meaning of the matrix M_int(A,B,k):
  
if ( ignore_single_HB ) then 
    
! B-->   o---x             o---x        ! Symbol "#" represents interactions 
!       /   # \           /             ! as hydrogen bonds
!           #  #         #  
!           \  #  /     o---x      
! A-->       o---x     /     \  
!        k = 4         k = 3
!        2 inter.      0 inter.         ! 0 interactions casuse is an unstable
!                                       ! strucure on GasPhase
       
! 
!          A = 1              A = 2                  A = 3                A = 4         
!B = 1     A,B,k
    M_int (1,1,1) = 0;   M_int (2,1,1) = 2;   M_int (3,1,1) = 2;   M_int (4,1,1) = 0
    M_int (1,1,2) = 0;   M_int (2,1,2) = 2;   M_int (3,1,2) = 0;   M_int (4,1,2) = 2
    M_int (1,1,3) = 0;   M_int (2,1,3) = 0;   M_int (3,1,3) = 0;   M_int (4,1,3) = 0
    M_int (1,1,4) = 0;   M_int (2,1,4) = 0;   M_int (3,1,4) = 0;   M_int (4,1,4) = 0
! B = 2
    M_int (1,2,1) = 0;   M_int (2,2,1) = 0;   M_int (3,2,1) = 0;   M_int (4,2,1) = 0
    M_int (1,2,2) = 0;   M_int (2,2,2) = 0;   M_int (3,2,2) = 0;   M_int (4,2,2) = 0
    M_int (1,2,3) = 2;   M_int (2,2,3) = 0;   M_int (3,2,3) = 2;   M_int (4,2,3) = 0
    M_int (1,2,4) = 2;   M_int (2,2,4) = 0;   M_int (3,2,4) = 0;   M_int (4,2,4) = 2
! B = 3
    M_int (1,3,1) = 0;   M_int (2,3,1) = 2;   M_int (3,3,1) = 2;   M_int (4,3,1) = 0
    M_int (1,3,2) = 0;   M_int (2,3,2) = 0;   M_int (3,3,2) = 0;   M_int (4,3,2) = 0
    M_int (1,3,3) = 2;   M_int (2,3,3) = 0;   M_int (3,3,3) = 2;   M_int (4,3,3) = 0
    M_int (1,3,4) = 0;   M_int (2,3,4) = 0;   M_int (3,3,4) = 0;   M_int (4,3,4) = 0
! B = 4
    M_int (1,4,1) = 0;   M_int (2,4,1) = 0;   M_int (3,4,1) = 0;   M_int (4,4,1) = 0
    M_int (1,4,2) = 0;   M_int (2,4,2) = 2;   M_int (3,4,2) = 0;   M_int (4,4,2) = 2
    M_int (1,4,3) = 0;   M_int (2,4,3) = 0;   M_int (3,4,3) = 0;   M_int (4,4,3) = 0
    M_int (1,4,4) = 2;   M_int (2,4,4) = 0;   M_int (3,4,4) = 0;   M_int (4,4,4) = 2

else  ! Simulate the real system
    
!                               /                
! B-->   o---x             o---x        ! Symbol "#" represents interactions 
!       /   # \           /             ! as hydrogen bonds
!           #  #         #  
!           \  #  /     o---x      
! A-->       o---x     /     \  
!        k = 4         k = 3
!        2 inter.      1 inter.         ! 1 interactions        
! 
! A: Central molecule state, but only the 1st four state are needed as others are
!    obtained by subtracting 4 or 8
! B: State of the neighbour of A it can interact with (only the first 4 are needed,
!    as others are obtained by subtracting 4 or 8
! k: number of the parallel neighbour (k=1...4); j is taken from direction(i) and
!    obtained by subtracting 0, 4 or 8 depending on A=occ(i).
! M_int(A,B,k) = number of H bonds between A and B in the "direction" k
!      
!
!          A = 1              A = 2                  A = 3                A = 4         
!B = 1     A,B,k
    M_int (1,1,1) = 1;   M_int (2,1,1) = 2;   M_int (3,1,1) = 2;   M_int (4,1,1) = 1
    M_int (1,1,2) = 1;   M_int (2,1,2) = 2;   M_int (3,1,2) = 1;   M_int (4,1,2) = 2
    M_int (1,1,3) = 1;   M_int (2,1,3) = 0;   M_int (3,1,3) = 1;   M_int (4,1,3) = 0
    M_int (1,1,4) = 1;   M_int (2,1,4) = 0;   M_int (3,1,4) = 0;   M_int (4,1,4) = 1
! B = 2
    M_int (1,2,1) = 0;   M_int (2,2,1) = 1;   M_int (3,2,1) = 1;   M_int (4,2,1) = 0
    M_int (1,2,2) = 0;   M_int (2,2,2) = 1;   M_int (3,2,2) = 0;   M_int (4,2,2) = 1
    M_int (1,2,3) = 2;   M_int (2,2,3) = 1;   M_int (3,2,3) = 2;   M_int (4,2,3) = 1
    M_int (1,2,4) = 2;   M_int (2,2,4) = 1;   M_int (3,2,4) = 1;   M_int (4,2,4) = 2
! B = 3
    M_int (1,3,1) = 1;   M_int (2,3,1) = 2;   M_int (3,3,1) = 2;   M_int (4,3,1) = 1
    M_int (1,3,2) = 0;   M_int (2,3,2) = 1;   M_int (3,3,2) = 0;   M_int (4,3,2) = 1
    M_int (1,3,3) = 2;   M_int (2,3,3) = 1;   M_int (3,3,3) = 2;   M_int (4,3,3) = 1
    M_int (1,3,4) = 1;   M_int (2,3,4) = 0;   M_int (3,3,4) = 0;   M_int (4,3,4) = 1
! B = 4
    M_int (1,4,1) = 0;   M_int (2,4,1) = 1;   M_int (3,4,1) = 1;   M_int (4,4,1) = 0
    M_int (1,4,2) = 1;   M_int (2,4,2) = 2;   M_int (3,4,2) = 1;   M_int (4,4,2) = 2
    M_int (1,4,3) = 1;   M_int (2,4,3) = 0;   M_int (3,4,3) = 1;   M_int (4,4,3) = 0
    M_int (1,4,4) = 2;   M_int (2,4,4) = 1;   M_int (3,4,4) = 1;   M_int (4,4,4) = 2

end if

end subroutine

end module
