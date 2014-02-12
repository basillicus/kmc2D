module moves
  implicit none
! attributes of the move
  type :: move_attr
     integer :: config ! Configuration of the molecule
     integer :: ini  ! global index of the initial site of the move
     integer :: fin  ! global index of the final site
     integer :: dir  ! direction of the move (from 1 to 4)
     real*8  :: barrier ! energy barrier for the move
     real*8  :: rate    ! rate of the move
     character*10 :: type ! type of the move
  end type move_attr
  type(move_attr), dimension(:), allocatable :: move
  integer No_moves,max_moves
  real*8 sum_rates
end module moves
