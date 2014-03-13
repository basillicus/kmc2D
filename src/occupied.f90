module occupied
  use info
  use moves
  use interactions
  implicit none
!======================================================================================
!... occ(i) - status of the site i=1,...,LL3 
!======================================================================================
!   occ =  0 - the site is empty
!   
!                 \     /        /   \      
!          o---x   o---x    o---x     o---x
!         /     \          /               \ 
!   occ =    1       2        3         4
!         ------CIS------|--------TRANS------
!                                           
!                                          
!               x--   --x     --x          x--
!              /       /        /         /
!             /       /        /         /
!            o--   --o        o--     --o 
!   occ =    5       6        7        8 
!                                           
!         ------CIS------|--------TRANS------
!                                          
!         x--      --x       --x        x--
!          \          \         \        \     
!           \          \         \        \   
!            o--      --o         o--    --o 
!   occ =    9       10        11       12
!  
!           o : Site used to describe the position of the molecule on the grid
!           x : Extreme of the molecule
!   
!======================================================================================

! occupancy of each site 
  integer, dimension(:), allocatable :: occ
  integer no_m
!
  integer, dimension(:), allocatable :: empty_info
!
  CONTAINS
 
  subroutine all_moves()
!..............................................................................
! all moves for all molecules are collected and numbered (by m)
!..............................................................................
    integer m,i,i1,j,inn,empty,oc,k,i2,i4,pivoting,isomerization
    real*8 barrier,rate
    real ran2
    logical  is_valid, oc1 
    if(testing) write(*,*) '<=================> MOVES <====================>'
    m=0
!
!==========================> Diffusion moves <=================
!
    DO i=1,LL3
       if(occ(i).gt.0) then
!     
!     
!                   (4)--------(3)--------(2)    
!                   ·.         · .         ·     
!                  ·  .       ·   .       ·      
!                 ·    .     ·     .     ·       
!                ·      .   ·       .   ·        
!               ·        · ·         . ·         
!             (5)--------(i)---------(1)          
!            · .         ·.          ·               
!           ·   .       ·  .        ·                
!          ·     .     ·    .      ·                 
!         ·       .   ·      .    ·                  
!        ·         . ·        .  ·                
!       (6)--------(7)---------(8)                  
!     
!       (j) : j-vector indicates the position of the neighbour site
!
          ! Horizontal configuration
          if (occ(i).lt.5) then
              inn = check_interactions(occ(i),i)
              barrier=elem_barriers(1) + inn*elem_barriers(2) 
              rate=exp(-beta*barrier+pref)
              ! Pivoting 1 and 2
              if(.not. is_occupied(neighb_M(i,8))) then   
                  pivoting=1
                  i1=neighb_M(i,8)
                  if (is_allowed(i,occ(i),pivoting)) then
                      call fill_in(m,occ(i),i,i1,pivoting,barrier,rate,'M diffusin')
                  endif
                  pivoting=2
                  if (is_allowed(i,occ(i),pivoting)) then
                      call fill_in(m,occ(i),i,i1,pivoting,barrier,rate,'M diffusin')
                  endif
              end if
              ! Pivoting 3
              pivoting=3
              i1=i
              if (is_allowed(i,occ(i),pivoting)) then
                  call fill_in(m,occ(i),i,i1,pivoting,barrier,rate,'M diffusin')
              endif
              if(.not. is_occupied(neighb_M(i,1))) then   
                  ! Pivoting 4 
                  pivoting=4
                  i1=neighb_M(i,1)
                  if (is_allowed(i,occ(i),pivoting)) then
                      call fill_in(m,occ(i),i,i1,pivoting,barrier,rate,'M diffusin')
                  end if
              end if
          ! Vertical-right case
          else if (occ(i).lt.9) then
             ! Compute the barriers
             inn = check_interactions(occ(i),i)
             barrier=elem_barriers(1) + inn*elem_barriers(2) 
             rate=exp(-beta*barrier+pref)
             ! Pivoting 1 
             if(.not. is_occupied(neighb_M(i,1))) then   
                 pivoting=1
                 i1=neighb_M(i,1)
                  if (is_allowed(i,occ(i),pivoting)) then
                     call fill_in(m,occ(i),i,i1,pivoting,barrier,rate,'M diffusin')
                 end if
             end if
             ! Pivoting 2
             pivoting=2
             i1=i
             if (is_allowed(i,occ(i),pivoting)) then
                 call fill_in(m,occ(i),i,i1,pivoting,barrier,rate,'M diffusin')
             end if
             ! Pivoting 3
             pivoting=3
             i1=i
             if (is_allowed(i,occ(i),pivoting)) then
                 call fill_in(m,occ(i),i,i1,pivoting,barrier,rate,'M diffusin')
             end if
             !Pivoting 4
             if(.not. is_occupied(neighb_M(i,4))) then   
                 pivoting=4
                 i1=neighb_M(i,4)
                 if (is_allowed(i,occ(i),pivoting)) then
                     call fill_in(m,occ(i),i,i1,pivoting,barrier,rate,'M diffusin')
                 end if
             end if
          ! Vertical-left case
          else if (occ(i).ge.9) then
             inn = check_interactions(occ(i),i)
             barrier=elem_barriers(1) + inn*elem_barriers(2) 
             rate=exp(-beta*barrier+pref)
             ! Pivoting 1
             if(.not. is_occupied(neighb_M(i,4))) then   
                 pivoting=1
                 i1=neighb_M(i,4)
                 if (is_allowed(i,occ(i),pivoting)) then
                     call fill_in(m,occ(i),i,i1,pivoting,barrier,rate,'M diffusin')
                 end if
             end if
             ! Pivoting 2
             pivoting=2
             i1=i
             if (is_allowed(i,occ(i),pivoting)) then
                 call fill_in(m,occ(i),i,i1,pivoting,barrier,rate,'M diffusin')
             end if
             ! Pivoting 3 and 4
             if(.not. is_occupied(neighb_M(i,5))) then   
                 pivoting=3
                 i1=neighb_M(i,5)
                 if (is_allowed(i,occ(i),pivoting)) then
                     call fill_in(m,occ(i),i,i1,pivoting,barrier,rate,'M diffusin')
                 end if
                 pivoting=4
                 if (is_allowed(i,occ(i),pivoting)) then
                     call fill_in(m,occ(i),i,i1,pivoting,barrier,rate,'M diffusin')
                 end if
             end if
          end if ! End of different configurations
      end if ! End of i is Occupied?
    END DO

!==========================> Isomerization <================================
!
!                     \                    /      
!       o---x   -->    o---x    or    o---x  
!      /     \              \        /       
!        
     do i=1,LL3
        if (occ(i) .gt. 0 ) then
            ! Isomerize the extreme "o"
            isomerization=1
            inn = check_inter_isomer(occ(i),i,isomerization)
            ! If the end is free,then check if is in a linker
            if (inn == 0 ) then
              ! Check if is in a linker and which kind it is.
              !------------------------------------------------------------------------------ 
              ! TODO:
              ! CHECK if is well called the subroutine
              !------------------------------------------------------------------------------ 
                 inn=check_inter_linker(occ(i),i,isomerization)
                 barrier = elem_barriers(3) + inn*elem_barriers(6)
            else
                 barrier = elem_barriers(3) + inn*elem_barriers(4)
            end if
            rate=exp(-beta*barrier+pref)
            call fill_in(m,occ(i),i,i,isomerization,barrier,rate,'M isomeriz')
            
            ! Isomerize the extreme "x"
            isomerization=2
            inn = check_inter_isomer(occ(i),i,isomerization)
            ! If the end is free,then check if is in a linker
            if (inn == 0 ) then
              ! Check if is in a linker and which kind is.
              !------------------------------------------------------------------------------ 
              ! TODO:
              ! CHECK if is well called the subroutine
              !------------------------------------------------------------------------------ 
                 inn=check_inter_linker(occ(i),i,isomerization)
                 barrier = elem_barriers(3) + inn*elem_barriers(6)
            else
                 barrier = elem_barriers(3) + inn*elem_barriers(4)
            end if
            rate=exp(-beta*barrier+pref)
            call fill_in(m,occ(i),i,i,isomerization,barrier,rate,'M isomeriz')
        end if
     end do
          
 

!==========================> Desorption <================================
!
    do i=1,LL3
       oc=occ(i)
       if(oc.ge.1) then
          inn = check_interactions(occ(i),i)
          barrier=elem_barriers(5) + inn*elem_barriers(2) 
          rate=exp(-beta*barrier+pref)
          call fill_in(m,occ(i),i,i,0,barrier,rate,'M desorptn')
      end if
   end do
!
!==============================> ADSORPTION <================================
!
   if(deposition) then
      empty=0
      do i=1,LL3
         if(occ(i).eq.0) then
            empty=empty+1   ;    empty_info(empty)=i  
         end if
      end do
      k= ran2(seed) * (empty+0.5) ; if(k.eq.0) k=1 ; if(k.gt.empty) k=empty
      i=empty_info(k) 
      call fill_in(m,occ(i),i,i,0,0.0d0,depos_rate,'M adsorbed')
   end if
!
!............. sum all legitimate rates
!
50 sum_rates=0.0d0
   do k=1,m
      sum_rates=sum_rates+move(k)%rate
   end do
!
!.................. FINISH counting moves
!
   No_moves=m
   if(testing) write(*,*)'..> Total number of moves found = ',No_moves

 end subroutine all_moves


! NOT NEEDED
!   integer function num_M_neighb(i)
!  !.... only the two sites along the In row are probed, above and below the molecule at site=i
!     integer i,i1,i0
!     i0=0
!     i1=neighb_M(i,2) ; if(occ(i1).gt.0) i0=i0+1
!     i1=neighb_M(i,4) ; if(occ(i1).gt.0) i0=i0+1
!     num_M_neighb=i0
!   end function num_M_neighb

 subroutine fill_in(m,config,i,i1,j,barrier,rate,typ)
!..............................................................................................
! fills in the array move for each possible step m with all necessary attributes
!..............................................................................................
   integer m,i,i1,j, config
   real*8 barrier,rate
   character*10 typ
   m=m+1
   if(m.gt.max_moves) stop 'max # of moves reached'
   move(m)=move_attr(config,i,i1,j,barrier,rate,typ)
   if(testing) write(*,'(i5,x,a10,x,3(i5,x),2(i3,x),f10.5,x,e12.6)') m,typ,config,i,i1,j,0,barrier,rate
 end subroutine fill_in

 logical function  is_occupied(i)
    implicit none
    integer :: i
    ! Check if the site i is occupied 
    is_occupied = .true.  
    ! if (allow_high_coverage) then
        if (occ(i) == 0 ) then
             is_occupied = .false.
        end if
    ! else ! Do not allow side-by-side molecules (Low Coverage Phase)
    !     if (occ(i) == 0 ) then
    !        if  (occ(neighb_M(i,8)).lt.9) then 
    !             if (occ(neighb_M(i,7)).gt.8 .or. occ(neighb_M(i,7)).lt.5) then
    !                 if (occ(neighb_M(i,5)).eq.0 .or. occ(neighb_M(i,5)).gt.4) then
    !                    is_occupied = .false.
    !                 end if
    !             end if
    !        end if
    !     end if
    ! end if
 end function is_occupied

 logical function is_allowed (i,config,pivoting)
 implicit none
 integer  :: i,config,pivoting,nc
 !--------------------------------------------------------------------------------
 ! i                : Position 'i' of the lattice (LL3)
 ! config           : Config of molecule on position 'i' (1..12)
 ! pivoting         : Kind of pivoting (1,2,3,4)
 ! nc (neighb. conf): Conformation of the neighbour of 'i' (1..12)
 !--------------------------------------------------------------------------------
 !--------------------------------------------------------------------------------
 ! If (allow_high_coverage): allows that difussion yields this kind of
 ! configurations:
 !
 !        \_____/ /
 !                \
 !                 \                  
 !                 /                  
 !
 ! else : Diffusion only can yield this kind of configurations:
 !
 !        \_____/  \_____/ 
 !
 !--------------------------------------------------------------------------------
                             

if (allow_high_coverage) then 
     is_allowed  = .true.
 else 
     is_allowed = .false.
     if (config >= 1 .and. config <= 4 ) then ! Horizontal case
        if (pivoting == 1 ) then   ! Pivoting #1
            ! Logic for checking the environment (Bit complex...Sorry!)
            nc = occ(neighb_M(i,8))
            if ( nc == 0) then
                nc = occ(neighb_M(i,1)) ! Looking at this new position
                if ( nc == 0 .or. nc >= 5 .and. nc <= 8 ) then
                    nc = occ(neighb_M(i,7))
                    if ( nc == 0 .or.  nc >= 5 ) then
                        nc = occ(neighb_M(neighb_M(i,8),1)) ! Looking for Second Neighbours 
                        if ( nc <= 8 )  then
                            nc = occ(neighb_M(neighb_M(i,8),8)) 
                            if ( nc <= 8 )  then
                                is_allowed = .true. 
                            end if
                        end if
                    end if
                end if 
            end if ! End logic
        else if (pivoting == 2) then  ! Pivoting #2
            ! Logic of checking the environment
            nc = occ(neighb_M(i,8)) 
            if (nc == 0) then
                nc = occ(neighb_M(i,5)) 
                if (nc == 0 .or. nc >= 5 ) then
                    nc = occ(neighb_M(i,7)) 
                    if ( nc == 0 .or. nc >= 9 ) then
                        nc = occ(neighb_M(neighb_M(i,8),7))   
                        if ( nc <= 4 .or. nc >= 9 )  then
                            is_allowed = .true. 
                        end if
                    end if
                end if 
            end if ! End logic
        else if ( pivoting == 3 ) then ! Pivoting #3
            ! Logic of checking the environment
            nc = occ(neighb_M(i,1))
            if ( nc <= 8 ) then
                nc = occ(neighb_M(i,3)) ! Looking at this new position
                if ( nc == 0 .or. nc >= 5 .and. nc <= 8 ) then
                    nc = occ(neighb_M(i,4))
                    if ( nc == 0 .or.  nc >= 5 ) then
                        nc = occ(neighb_M(i,5))
                        if ( nc == 0 .or.  nc >= 5 ) then
                            nc = occ(neighb_M(i,8))
                            if ( nc <= 8 )  then
                                is_allowed = .true. 
                            end if
                        end if
                    end if
                end if 
            end if ! End logic
        else if ( pivoting == 4 ) then ! Pivoting #4
            ! Logic of checking the environment
            nc = occ(neighb_M(i,1)) 
            if (nc == 0) then
                nc = occ(neighb_M(i,3)) 
                if (nc == 0 .or. nc >= 9 ) then
                    nc = occ(neighb_M(i,4)) 
                    if ( nc == 0 .or. nc >= 5 ) then
                        nc = occ(neighb_M(i,8)) 
                        if ( nc <= 4 .or. nc >= 9 ) then
                            is_allowed = .true. 
                        end if
                    end if
                end if 
            end if ! End logic
        end if
     else if (config >= 5 .and. config <= 8 ) then ! Vertical-right case
        if (pivoting == 1 ) then   ! Pivoting #1
            ! Logic of checking the environment
            nc = occ(neighb_M(i,1)) 
            if (nc == 0) then
                nc = occ(neighb_M(i,3)) 
                if (nc == 0 .or. nc >= 9 ) then
                    nc = occ(neighb_M(i,4)) 
                    if ( nc == 0 .or. nc >= 5 ) then
                        nc = occ(neighb_M(i,8)) 
                        if ( nc <= 4 .or. nc >= 9 ) then
                            is_allowed = .true. 
                        end if
                    end if
                end if 
            end if ! End logic
        else if ( pivoting == 2 ) then 
            ! Logic of checking the environment
            nc = occ(neighb_M(i,1))
            if ( nc <= 4 ) then
                nc = occ(neighb_M(i,7)) ! Looking at this new position
                if ( nc <= 4 .or. nc >= 9) then
                    nc = occ(neighb_M(i,8))
                    if ( nc <= 4 ) then
                        nc = occ(neighb_M(neighb_M(i,8),1)) ! Looking for Second Neighbours 
                        if ( nc <= 8 )  then
                            is_allowed = .true. 
                        end if
                    end if
                end if 
            end if ! End logic
        else if ( pivoting == 3 ) then 
            ! Logic of checking the environment
            nc = occ(neighb_M(i,4)) 
            if (nc == 0 .or. nc >= 9 ) then
                nc = occ(neighb_M(neighb_M(i,4),5)) 
                if (nc == 0 .or. nc >= 5 ) then
                    nc = occ(neighb_M(i,5)) 
                    if ( nc == 0 .or. nc >= 9 ) then
                        nc = occ(neighb_M(i,7))   
                        if ( nc <= 4 .or. nc >= 9 )  then
                            is_allowed = .true. 
                        end if
                    end if
                end if 
            end if ! End logic
        else if ( pivoting == 4 ) then 
            ! Logic of checking the environment
            nc = occ(neighb_M(i,1))
            if ( nc <= 8 ) then
                nc = occ(neighb_M(i,3))
                if ( nc <= 4 ) then
                    nc = occ(neighb_M(i,4)) ! Looking at this new position
                    if ( nc == 0 ) then
                        nc = occ(neighb_M(i,5))
                        if ( nc <= 4 .or. nc >= 9 ) then
                            is_allowed = .true. 
                        end if
                    end if
                end if 
            end if ! End logic
        end if 
     else if (config >= 9 ) then ! Vertical-left case
        if (pivoting == 1 ) then   ! Pivoting #1
            ! Logic of checking the environment
            nc = occ(neighb_M(i,1))
            if ( nc <= 8 ) then
                nc = occ(neighb_M(i,3))
                if ( nc <= 4 ) then
                    nc = occ(neighb_M(i,4)) ! Looking at this new position
                    if ( nc == 0 ) then
                        nc = occ(neighb_M(i,5))
                        if ( nc <= 4 .or. nc >= 9 ) then
                            is_allowed = .true. 
                        end if
                    end if
                end if 
            end if ! End logic
        else if ( pivoting == 2 ) then 
            ! Logic of checking the environment
            nc = occ(neighb_M(i,1))
            if ( nc <= 8 ) then
                nc = occ(neighb_M(i,3)) ! Looking at this new position
                if ( nc == 0 .or. nc >= 5 .and. nc <= 8 ) then
                    nc = occ(neighb_M(i,4))
                    if ( nc == 0 .or.  nc >= 5 ) then
                        nc = occ(neighb_M(i,5))
                        if ( nc == 0 .or.  nc >= 5 ) then
                            nc = occ(neighb_M(i,8))
                            if ( nc <= 8 )  then
                                is_allowed = .true. 
                            end if
                        end if
                    end if
                end if 
            end if ! End logic
        else if ( pivoting == 3 ) then 
            ! Logic of checking the environment
            nc = occ(neighb_M(i,5))
            if ( nc == 0 ) then
                nc = occ(neighb_M(i,6)) ! Looking at this new position
                if ( nc <= 4 .or. nc >= 9) then
                    nc = occ(neighb_M(i,7))
                    if ( nc <= 4 ) then
                        nc = occ(neighb_M(i,8)) 
                        if ( nc <= 4 .or. nc >= 9 )  then
                            is_allowed = .true. 
                        end if
                    end if
                end if 
            end if ! End logic
        else if ( pivoting == 4 ) then 
            ! Logic of checking the environment
            nc = occ(neighb_M(i,4))
            if ( nc == 0 .or. nc >= 5 .and. nc <= 8 ) then
                nc = occ(neighb_M(i,5)) ! Looking at this new position
                if ( nc == 0 ) then
                    nc = occ(neighb_M(neighb_M(i,5),5))
                    if ( nc == 0 .or.  nc >= 5 ) then
                        nc = occ(neighb_M(i,7))
                        if ( nc <= 8 ) then
                            is_allowed = .true. 
                        end if
                    end if
                end if 
            end if ! End logic
        end if ! End pivotings
    end if ! End three main cases
  end if ! End (allow_high_coverage)

 end function

 integer function check_interactions(A,pos)
 !---------------------------------------------------------------
 ! Input:  A - occ of the molecule to be checked
 !         pos - site number of that molecule
 ! Output: number of all H bonds it created with its neighbours
 !---------------------------------------------------------------
     use info
     implicit none
     integer  :: A,B,j,pos,i,inn
     
 ! Check the interactions for diffusion and desorption movements
     inn = 0
     if (A <= 4) then
         ! Horizontal case
         do i=1,4
             j = direction(i)
             B = occ(neighb_M (pos,j))
             if ( B > 0 .and. B <= 4 ) then 
                 inn = inn + M_int(A,B,i)
             end if
         end do
     else if ( A <= 8 ) then
         ! Vertical-right case
         do i=5,8
             j = direction(i)
             B = occ(neighb_M(pos,j)) - 4
             if ( B > 0 .and. B <= 4 ) then 
                 inn = inn + M_int(A-4,B,i-4)
             end if
         end do
     else ! if ( A <= 12 ) then
         ! Vertical-left case
         do i=9,12
             j = direction(i)
             B = occ(neighb_M(pos,j)) - 8
             if ( B > 0 .and. B <= 4 ) then 
                 inn = inn + M_int(A-8,B,i-8)
             end if
         end do
     end if
     check_interactions = inn
 end function
 
 integer function check_inter_isomer (A,pos,kind)
     use info
     ! Kind of isomerization (1 or 2)
     implicit none
     integer  :: A,B,pos,kind,j,inn
 
     ! Check the environment for the isomeriztaion movement
     inn = 0
     if (A <= 4) then 
         j = dir_isomer(A,kind)
         B = occ(neighb_M(pos,j))
     else if (A <= 8) then 
         j = dir_isomer(A,kind)
         B = occ(neighb_M(pos,j)) - 4
     else ! if (A <= 12 ) then 
         j = dir_isomer(A,kind)
         B = occ(neighb_M(pos,j)) - 8
     end if
     ! Molecule A interacts with molecules B with the same direction as 
     ! A (parallel molecules) regardless of the configuration of B.
     if ( B > 0 .and. B <= 4 ) inn = 1
     check_inter_isomer = inn
 end function
 
 integer function check_inter_linker (A,pos,kind)
     use info
     ! Kind of isomerization (1 or 2)
     implicit none
     integer  :: A,B,C,pos,kind,inn
 
     ! Check if the molecule is in a linker, and if yes
     ! which kind of linker
     inn = 0
     ! --------------------------------------------------------------------------------
     ! TODO: Create a better routine. Up to now we only increase the barrier if
     ! the monomer to isomerize is trans, regardless of the environment. The
     ! increase of the barrier if the monomer is part of a linker is too low
     ! (0.020 eV) so in principle has a very little effect, and we can consider
     ! the isomerization of the "tail of a linker" as an isomerization of an
     ! isolated monomer.
     ! --------------------------------------------------------------------------------
     if (A == 3  .or. A == 4 .or.  A == 7  .or. A == 8 .or.  A == 11 .or. A == 12 ) then
         inn = 1
     end if
 end function
end module occupied


