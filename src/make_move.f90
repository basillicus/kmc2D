  subroutine which_move(m)
    use param
    use moves
!.........................................................................
! From all the moves possible it take one using KMC rules:
! m    - move number to be taken
! time - time for that move to take place (in units of the rates)
!.........................................................................
    implicit none
    integer m,m1
    real*8 R,x,rx,x1
    real ran2
!
!__________________ calculate time
    R=sum_rates ; x1=ran2(seed) ; dt=-1.0/R*log(x1)  
!
!_____________ trap the transition number
!    
    x=ran2(seed) ; rx=0.0
    do m1=1,No_moves
       rx=rx+move(m1)%rate/R 
       if(x.lt.rx) then
          m=m1 ; return
       end if
    end do
  end subroutine which_move

  subroutine perform_move(m,occ,kmc)
    use info
    use moves
    use change
    use interactions
    implicit none
    integer m,i,j,i1,kmc,conf
    integer occ(0:LL3)
    character*10 type
!........... get the attributes of the successful move m
    conf=move(m)%config; i=move(m)%ini ; i1=move(m)%fin 
    ! j = direction of the pivoting
    j=move(m)%dir ; type=move(m)%type 

    select case (type)
!
!................. diffusion
!
    case ('M diffusin')
        ! Logic of the diffusion
        occ(i) = 0
        ! Computes the final configuration depending on the 
        ! direction of the pivoting
        !
        ! Horizontal configuration
        select case (conf)
        case (1) 
            if (j.eq.1 .or. j.eq.3) then
                  occ(i1) = 5
            else if (j.eq.2 .or. j.eq.4) then
                  occ(i1) = 10 
            end if
        case(2)
            if (j.eq.1 .or. j.eq.3) then
                 occ(i1) = 6
            else if (j.eq.2 .or. j.eq.4) then
                 occ(i1) = 9 
            end if
        case(3)
            if (j.eq.1 .or. j.eq.3) then
                 occ(i1) = 7 
            else if (j.eq.2 .or. j.eq.4) then
                 occ(i1) = 11 
            end if
        case(4)
            if (j.eq.1 .or. j.eq.3) then
                 occ(i1) = 8 
            else if (j.eq.2 .or. j.eq.4) then
                 occ(i1) = 12 
            end if
            ! Vertical-Right Configuration
        case(5 )
            if (j.eq.1 .or. j.eq.3) then
                 occ(i1) = 9 
            else if (j.eq.2 .or. j.eq.4) then
                 occ(i1) = 1 
            end if
        case(6 )
            if (j.eq.1 .or. j.eq.3) then
                 occ(i1) = 10 
            else if (j.eq.2 .or. j.eq.4) then
                 occ(i1) = 2 
            end if
        case(7 )
            if (j.eq.1 .or. j.eq.3) then
                 occ(i1) = 11 
            else if (j.eq.2 .or. j.eq.4) then
                 occ(i1) = 3 
            end if
        case(8 )
            if (j.eq.1 .or. j.eq.3) then
                 occ(i1) = 12 
            else if (j.eq.2 .or. j.eq.4) then
                 occ(i1) = 4 
            end if
            ! Vertical-Left Configuration
        case (9 )
            if (j.eq.1 .or. j.eq.3) then
                 occ(i1) = 2 
            else if (j.eq.2 .or. j.eq.4) then
                 occ(i1) = 5 
            end if
        case (10 )
            if (j.eq.1 .or. j.eq.3) then
                 occ(i1) = 1 
            else if (j.eq.2 .or. j.eq.4) then
                 occ(i1) = 6 
            end if
        case (11) 
            if (j.eq.1 .or. j.eq.3) then
                 occ(i1) = 3 
            else if (j.eq.2 .or. j.eq.4) then
                 occ(i1) = 7 
            end if
        case (12)
            if (j.eq.1 .or. j.eq.3) then
                 occ(i1) = 4 
            else if (j.eq.2 .or. j.eq.4) then
                 occ(i1) = 8 
            end if
        end select
            
       ! occ(i1)=1 ; occ(i)=0 
       i_chg(1)=i  ; occ_chg(1)=occ(i) ; i_chg(2)=i1 ; occ_chg(2)=occ(i1) ; no_chg=2 
!
!....... Isomerization
!
    case ('M isomeriz')
        if (conf <= 4 ) then
            occ(i) = isomer_M (conf,j)
        else if (conf <= 8 ) then
            occ(i) = isomer_M (conf-4 ,j) + 4
        else  ! if (conf <= 12 )
            occ(i) = isomer_M (conf-8 ,j) + 8
        end if

!................. desorption
!
    case ('M desorptn')
       occ(i)=0
       no_chg=1 ; i_chg(1)=i; occ_chg(1)=0
!
!................. adsorption
!
    case ('M adsorbed')
       occ(i)=int(ran(seed)*12)+1
       no_chg=1 ; i_chg(1)=i; occ_chg(1)=occ(i)
    end select

  end subroutine perform_move

!   subroutine move_hor_to_vert (kind_of_pivoting, new_conf, fin )
! 
!      implicit none
!      integer :: kind_of_pivoting, new_conf, fin
! 
!      if  ( kind_of_pivoting == 1 ) then
!          occ(fin) = new_conf
! 
! 
! 
