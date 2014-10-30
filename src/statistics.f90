module statistics
    use param
    use moves
    
    implicit none
    integer :: nCis,nLtrans,nDtrans,nTotal
    integer :: nHbond,nNeighb  ! Environment analisys
    integer, save :: mDiff_free=0,mIsomer_free=0,mAdsorb=0,mDesor=0 ! Free moves
    integer, save :: mDiff_cluster=0,mIsomer_cluster=0 ! Moves from clusters
    integer, parameter :: lust=34    ! Logic Unit for statistcs

CONTAINS

subroutine init_statistics () 
    implicit none

    if (restart) then 
        open(lust,file='statistics.kmc',position='append')
    else
        open(lust,file='statistics.kmc',status='new')
        write(lust,'(a)') '# Statistics of KMC'
        write(lust,'(a)') '#1  : Time(micro secs) ' ! time
        write(lust,'(a)') '#2  : T (K) '            ! temp
        write(lust,'(a)') '#3  : Cis  '             ! nCis
        write(lust,'(a)') '#4  : L-trans '          ! nLtrans
        write(lust,'(a)') '#5  : D-Trans '          ! nDtrans
        write(lust,'(a)') '#6  : Total Molecs '     ! nCis + nLtrans + nDtrans
        write(lust,'(a)') '#7  : Diffusions free  ' ! mDiff_free
        write(lust,'(a)') '#8  : Diffusions from a cluster '   ! mDiff_cluster
        write(lust,'(a)') '#9  : Isomerizations free'          ! mIsomer_free
        write(lust,'(a)') '#10 : Isomerizations double H bond' !  mIsomer_cluster
        write(lust,'(a)') '#11 : N-Hbond  '         ! nHbond
        write(lust,'(a)') '#12 : N-Neighbors'       ! nNeighb
    end if

end subroutine init_statistics  

subroutine perform_statistics (m,kmc,occ,time,temp)
    implicit none
    integer, intent(in) :: m,kmc
    integer, intent(in) :: occ(0:LL3) 
    real*8 , intent(in) :: time, temp
    real*8 :: barrier
    integer  :: counter=0, i
    character*10  :: type
    
    type=move(m)%type ; barrier=move(m)%barrier 

    select case (type)
    case ('M diffusin')
        if (barrier > elem_barriers(1) ) then
            mDiff_cluster = mDiff_cluster + 1
        else
            mDiff_free = mDiff_free + 1
        end if
    case ('M isomeriz')
        if (barrier > elem_barriers(3) ) then
            mIsomer_cluster = mIsomer_cluster + 1
        else
            mIsomer_free = mIsomer_free + 1
        end if
    case ('M desorptn')
        mDesor = mDesor + 1
    case ('M adsorbed')
        mAdsorb = mAdsorb + 1
    end select

    if (counter >= freq_statistics ) then
        nHbond=0; nNeighb=0
        do i=1,LL3
        ! TODO: Include Count the Cis and Trans
            if (occ(i) > 0 ) call check_environment (occ(i),i,nHbond,nNeighb)
        end do
        nTotal = nCis + nLtrans + nDtrans
        nHbond = nHbond/2  ! Because A interacts with B and B interacts with A
        write(lust,'(f18.6,f8.2,4(i8),4(i11),2(i8))') time*1d-6, temp, nCis, nLtrans, nDtrans, & 
            nTotal, mDiff_free,mDiff_cluster, mIsomer_free,mIsomer_cluster, nHbond, nNeighb 
        counter=0
    end if
    counter=counter+1

end subroutine perform_statistics 

subroutine check_environment(A,pos,nNeighb,nHbond)
!---------------------------------------------------------------
! Input:  A - occ of the molecule to be checked
!         pos - site number of that molecule
! Output: number of all H bonds it created with its neighbours
!         number of neighbours, regardless are forming an Hbond
!---------------------------------------------------------------
    !use info
    use occupied
    implicit none
    integer, intent(in) :: A,pos
    integer, intent(inout):: nNeighb, nHbond
    integer  :: B,j,i
    real     :: inn
    
! Counts the interactions and the number of neighbours
    inn = 0
    ! nNeighb=0; nHbond=0
    if (A <= 4) then
        ! Horizontal case
        do i=1,4
            j = direction(i)
            B = occ(neighb_M (pos,j))
            if ( B > 0 .and. B <= 4 ) then 
                nNeighb = nNeighb + 1
                inn = inn + M_int(A,B,i)
                nHbond = nHbond + int(inn)
            end if
        end do
    else if ( A <= 8 ) then
        ! Vertical-right case
        do i=5,8
            j = direction(i)
            B = occ(neighb_M(pos,j)) - 4
            if ( B > 0 .and. B <= 4 ) then 
                nNeighb = nNeighb + 1
                inn = inn + M_int(A-4,B,i-4)
                nHbond = nHbond + int(inn)
            end if
        end do
    else ! if ( A <= 12 ) then
        ! Vertical-left case
        do i=9,12
            j = direction(i)
            B = occ(neighb_M(pos,j)) - 8
            if ( B > 0 .and. B <= 4 ) then 
                nNeighb = nNeighb + 1
                inn = inn + M_int(A-8,B,i-8)
                nHbond = nHbond + int(inn)
            end if
        end do
    end if
end subroutine

subroutine finish_statistics () 
    implicit none
    close(lust)
end subroutine finish_statistics  

end module statistics

