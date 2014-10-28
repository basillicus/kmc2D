module statistics
    use param
    use moves
    
    implicit none
    integer :: nCis,nLtrans,nDtrans,nTotal
    integer :: nHbond,nNeighb
    integer, save :: mDiff=0,mIsomer=0,mAdsorb=0
    integer, parameter :: lust=34    ! Logic Unit for statistcs

CONTAINS

subroutine init_statistics () 
    implicit none

    if (restart) then 
        open(lust,file='statistics.kmc',position='append')
    else
        open(lust,file='statistics.kmc',status='new')
        write(lust,'(a)') '# Statistics of KMC'
        write(lust,'(a)') '# Time    T     Cis    L-trans    D-Trans    Total M   Diffusions    Isomerizations    N-Hbond    N-Neighbors'

    end if

end subroutine init_statistics  

subroutine perform_statistics (m,kmc,occ,time,temp)
    implicit none
    integer :: m,kmc
    integer :: occ(0:LL3) 
    real*8  :: time, temp

end subroutine perform_statistics 

subroutine finish_statistics () 
    implicit none
    close(lust)
end subroutine finish_statistics  

end module statistics

