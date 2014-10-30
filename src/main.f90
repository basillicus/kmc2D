program KMC_main
  use occupied ! occupied uses: param, moves and interactions
  use change
  use statistics
  ! use interactions
  implicit none
!
!... L      - surface dimension
!... working
  integer kmc,m,j,kmc0,l0
  character cha*10
  real*8 time

  call flush ()
!
!........ open the main output file
  open(9,file='OUTPUT.kmc')
!
  call hat()
!
!....... read input info
!
  call input()
!
!....... setup the lattice
!
  call setup()
!
!.........................................................................
!............ (i) construct the lattice and all info arrays ..............
!............ (ii) find all neighbouring sites for PTCDA's
!.........................................................................
!
  allocate(occ(0:LL3))  
  call do_info()
!
!
!............. initialice the Interaction's matrix
!
  call init_interactions(ignore_single_HB)
!
!............. initialice the statistics module
!
  if (do_statistics) call init_statistics ()

  if(restart) then
     open(1,file=trim(restart_filename),err=100)
     read(1,*) kmc0,time ; read(1,*) occ  
     close(1)
     write(*,'(a)')'WARNING: file ['//trim(restart_filename)//'] will be overwritten'
     write(*,'(a)')'/bin/cp '//trim(restart_filename)//' '//trim(restart_filename)//'.bak'
     call system('/bin/cp '//trim(restart_filename)//' '//trim(restart_filename)//'.bak' )
     write(cha,'(i10)') No_kmc+kmc0
     call posit(cha,10,l0)
     len_No_kmc=10-l0+1
  else
     occ=0  ; kmc0=0 ; time=0.0
  end if
! 
!.........................................................................
!...................KMC starts here ......................................
!.........................................................................
!
  allocate(empty_info(LL3)) 
!...... maximum # of possible PTMDC moves:
!  - 4*LL3 diffusion moves
!  - 2*LL3 isomerizations
!  - 1 adsorption move (at random)
!  - LL3 desorption moves
  max_moves=LL3*7 + 1
  if(.not.standard_kmc) max_moves=max_moves+1
  allocate(move(max_moves)) 
!
!................ if to do TEST
!
  if(testing) then
     call do_test(kmc0) ; stop
  end if
!
!................................................................................
!................ if to do real KMC .............................................
!................................................................................
!
  if(.not.write_formatted) then
     if(restart) then
        open(10,file='kmc.out',form='unformatted',err=100,position='append')
     else
        open(10,file='kmc.out',form='unformatted',err=100)
     end if
  else
     if(restart) then
        open(10,file='kmc.out',form='formatted',err=100,position='append')
     else
        open(10,file='kmc.out',form='formatted',err=100)
     end if
  end if
!
  write(*,'(a,i10,a,f10.3,a,f18.3)')'.......INITIAL KMC = ',kmc,', T = ',Temp,' K .....Time Interval =.', time_interval
  DO kmc=1+kmc0,No_kmc+kmc0
     if(mod(kmc,No_draw).eq.0) &
         write(*,'(a,i10,a,f10.3,a)')'....... KMC = ',kmc,', T = ',Temp,' K ......'
         ! write(*,'(a,f18.3,a,f18.3,a,f18.3)')'Time ', time,',Temperature ', temp, 'K, Next Time: ', next_time
     if(iPrnt.ge.2) write(9,*)'............. KMC = ',kmc,' ...............'
!
!_________________ create/update moves
     call check_coverage(kmc,occ,time)
     call all_moves()      

!_________________ chosee a move
     call which_move(m)         

!_________________ perform the move
     call perform_move(m,occ,kmc) 
!_________________ do statistics
     if (do_statistics) call perform_statistics (m,kmc,occ,time,temp)

!_________________ draw
     if(mod(kmc,No_draw).eq.0) then
        kmc_to_draw=kmc
        call draw(kmc,occ,.false.,time,Temp)
     end if
    !  time=time+dt

    call updateTime (time)
    call updateTemperature (kmc)
!
!________________________ simple error check
     call error_check(occ,kmc)
!
!________________________ save to the file the current configuration
     if (write_kmcout) then
         if(.not.write_formatted) then
            write(10) kmc,dt,time,no_chg,(i_chg(j),occ_chg(j),j=1,no_chg), &
                 m,move(m)%type,move(m)%ini,move(m)%fin,move(m)%dir
         else
            write(10,'(i20,2(x,e12.6))') kmc,dt,time
            write(10,'(i10,x,a10,3(x,i10))') m,move(m)%type,move(m)%ini,move(m)%fin, &
                 move(m)%dir
            write(10,'(7(i10,x))') no_chg,(i_chg(j),occ_chg(j),j=1,no_chg)
         end if
     end if
!
     if(iPrnt.ge.2) then
        write(9,'(i20,2(x,e12.6))') kmc,dt,time
        write(9,'(i10,x,a10,3(x,i10),x,e12.6)') m,move(m)%type,move(m)%ini,move(m)%fin, &
             move(m)%dir,move(m)%rate
     end if
  END DO
!
  kmc_to_draw=kmc
  call draw(kmc,occ,.false.,time,temp)
  close (10) 
  write(9,'(a)')'_____}> going to [finish] to STOP <{________ normal END'
  call finish(occ,No_kmc+kmc0,time)
!
!............ error
!
100 write(*,*) 'FATAL! Cannot open restart file '//trim(restart_filename)
  stop 'stop in main.f90 by error'
end program KMC_main

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!< T E S T >!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine do_test(kmc0)
    use occupied
    implicit none
    integer,dimension(:),allocatable:: occ_old
    character cha,answer,cha5*5,filename*15
    integer kmc,i,m,oc,m1,kmc_jump,kmc0
    real*8 time
!
!..................... drop the molecules first
!
    occ=0 ; time=0.0 ; i = 1
    write(*,*) ' INSRUCTIONS:'
    write(*,*) ' -1 : End droping molecules'
    write(*,*) ' -2 : Save current configuration to a file'
    write(*,*) ' -3 : Load old configuration from a file'

    do while (i >= 0 )
    write(*,*)'Drop: give position i and occ(i): (Type -1 for ending)'
        read(*,*,err=10) i
        if (i > 0 ) then
             read(*,*,err=10) oc
             if(i.ge.1 .and. i.le.LL3 .and. oc.ge.0 .and. oc.le.12) then
                 occ(i)=oc ;     kmc_to_draw=1
                 call draw(1,occ,testing,time,temp)
             else 
                 write(*,*) "Wrong i or occ(i). Molecule not added"
             end if
         else if (i == -2 ) then
             write(*,*) 'Give a filename to save the current configuration:'
             read(*,'(15a)') filename 
             open(33, file=trim(filename),form='formatted',status='new')
             write(33,*)  " 0    0.000" ! kmc and time
             write(33,'(40i3)') occ
             close(33)
             i = 0  ! Continue dopositing molecules
         else if (i == -3 ) then
             write(*,*) 'Give a filename to read from the configuration:'
             read(*,'(15a)') filename 
             open(33, file=trim(filename),form='formatted',status='old')
             read(33,*) kmc, time
             read(33,'(40i3)') occ
             close(33)
             call draw(1,occ,testing,time,temp)
             i = 0  ! Continue dopositing molecules
        end if
    end do
!
!.................... make moves (including drops of molecules)
!
10  allocate(occ_old(0:LL3))
    occ_old=occ
    write(*,'(/a)') ' Jump to specific KMC step number (0 for NOT):'
    read(*,*) kmc_jump
    do kmc=1+kmc0,No_kmc+kmc0
       if(kmc.ge.kmc_jump) then
          if(iPrnt.ge.2) then
             write(*,*)'.................. occupancies ..................'
             write(9,*)'.................. occupancies ..................'
             do i=1,LL3
                if(occ(i).ne.0 .or. occ_old(i).ne.0) then
                   cha=' ' ; if(occ(i).ne.occ_old(i)) cha='*' 
                   write(*,'(a,i5,a,i2,3x,a)') 'occ(',i,') = ',occ(i),cha
                   write(9,'(a,i5,a,i2,3x,a)') 'occ(',i,') = ',occ(i),cha
                end if
             end do
          end if
       end if
       write(*,*)'............. KMC = ',kmc,' ...............'
       call all_moves()
       kmc_to_draw=kmc
       call draw(kmc,occ,.true.,time,temp)
       call which_move(m1)
       if(kmc.ge.kmc_jump) then
1         write(*,'(/a)') 'Choose which move to accept: '
          write(*,'(a)') ' - Enter move number from the list'
          write(*,'(a)') ' OR'
          write(*,'(a)') ' - Enter any letter to choose move number at random via KMC algorithm'
          write(*,'(a)') '=================>'
          read(*,'(a)') cha5
          read(cha5,'(i5)',err=3) m1
          if(m1.lt.0 .or. m1.gt.No_moves) go to 1
          go to 3
       end if
3      m=m1
       write(*,'(a,i5)')'===> The move chosen = ',m
       occ_old=occ
       call perform_move(m,occ,kmc)
       call updateTime (time)
       kmc_to_draw=kmc
       call draw(kmc,occ,testing,time,temp)
    end do
  end subroutine do_test

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine updateTime (time)
      use param
      
      implicit none
      real*8, intent(inout) :: time

      time = time + dt 

  end subroutine

  subroutine updateTemperature (kmc)
      use param

      implicit none
      ! real*8 :: time
      integer :: kmc
      
      
      ! Update by time step
      ! if (time >= next_time) then
      !     next_time = time + time_interval
      !     Temp = Temp + temperature_ramp
      !     beta = 1.0/(Boltzm*Temp)
      ! end if

      ! Update by kmc steps
      if (kmc >= next_step) then
          next_step = kmc + steps_interval 
          write(*,*) 'Next T update at ', next_step
          Temp = Temp + temperature_ramp
          beta = 1.0/(Boltzm*Temp)
      end if

  end subroutine 

  subroutine finish(occ,kmc,time)
    use param
    implicit none
    integer occ(0:LL3),kmc
    real*8 time
    open(3,file='restart.kmc')
    write(3,*) kmc,time
    write(3,'(40i3)') occ
    close (3)
    write(9,'(a,i20)')'______]> STOP in [finish] at KMC=',kmc
    close (9) ; stop 'STOP in finish'
  end subroutine finish
  
  subroutine check_coverage(kmc,occ,time)
    use info
    implicit none
    integer kmc,no_m,i,occ(0:LL3)
    real*8 time
    if(deposition) then
       no_m=0
       do i=1,LL3
          if(occ(i).gt.0) no_m=no_m+1
       end do
       if (no_m.eq.Max_no_of_M .or. no_m.ge.LL3) then
          write(*,*)'WARNING: The coverage target reached at KMC step ',kmc
          write(*,*)'         NO more deposition!'
          write(9,*)'WARNING: The coverage target reached at KMC step ',kmc
          write(9,*)'         NO more deposition!'
          deposition=.false.
!____________ save the current configuration
          call save_configuration(occ,kmc,time)
       end if
    end if
  end subroutine check_coverage
  
  subroutine save_configuration(occ,kmc,time)
    use param
    implicit none
    integer occ(0:LL3),kmc
    real*8 time
    open(3,file='cover_reached.kmc')
    write(3,*) kmc,time
    write(3,'(40i3)') occ
    close (3)
  end subroutine save_configuration

  subroutine error_check(occ,kmc)
    use param
    integer occ(0:LL3),i,kmc
    do i=1,LL3
       if(occ(i).gt.12 .or. occ(i).lt.0) then
          write(9,*) 'error [KMC main] kmc=',kmc,i,occ(i)
          stop 'STOP in main: error in occ found'
       end if
    end do
  end subroutine error_check

  subroutine do_info()
    use info
    implicit none
!
!.........................................................................
!............ (i) construct the lattice and all info arrays ..............
!............ (ii) find all neighbouring sites for PTCDA's
!.........................................................................
!
    allocate(x_lat(LL3,2))  ; allocate(x_lat_ext(LL4,2)) ;  allocate(n1(LL3))
    allocate(n2(LL3)) ; allocate(neighb_M(1:LL3,1:8)) 

!....... attempt
    write(*,'(a)')'... Attempting to read the info data from the file [info.dat] ...'
    open(11,file='info.dat',form='unformatted',status='old',err=1)
    read(11,err=1,end=1) AL,BL,x_lat_ext,x_lat,n1,n2,neighb_M
    close(11)
    write(*,'(a)')'... Done! Proceed ..'
    return

!....... if unsuccessful, calculate the info and replace the file
1   write(*,*)'... File [info.dat] is absent or bad; calculating the info NOW ...'
    close(11,status='delete')
    open(11,file='info.dat',form='unformatted',status='new')
    call build_latt_pbc()
    close(11)
    write(*,'(a)')'... Done! Proceed ..'
  end subroutine do_info
