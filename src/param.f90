module param
  use strings
  implicit none
!... L      - grid dimension
  integer L,LL3,LL4
!
!... No_kmc - number of KMC steps
  integer No_kmc
!
!... Periodic Boundary Conditions (PBC)
  logical pbc
!
!... elementary barriers (in eV) and inverse temperature (in 1/eV), temperature (K)
!     1 - PTMDC diffusion via pivoting mechanism
!     2 - PTMDC-PTMDC detachement (a dimer); then additive
!     3 - PTMDC single isomerization
!     4 - PTMDC assisted isomerization
!     5 - PTMDC desorption from surface
!     6 - PTMDC extra isomerization barrier because in a linker
  real*8 elem_barriers(20),beta,Temp,pref,beta1,Temp1
  logical ignore_single_HB
!
!... 
  real*8,parameter :: Boltzm=  8.617343e-5  ! in eV/K
!
!... deposition
  real*8 depos_rate,coverage,max_coverage
  integer Max_no_of_M
  logical deposition
!... High Coverage or Low Coverage Phase
  logical allow_high_coverage
!... useful global-types
  real*8, parameter :: tiny=0.000001
  real*8   :: real_seed
  integer seed, iPrnt
!
!... for drawing
  logical images,rotate,show_site_num,multicolor,show_prop,create_jpgs
  integer No_draw,len_No_kmc,scale,kmc_to_draw,color
  integer tsize
!
!... for testing (interactive  regime)
  logical testing
!... what type of KMC: standard or fixed step
  logical standard_kmc
!... time step (in ps)
  real*8 dt
!... restart
  logical restart
  character restart_filename*50
!... kmc.out file formatting
  logical write_formatted, write_kmcout
!... strings for working with input
  character Line*200
  integer LinEnd(100),LinPos(100),NumLin,iErr
 
  CONTAINS

  subroutine input()
!..............................................................................
!....... read input info ......................................................
!..............................................................................
    integer i,l0
    character cha*10
    real*8 rat(10)
!
    open(1,file='input.dat', form='formatted')

    standard_kmc=.true.
    write(*,'(a)')'============================================='
    write(9,'(a)')'============================================='
    write(*,'(a)')'==> STYLE: standard KMC -> variable time-step'
    write(9,'(a)')'==> STYLE: standard KMC -> variable time-step'
    write(*,'(a)')'============================================='
    write(9,'(a)')'============================================='
!
!______________ restart from a file
    restart=.false.
    call find_string('restart',7,Line,1,.true.,iErr)
    if(iErr.ne.0) then
       write(*,'(a)')'... This is a calculation from scratch'
    else
       call CutStr(Line,NumLin,LinPos,LinEnd,0,0,iErr)
       if(NumLin.lt.2) go to 10
       restart_filename=Line(LinPos(2):LinEnd(2))
       write(*,'(2a)')'... restart from file ',trim(restart_filename)
       write(9,'(2a)')'... restart from file ',trim(restart_filename)
       restart=.true.
    end if
!
!______________ whether to write kmc.out, formatted or not
    write_kmcout = .false.
    call find_string('write_kmcout',12,Line,1,.true.,iErr)
    if (iErr == 0 ) then 
        write_kmcout = .true.
        call find_string('formatted',9,Line,1,.true.,iErr)
        if(iErr.ne.0) then
           write(*,'(a)')'... [kmc.out] will be UNFORMATTED'
           write(9,'(a)')'... [kmc.out] will be UNFORMATTED'
           write_formatted=.false.
        else
           write(*,'(a)')'... [kmc.out] will be FORMATTED'
           write(9,'(a)')'... [kmc.out] will be FORMATTED'
           write_formatted=.true.
        end if
    else 
       write(*,'(a)')'... [kmc.out] will not be written'
       write(9,'(a)')'... [kmc.out] will not be written'
    end if
!
!_______________ only PBC
!
    pbc=.true.
    write(*,'(a)') '... Periodic Boundary Conditions = ON'
    write(9,'(a)') '... Periodic Boundary Conditions = ON'
    write(9,'(a)') '... Molecules come from the opposite wall when crossing'
!
!_______________ supercell dimensions
!
    call find_string('dimension',9,Line,1,.true.,iErr)
    if(iErr.ne.0) go to 10
    call CutStr(Line,NumLin,LinPos,LinEnd,0,0,iErr)
    if(NumLin.lt.2) go to 10
    read(Line(LinPos(2):LinEnd(2)),*,err=10) L
    if(L.lt.2) then
       stop "wrong supercell size L"
    else
       write(*,'(2(a,i5))') '... Supercell size  = ',L,' x ',L
       write(9,'(2(a,i5))') '... Supercell size  = ',L,' x ',L
       LL3=L*L ;  LL4=(L+2)*(L+2) 
    end if
    write(9,'(a,i10)')'... Total number of sites available to PTMDC => ',LL3
    write(9,'(a,i10)')'... Total number of sites available for drawing => ',LL4
!
!_________________ deposition rate: in 1/ps and coverage
!
    deposition=.true.
    call find_string('deposition_rate',15,Line,1,.true.,iErr)
    if(iErr.ne.0) go to 10
    call CutStr(Line,NumLin,LinPos,LinEnd,0,0,iErr)
    if(NumLin.lt.2) go to 10
    read(Line(LinPos(2):LinEnd(2)),*,err=10) depos_rate
    if(depos_rate.eq.0.0) deposition=.false.
    write(*,'(a,e12.6,a)')'... deposition rate = ',depos_rate,' ps^{-1}'
    write(9,'(a,e12.6,a)')'... deposition rate = ',depos_rate,' ps^{-1}'
!
    max_coverage=1.0
    call find_string('coverage',8,Line,1,.true.,iErr)
    if(iErr.eq.0) then
       call CutStr(Line,NumLin,LinPos,LinEnd,0,0,iErr)
       if(NumLin.lt.2) go to 10
       read(Line(LinPos(2):LinEnd(2)),*,err=10) max_coverage
    end if
    Max_no_of_M=max_coverage*LL3
    write(*,'(a,i10)')'... Maximum number of M molecules to deposit = ',Max_no_of_M
    write(9,'(a,f10.5)')'... coverage = ',max_coverage
    write(9,'(a,i10)')'... Maximum number of M molecules to deposit = ',Max_no_of_M
!
!_________________ No of KMC steps
!
    call find_string('number_of_kmc_steps',19,Line,1,.true.,iErr)
    if(iErr.ne.0) go to 10
    call CutStr(Line,NumLin,LinPos,LinEnd,0,0,iErr)
    if(NumLin.lt.2) go to 10
    read(Line(LinPos(2):LinEnd(2)),*,err=10) No_kmc
    write(9,'(a,i10)') '... Number of KMC steps found = ',No_kmc
    write(cha,'(i10)') No_kmc
    call posit(cha,10,l0)
    len_No_kmc=10-l0+1
!
!_________________ if to do testing
!
    testing=.false.
    call find_string('do_testing',10,Line,1,.true.,iErr)
    if(iErr.eq.0) testing=.true.
    write(9,'(a,l1)') '... TESTING = ',testing
    ! TO CHECK: If read correctly the barriers
!
!_________________ elementary barrier (in eV)
!
!    elem_barriers=0.0
    call find_string('barriers',8,Line,1,.true.,iErr)
    if(iErr.ne.0) go to 10
    call CutStr(Line,NumLin,LinPos,LinEnd,0,0,iErr)
    if(NumLin.lt.5) go to 10
    ! read(Line(LinPos(2):LinEnd(5)),*,err=10) (elem_barriers(i),i=1,4)
    read(Line(LinPos(2):LinEnd(7)),*,err=10) (elem_barriers(i),i=1,6)
!
!_________________ temperature 
!
    call find_string('temperature',11,Line,1,.true.,iErr)
    if(iErr.ne.0) go to 10
    call CutStr(Line,NumLin,LinPos,LinEnd,0,0,iErr)
    if(NumLin.lt.2) go to 10
    read(Line(LinPos(2):LinEnd(2)),*,err=10) Temp
    beta=1.0/(Boltzm*Temp)
!
!________ prefix to the rate, i.e. the exp prefactor to the rate:
!         - in inverse ps (1 ps=10^{-12} s)
    pref= +1.0 * log(10.0)
!
!...................... RATES TABLE ...............................................................
!
    do i=1,5
       rat(i)=exp(-beta*elem_barriers(i)+pref)
    end do
    rat(8)=depos_rate   ;    elem_barriers(8)=0.0d0
    write(9,'(/a)')'============|  Barriers (in Ev)  and Rates (in ps^-1) |============'
    write(9,'(a,f10.3,a)') '   Temperature = ',temp,' K'
    write(9,'(a,f10.5,x,e12.6)') ' 1. PTMDC diffusion via pivoting:  ',elem_barriers(1),rat(1)
    write(9,'(a,f10.5,x,e12.6)') ' 2. PTMDC-PTMDC single detachment: ',elem_barriers(2),rat(2)
    write(9,'(a,f10.5,x,e12.6)') ' 3. PTMDC single isomerization:     ',elem_barriers(3),rat(3)
    write(9,'(a,f10.5,x,e12.6)') ' 4. Extra barrier due assisted isomer.:',elem_barriers(4),rat(4)
    write(9,'(a,f10.5,x,e12.6)') ' 5. PTMDC desorption from surface: ',elem_barriers(5),rat(5)
    write(9,'(a,f10.5,x,e12.6)') ' 6. Extra end-isomer. barrier in a linker: ',elem_barriers(6),rat(6)
    write(9,'(a,f10.5,x,e12.6)') ' 8. PTMDC deposition on surface:   ',elem_barriers(8),rat(8)
    write(9,'(a/)')'==================================================================='
!
!_________________ Ignoring Single HB 
!
    ignore_single_HB = .false.
    call find_string('ignore_single_HB',16,Line,1,.true.,iErr)
    if(iErr.eq.0) ignore_single_HB = .true.
    write(9,'(a,l1)')'... Ignore Single H-bond = ',ignore_single_HB
!
!________________ level of print: 0 - very little; 5 - a lot
! 
   iPrnt=0
    call find_string('printing',8,Line,1,.true.,iErr)
    if(iErr.eq.0) then
       call CutStr(Line,NumLin,LinPos,LinEnd,0,0,iErr)
       if(NumLin.lt.2) go to 10
       read(Line(LinPos(2):LinEnd(2)),*,err=10) iPrnt
    end if
    write(9,'(a,i1)')'... Printing level = ',iPrnt

!_________________ if to show site numbers on frames
    show_site_num=.false.
    call find_string('show_site_numbers',17,Line,1,.true.,iErr)
    if(iErr.eq.0) show_site_num=.true.
    if(show_site_num) then
       write(9,'(a,l1)') '... Site numbers WILL be shown'
    else
       write(9,'(a,l1)') '... Site numbers will NOT be shown'
    end if

!__________________ color of PTMDC
    ! color=4
    multicolor=.false.
    call find_string('multicolor',10,Line,1,.true.,iErr)
    if(iErr.eq.0) multicolor = .true. 
    if (multicolor) then
         write(9,'(a,i10,a)')'... You choose Multicolor: (Inspired on the Day Of The Tentacle)' 
         write(9,'(a,i10,a)')'....... Cis : Blue rectangles, pink circles ' 
         write(9,'(a,i10,a)')'....... L-trans: Magenta rectangles, green circles' 
         write(9,'(a,i10,a)')'....... D-trans: Green rectangles, yellow circles' 
    else 
         write(9,'(a,i10,a)')'...  Multicolor not used' 
    end if

!_________________ after how many KMC steps to print frames
    No_draw=100
    call find_string('drawing_frequency',17,Line,1,.true.,iErr)
    if(iErr.eq.0) then
       call CutStr(Line,NumLin,LinPos,LinEnd,0,0,iErr)
       if(NumLin.lt.2) go to 10
       read(Line(LinPos(2):LinEnd(2)),*,err=10) No_draw
    end if
    write(9,'(a,i10,a)')'... Printing every ',No_draw,' frames'
          
!_________________ drawing scale
    scale=800
    call find_string('drawing_scale',13,Line,1,.true.,iErr)
    if(iErr.eq.0) then
       call CutStr(Line,NumLin,LinPos,LinEnd,0,0,iErr)
       if(NumLin.lt.2) go to 10
       read(Line(LinPos(2):LinEnd(2)),*,err=10) scale
    end if
    write(9,'(a,i5)')'... Scale for drawing frames = ',scale

!_________________ size for drawing the title
    tsize=20
    call find_string('drawing_title_size',18,Line,1,.true.,iErr)
    if(iErr.eq.0) then
       call CutStr(Line,NumLin,LinPos,LinEnd,0,0,iErr)
       if(NumLin.lt.2) go to 10
       read(Line(LinPos(2):LinEnd(2)),*,err=10) tsize
    end if
    write(9,'(a,i5)')'... Size for drawing titles on frames = ',tsize

!_________________ if to show images
    images=.false.
    call find_string('show_images',11,Line,1,.true.,iErr)
    if(iErr.eq.0) images=.true.
    write(9,'(a,l1)') '... Drawing: images=',images
!_________________ if to create JPGs images
    create_jpgs=.true.
    call find_string('no_jpgs',7,Line,1,.true.,iErr)
    if(iErr.eq.0) create_jpgs=.false.
    write(9,'(a,l1)') '... Drawing: Create JPGs=', create_jpgs
!_________________ if to show proportions
    show_prop =.false.
    call find_string('show_proportions',16,Line,1,.true.,iErr)
    if(iErr.eq.0) show_prop=.true.
    write(9,'(a,l1)') '... Show isomers proportion= ',show_prop
!
!_________________ if allow for High Coverage Phase
    allow_high_coverage =.false.
    call find_string('allow_high_coverage',19,Line,1,.true.,iErr)
    if(iErr.eq.0) allow_high_coverage=.true.
    write(9,'(a,l1)') '... Allow High Coverage Phase = ',allow_high_coverage 
!
!___________ Choose fixed or random seed for random number generation
    seed= 0
    call find_string('seed',4,Line,1,.true.,iErr)
    if(iErr.eq.0) then
       call CutStr(Line,NumLin,LinPos,LinEnd,0,0,iErr)
       if(NumLin.lt.2) go to 10
       read(Line(LinPos(2):LinEnd(2)),*,err=10) seed
    end if
    if (seed > 0 ) then 
        seed = -seed
        write(9,'(a,i6)')'... Seed used for Random Number Generation = ', -seed
    else 
!
!__________ Change the seed for random numbers generation
!
!__________ initialise seed
        call init_random_seed()
        call random_number(real_seed)
        seed = -int(real_seed*10000)-1
        write(9,'(a,i6)')'... Random seed used for Random Number Generation = ', -seed
    end if 
!
    close (1)
    write(*,'(/a/)')'... Input file has been read in successfully ...'
    write(9,'(/a/)')'... Input file has been read in successfully ...'

    return
!................ errors in reading input
10  write(*,*)'FATAL errors in reading the input file!'
    stop 'FULL STOP!'
  end subroutine input

end module param
