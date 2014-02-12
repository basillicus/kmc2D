subroutine hat()
  implicit none
  character line*30

      write(*,*)'.....................................................'
      write(*,*)'.....................................................'
      write(*,*)'.......         Kinetic Monte Carlo            ......'
      write(*,*)'.......   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^   ......'
      write(*,*)'.......         Hexagonal lattice              ......'
      write(*,*)'.......   Designed for PTMDC@Ag(111) problem   ......'
      write(*,*)'.......   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^   ......'
      write(*,*)'.....................................................'
      write(*,*)'>>>>>>>>  1st release: 15.11.2010'
      write(9,*)'.....................................................'
      write(9,*)'.....................................................'
      write(9,*)'.......         Kinetic Monte Carlo            ......'
      write(9,*)'.......   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^   ......'
      write(9,*)'.......         Hexagonal lattice              ......'
      write(9,*)'.......     Designed for PTMDC@Ag(111) problem ......'
      write(9,*)'.......   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^   ......'
      write(9,*)'.....................................................'
      write(9,*)'>>>>>>>>  1st release: 15.11.2010'
!
      include 'build.hat'
      call system('date > tmp.date')
      open(1,file='tmp.date')
      read(1,'(a)') line
      close (1,status='delete')
      write(*,*)'>>>>>>>>  today      : '//line
      write(*,*)'.....................................................'
      write(*,*)'......      inquiries to L. Kantorovich     .........'
      write(*,*)'......      lev.kantorovitch@kcl.ac.uk      .........'
      write(*,*)'.....................................................'
!
      write(9,*)'>>>>>>>>  today      : '//line
      write(9,*)'.....................................................'
      write(9,*)'......      inquiries to L. Kantorovich     .........'
      write(9,*)'......      lev.kantorovitch@kcl.ac.uk      .........'
      write(9,*)'.....................................................'

end subroutine hat



