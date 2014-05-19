  subroutine draw(kmc,occ,test,time)
    use info
!.............................................................................
! Draw the current frame of the system into frame#.fig file
! and previe, where # is the kmc number
!.............................................................................
    implicit none
    integer kmc,l0,k,i,i1,i2,occ(0:LL3),len,diff
    ! Number of the molecules adsorbed
    integer :: nCis=0, nLtrans=0, nDtrans=0, nTotal=0
    character cha*10
    real*8 r(2),time
    real*8 :: x1_h,y1_h,x2_h,y2_h, x1_v,y1_v,x2_v,y2_v
    logical if_in_latt,test
!
    if(test) then
       open(1,file='frame.fig')
    else
       write(cha,'(i10)') kmc
       call posit(cha,10,l0)  ; len=10-l0+1 ; diff=len_No_kmc-len
       do k=1,diff
          cha(10-len_No_kmc+k:10-len_No_kmc+k)='0'
       end do
       l0=l0-diff 
       open(1,file='frame'//cha(l0:)//'.fig')
    end if
    write(1,'(a)') '#FIG 3.2  Produced by xfig version 3.2.5'
    write(1,'(7(a/),a)') 'Portrate','Center','Inches','Letter','100.00','Single','-2','1200 2'
!
!............. write the KMC step and the time in ns
!    
    call title(time,AJ)
!
!............. plot the lattice first
!    
    do k=-1,L+1
       ! Vertical lines
        x1_v = k*AJ(1,1) - AJ(2,1)
        y1_v = AJ(2,2)
        x2_v = k*(AJ(1,1)) + (L+1)*AJ(2,1)
        y2_v = -(L+1)*AJ(2,2)
       call draw_line(x1_v,y1_v,x2_v,y2_v)

       ! Horizontal Lines
       x1_h = k*AJ(2,1)-1 
       y1_h = -k*AJ(2,2)
       x2_h = k*AJ(2,1) + (L+1)*AJ(1,1)
       y2_h = -k*AJ(2,2) 
       call draw_line(x1_h,y1_h,x2_h,y2_h)

       ! Diagonal Lines
       call draw_line(x1_h,y1_h,x1_v,y1_v )
       call draw_line(x2_h,y2_h,x2_v,y2_v )
    end do
!
!............. indicate the cell
!
    call cell(AL)
!
!............. plot the molecules and counts the number of isomers
!
    nTotal  = 0
    nCis    = 0
    nLtrans = 0
    nDtrans = 0

    do i=1,LL3
       if(show_site_num) call number(x_lat(i,1),-x_lat(i,2),i)
       if(occ(i).gt.0) then
          call ptpmd(x_lat(i,1),-x_lat(i,2),occ(i))
          ! Count the number of each isomer
          if (show_prop) then 
              if (occ(i) == 3 .or. occ(i) == 7 .or. occ(i) == 11) then
                    nLtrans = nLtrans + 1        
              else if (occ(i) == 4 .or. occ(i) == 8 .or. occ(i) == 12) then
                    nDtrans = nDtrans + 1        
              else 
                    nCis = nCis + 1
              end if 
          nTotal = nTotal + 1
          end if
          if(test) call number(x_lat(i,1),-x_lat(i,2),i)
       end if
    end do

    if (show_prop) call draw_proportions (nCis,nLtrans,nDtrans,nTotal, AJ )

!
!............ plot images of molecules and atoms within the extended lattice
!
    if(images) then
       do i1=-1,L+1
          do i2=-1,L+1
             r(:)=i1*AJ(1,:)+i2*AJ(2,:)
             if( .not. if_in_latt(i1,i2) ) then
                call find_site(r,i)
                if(occ(i).gt.0) then
                   ! write a molecule of para-terphenyl-m-Dicarbonitrile(PTPMD)
                   call ptpmd(r(1),-r(2),occ(i))
                   if(test) call number(r(1),r(2),i)
                end if
             end if
          end do
       end do
    end if
10  close(1)
    if(.not.test)  write(*,*)'File frame'//cha(l0:)//'.fig written ...'
!
!............. quick preview
!
    if(test) then
       call system('fig2dev -L jpeg frame.fig frame.jpg')
!       call system('fig2dev -L jpeg -P frame.fig frame.png')
!       call system('xv frame'//cha(l0:)//'.jpg')
    else if (create_jpgs) then
       call system('fig2dev -L jpeg -q 50 frame'//cha(l0:)//'.fig frame'//cha(l0:)//'.jpg')
!       call system('fig2dev -L jpeg frame'//cha(l0:)//'.fig frame'//cha(l0:)//'.png')
    end if
!    if(.not.test) call system('rm -rf frame*fig')
  end subroutine draw

  subroutine ptcda(x,y,oc)
    implicit none
!.............................................................................
! Plot PTCDA of occupany oc with the centre at (x,y) 
!.............................................................................
    real*8 x,y
    integer oc
!_________ draw the rectangular
    call rectangular(x,y)
    if(oc.eq.2) then
       call circle(x,y-0.4)
    else if(oc.eq.3) then   
       call circle(x,y+0.4)
    else if(oc.eq.4) then   
       call circle(x,y-0.4)
       call circle(x,y+0.4)
    end if
  end subroutine ptcda
  
  subroutine ptpmd(x,y,config)
    implicit none
!.............................................................................
! Plot PTPMD of occupany oc with the first point at (x,y) 
!.............................................................................
    real*8 x,y
    integer oc, config
    call molecule (x,y,config)
  end subroutine ptpmd

  subroutine molecule (x,y,config)
!.............................................................................
!  draw a molecule placed  at (x,y)
!.............................................................................
    use param
    use lattice
    implicit none
    real*8 x,y, long
    integer config
    ! coordinates for the rectangle
    integer i1(2),i2(2),i3(2),i4(2) 
    real :: i5(2),i6(2) 
    ! coordinates for the circles
    integer ic1,jc1,ic2,jc2,k
    real*8,parameter :: rad=0.10
    ! Colors of Rectangles and Circles
    integer, parameter :: blue=1, yellow=6,green=14, magenta=21, pink=27
    integer ::   col_Rect=5, col_Circ=14

    ! Radius of the circles
    k=rad*scale 

    ! Horizontal configuration
    if ( config <= 4 ) then 
        !
        !----------+-------------+-------------+---             
        !         ·             ·             ·                                   
        !        ·             ·             ·                                    
        !      1            2 ·             ·                  
        !      ///////////// ·             ·                   
        !-site>X////////////+-------------+---             
        !      /////////////            ·                      
        !      /////////////           ·                       
        !    · 4            3         ·                        
        !   ·             ·          ·                                      
        !  ·             ·          ·                                       
        !------+------------+------+----------             
        !                                                                          
        ! Rectangle position         
        ! i1(1)=(x)*scale   ; i1(2)=(y-0.1)*scale 
        ! i2(1)=(x+1)*scale ; i2(2)=(y-0.1)*scale 
        ! i3(1)=(x+1)*scale ; i3(2)=(y+0.1)*scale 
        ! i4(1)=(x)*scale   ; i4(2)=(y+0.1)*scale 
        ! Draw smaller rectangles
        i1(1)=(x)*scale   ; i1(2)=(y-0.1)*scale 
        i2(1)=(x+0.7)*scale ; i2(2)=(y-0.1)*scale 
        i3(1)=(x+0.7)*scale ; i3(2)=(y+0.1)*scale 
        i4(1)=(x)*scale   ; i4(2)=(y+0.1)*scale 
        ! Circles Positions as diferent isomer
        if ( config == 1 ) then; 
            ic1=i4(1)  ; jc1=i4(2)  
            ic2=i3(1)  ; jc2=i3(2)  
        else if ( config == 2 ) then 
            ic1=i1(1)  ; jc1=i1(2)  
            ic2=i2(1)  ; jc2=i2(2)  
        else if ( config == 3 ) then 
            ic1=i4(1)  ; jc1=i4(2)  
            ic2=i2(1)  ; jc2=i2(2)  
        else if ( config == 4 ) then 
            ic1=i1(1)  ; jc1=i1(2)  
            ic2=i3(1)  ; jc2=i3(2)  
        end if
    !
    ! Vertical-right configuration
        !--------4---(5)--3-------+-------------+---             
        !        /////////       ·             ·                          
        !       /////////       ·             ·                         
        !      /////////       ·             ·                        
        !     /////////       ·             ·                       
        !    /////////       ·             ·                      
        !   /////////       ·             ·                     
        !-site(X)--2-------+-------------+---             
        !  1              ·             ·                           
        !                ·             ·                            
        !               ·             ·                             
        !              ·             ·                              
        !             ·             ·                               
        !            ·             ·                                
        !-----------+-------------+---             
        !                                                                      
     else if ( config <= 8 ) then
        ! Rectangle position
        ! i5: The center position of the other extrem
        ! i5(1)=x+AJ(2,1); i5(2)=y-AJ(2,2)
        i5(1)=x+0.7*AJ(2,1); i5(2)=y-0.7*AJ(2,2)
        i1(1)=(x-0.1)*scale ; i1(2)=(y)*scale 
        i2(1)=(x+0.1)*scale ; i2(2)=(y)*scale 
        i3(1)=(i5(1)+0.1)*scale ; i3(2)=i5(2)*scale 
        i4(1)=(i5(1)-0.1)*scale ; i4(2)=i5(2)*scale 

        if ( config == 5 ) then; 
            ic1=i3(1)  ; jc1=i3(2)  
            ic2=i2(1)  ; jc2=i2(2)  
        else if ( config == 6 ) then 
            ic1=i1(1)  ; jc1=i1(2)  
            ic2=i4(1)  ; jc2=i4(2)  
        else if ( config == 7 ) then 
            ic1=i4(1)  ; jc1=i4(2)  
            ic2=i2(1)  ; jc2=i2(2)  
        else if ( config == 8 ) then 
            ic1=i1(1)  ; jc1=i1(2)  
            ic2=i3(1)  ; jc2=i3(2)  
        end if
     else if ( config >= 9 ) then
    ! Vertical-left configuration
        !-4---(5)--3----------------+-------------+---             
        !  \\\\\\\\\             ·             ·                          
        !   \\\\\\\\\           ·             ·                         
        !    \\\\\\\\\         ·             ·                        
        !     \\\\\\\\\       ·             ·                       
        !      \\\\\\\\\     ·             ·                      
        !       \\\\\\\\\   ·             ·                     
        !--------1--(X)--2-+-------------+---             
        !                 ·             ·                           
        !                ·             ·                            
        !               ·             ·                             
        !              ·             ·                              
        !             ·             ·                               
        !            ·             ·                                
        !-----------+-------------+---             
        !                                                                      
        ! Rectangle position
        ! i5: The center position of the other extrem
        i5(1)=x-0.7*AJ(2,1); i5(2)=y-0.7*AJ(2,2)
        i1(1)=(x-0.1)*scale ; i1(2)=(y)*scale 
        i2(1)=(x+0.1)*scale ; i2(2)=(y)*scale 
        i3(1)=(i5(1)+0.1)*scale ; i3(2)=i5(2)*scale 
        i4(1)=(i5(1)-0.1)*scale ; i4(2)=i5(2)*scale 

        if ( config == 9 ) then; 
            ic1=i3(1)  ; jc1=i3(2)  
            ic2=i2(1)  ; jc2=i2(2)  
        else if ( config == 10 ) then 
            ic1=i1(1)  ; jc1=i1(2)  
            ic2=i4(1)  ; jc2=i4(2)  
        else if ( config == 11 ) then 
            ic1=i4(1)  ; jc1=i4(2)  
            ic2=i2(1)  ; jc2=i2(2)  
        else if ( config == 12 ) then 
            ic1=i1(1)  ; jc1=i1(2)  
            ic2=i3(1)  ; jc2=i3(2)  
        end if
      end if

    if (multicolor) then
        ! L-trans colors (Tentacle 1)
        if (config == 3 .or. config == 7 .or. config == 11) then
            col_Rect = magenta
            col_Circ = green
        ! D-trans colors (Tentacle 2)
        else if (config == 4 .or. config == 8 .or. config == 12) then
            col_Rect =  green
            col_Circ =  yellow
        ! Cis colors  (Benjamin Franklin)
        else 
            col_Rect = blue
            col_Circ = pink
        end if 
    end if 
    !  ! IT WORKS----------------------------
    !  ! Write the rectangle
    !  write(1,'(a)') '2 1 0 2 0 5 50 -1 20 0.000 0 0 -1 0 0 5'
    !  write(1,'(9x,10i10)') i1(1),i1(2),i2(1),i2(2),i3(1),i3(2),i4(1),i4(2),i1(1),i1(2)
    !  ! Write the Circles
    !  write(1,'(a,8(i10))')'1 3 0 2 0 14 50 -1 20 0.000 1 0.0000 ',ic1,jc1,k,k,ic1,jc1,ic1+k,jc1
    !  write(1,'(a,8(i10))')'1 3 0 2 0 14 50 -1 20 0.000 1 0.0000 ',ic2,jc2,k,k,ic2,jc2,ic2+k,jc2
    !  ! IT WORKS----------------------------

    write(1,100) col_Rect
    100 format ( '2 1 0 2 0', I3 , ' 50 -1 20 0.000 0 0 -1 0 0 5')
    write(1,'(9x,10i10)') i1(1),i1(2),i2(1),i2(2),i3(1),i3(2),i4(1),i4(2),i1(1),i1(2)

    ! Write the Circles
    write(1,101) col_Circ,ic1,jc1,k,k,ic1,jc1,ic1+k,jc1
    101 format ('1 3 0 2 0',I3,' 50 -1 20 0.000 1 0.0000 ' , 8I10 )
    write(1,102) col_Circ,ic2,jc2,k,k,ic2,jc2,ic2+k,jc2
    102 format ('1 3 0 2 0',I3,' 50 -1 20 0.000 1 0.0000 ' , 8I10 )
  
  end subroutine molecule 
