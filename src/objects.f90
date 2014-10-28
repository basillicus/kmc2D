  subroutine title(time,AJ)
!.............................................................................
!  write the number ii
!.............................................................................
    use param
    implicit none
    real*8 time,AJ(2,2),tm
    integer i,j,l0,l1,j1
    character cha*10,cha1*5,cha2*16,ch*4

!... to print on the right side
    ! i= ((L+1)*AJ(1,1)+2)*scale ; j1=1*AJ(2,2)*scale
    i= ((L+1)*AJ(1,1)+5)*scale ; j1=-5*AJ(2,2)*scale
    write(cha,'(i10)') kmc_to_draw ;     call posit(cha,10,l0)
    write(cha1,'(i5)') tsize    ;     call posit(cha1,5,l1)
    write(1,'(a,2(i10,x),a)') '4 0 0 44 -1 26 '//cha1(l1:)//' 0.0000 4 225 390 ',&
             i,j1,'Step '//cha(l0:)//char(92)//'001'

!..... work out the time
    if(time.lt.1.d+3) then
       tm=time ; ch=' ps '
    else if(time.lt.1.d+6) then
       tm=time*1.d-3 ; ch=' ns '
    else if(time.lt.1.d+9) then
       tm=time*1.d-6 ; ch=' mcs'
    else if(time.lt.1.d+12) then
       tm=time*1.d-9 ; ch=' ms '
    else if(time.lt.0.6d+14) then
       tm=time*1.d-12 ; ch='  s '
    else if(time.lt.3.6d+15) then
       tm=time*1.6667d-14 ; ch=' min'
    else
       tm=time*2.77778d-16 ; ch='  h '
    end if
    write(cha2,'(f16.4)') tm ; call posit(cha2,16,l0) 

!... to print on the right side
    j=j1-2*AJ(1,1)*scale
    write(1,'(a,2(i10,x),a)') '4 0 0 44 -1 26 '//cha1(l1:)//' 0.0000 4 225 390 ',&
             i,j,'Time: '//cha2(l0:)//ch//char(92)//'001'

!.... draw a box around the title if printing on the right side
    write(1,'(a)') '2 1 0 1 7 7 50 -1 20 0.000 0 0 -1 0 0 5'
    write(1,'(9x,10i10)') i-scale,j1+scale,i+10*scale,j1+scale,i+10*scale,j-scale,&
                           i-scale,j-scale,i-scale,j1+scale
  end subroutine title

  subroutine draw_temperature(Temperature,AJ)
    use param

    implicit none
    real*8 Temperature,AJ(2,2)
    integer i,j,l0,l1,j1
    character cha*10,cha1*5,cha2*16,ch*4

!... to print on the right side
    !i= ((L+1)*AJ(1,1)+2)*scale ; j1=0*AJ(2,2)*scale
    i= ((L+1)*AJ(1,1)+2)*scale ; j1=-3*AJ(2,2)*scale
    write(cha,'(f10.2)') Temperature ;     call posit(cha,10,l0)
    write(cha1,'(i5)') tsize    ;     call posit(cha1,5,l1)
    write(1,'(a,2(i10,x),a)') '4 0 0 44 -1 26 '//cha1(l1:)//' 0.0000 4 225 390 ',&
             i,j1,'Temperature '//cha(l0:)//' K'//char(92)//'001'
  end subroutine draw_temperature

  subroutine draw_proportions(nCis,nLtrans,nDtrans,nTotal,AJ)
    use param

    implicit none
    real*8    :: AJ(2,2)
    integer   :: i,j,l0,l1,j1
    integer   :: nCis,nLtrans,nDtrans,nTotal
    real*8    :: proportion
    character cha*10,cha1*5,cha2*16,ch*4

    ! Proportion of Cis
    i= ((L+1)*AJ(1,1)-28)*scale ; j1=2*AJ(2,2)*scale
    proportion = float(nCis)/float(nTotal) * 100
    write(cha,'(f10.2)') proportion ;     call posit(cha,10,l0)
    write(cha1,'(i5)') tsize    ;     call posit(cha1,5,l1)
    write(1,'(a,2(i10,x),a)') '4 0 12 44 -1 26 '//cha1(l1:)//' 0.0000 4 225 390 ',&
             i,j1,'Cis: '//cha(l0:)//char(92)//'001'

     ! Proportion of L-trans
     i= ((L+1)*AJ(1,1)-8)*scale ; j1=2*AJ(2,2)*scale
     proportion = float(nLtrans)/float(nTotal) * 100
     write(cha,'(f10.2)') proportion ;     call posit(cha,10,l0)
     write(cha1,'(i5)') tsize    ;     call posit(cha1,5,l1)
     write(1,'(a,2(i10,x),a)') '4 0 5 44 -1 26 '//cha1(l1:)//' 0.0000 4 225 390 ',&
              i,j1,'L-Trans: '//cha(l0:)//char(92)//'001'
 
    ! Proportion of D-trans
    i= ((L+1)*AJ(1,1)+8)*scale ; j1=2*AJ(2,2)*scale
    proportion = float(nDtrans)/float(nTotal) * 100
    write(cha,'(f10.2)') proportion ;     call posit(cha,10,l0)
    write(cha1,'(i5)') tsize    ;     call posit(cha1,5,l1)
    write(1,'(a,2(i10,x),a)') '4 0 1 44 -1 26 '//cha1(l1:)//' 0.0000 4 225 390 ',&
             i,j1,'D-Trans: '//cha(l0:)//char(92)//'001'

  end subroutine draw_proportions 



  subroutine number(x,y,ii)
!.............................................................................
!  write the number ii
!.............................................................................
    use param
    implicit none
    real*8 x,y
    integer i,j,l0,ii
    character cha*10
    i=x*scale-300/2. ; j=y*scale+225/2.
    write(cha,'(i10)') ii
    call posit(cha,10,l0)
    write(1,'(a,2(i10,x),a)')'4 0 0 44 -1 26 15 0.0000 4 225 390 ',&
             i,j,cha(l0:)//char(92)//'001'
  end subroutine number

  subroutine draw_line(x1,y1,x2,y2)
!.............................................................................
!  draw a line between points (x1,y1) and (x2,y2)
!.............................................................................
    use param
    implicit none
    real*8 x1,x2,y1,y2
    integer i1(2),i2(2)
    i1(1)=x1*scale ; i1(2)=y1*scale ;  i2(1)=x2*scale ; i2(2)=y2*scale 
    write(1,'(a)') '2 1 0 1 0 26 50 -1 20 0.000 0 0 -1 0 0 2'
    write(1,'(9x,8i10)') i1(1),i1(2),i2(1),i2(2)
  end subroutine draw_line

  subroutine cell(AL)
!.............................................................................
!  indicate the cell
!.............................................................................
    use param
    implicit none
    real*8 AL(2,2)
    integer i1(2),i2(2),i3(2),i4(2)
    i1(:)=0 ; i2(:)=AL(1,:)*scale ; i4(:)=AL(2,:)*scale ; i3(:)=i2(:)+i4(:)
    write(1,'(a)') '2 1 1 3 4 7 60 -1 -1 4.000 0 0 -1 0 0 5'
    write(1,'(10(i10))') i1(1),-i1(2),i2(1),-i2(2),i3(1),-i3(2),i4(1),-i4(2),i1(1),-i1(2)
  end subroutine cell

   subroutine circle(x,y)
!.............................................................................
!  draw a sphere centred at (x,y)
!.............................................................................
    use param
    implicit none
    real*8,parameter :: rad=0.1
    real*8 x,y
    integer i,j,k
    i=x*scale  ; j=y*scale  ; k=rad*scale
    write(1,'(a,8(i10))')'1 3 0 1 0 1 40 -1 20 0.000 1 0.0000 ',i,j,k,k,i,j,i+k,j
  end subroutine circle

  subroutine rectangular(x,y)
!.............................................................................
!  draw a rectangular centred at (x,y)
!.............................................................................
    use param
    implicit none
    real*8 x,y
    integer i1(2),i2(2),i3(2),i4(2)
    i1(1)=(x-0.2)*scale ; i1(2)=(y-0.4)*scale 
    i2(1)=(x+0.2)*scale ; i2(2)=(y-0.4)*scale 
    i3(1)=(x+0.2)*scale ; i3(2)=(y+0.4)*scale 
    i4(1)=(x-0.2)*scale ; i4(2)=(y+0.4)*scale 
    write(1,'(a)') '2 1 0 1 5 27 50 -1 20 0.000 0 0 -1 0 0 5'
    write(1,'(9x,10i10)') i1(1),i1(2),i2(1),i2(2),i3(1),i3(2),i4(1),i4(2),i1(1),i1(2)
  end subroutine rectangular

