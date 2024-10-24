      subroutine pbcdftd3(iz,lat,n,edisp,gdisp,sdisp,
     .              grad,method,echo,nversion,s6x,rs6x,s8x,rs8x,alp6,
     .              rthr,rthr2,rthr3,noabc)
      USE NUMBERS
      USE PARINF_MODULE
      USE BASATO_MODULE,ONLY:XA
      USE MEMORY_USE
      IMPLICIT REAL(FLOAT)(A-H,O-Z)
      PARAMETER(PARK1=16._FLOAT,PARK3=-4._FLOAT)
      PARAMETER(autokcal=627.509541_FLOAT)
      integer,parameter :: max_elem=94,maxc=5
      REAL(FLOAT),DIMENSION(:,:),ALLOCATABLE :: r0ab
      REAL(FLOAT),DIMENSION(:,:,:,:,:),ALLOCATABLE :: c6ab
      REAL(FLOAT),DIMENSION(MAX_ELEM) :: R2R4,RCOV,IDA
      REAL(FLOAT)::gdisp(3,n),LAT(3,3),sdisp(3,3),stress(3,3),dum_vec(3)
      INTEGER,DIMENSION(:),ALLOCATABLE :: mxc
      integer :: iz(n),z
      REAL(FLOAT) :: cn(n)
      REAL(FLOAT),optional :: s6x,s8x,rs6x,rs8x
      character(LEN=12) :: method,filename
      logical grad,echo,noabc
      integer rep_vdw(3),rep_cn(3)
      DATA (R2R4(I),I=1,72)/2.00734898_FLOAT,1.56637132_FLOAT,5.01986934
     *_FLOAT,3.85379032_FLOAT,3.64446594_FLOAT,3.10492822_FLOAT,2.711752
     *47_FLOAT,2.59361680_FLOAT,2.38825250_FLOAT,2.21522516_FLOAT,6.5858
     *5536_FLOAT,5.46295967_FLOAT,5.65216669_FLOAT,4.88284902_FLOAT,4.29
     *727576_FLOAT,4.04108902_FLOAT,3.72932356_FLOAT,3.44677275_FLOAT,
     *7.97762753_FLOAT,7.07623947_FLOAT,6.60844053_FLOAT,6.28791364_FLOA
     *T,6.07728703_FLOAT,5.54643096_FLOAT,5.80491167_FLOAT,5.58415602_FL
     *OAT,5.41374528_FLOAT,5.28497229_FLOAT,5.22592821_FLOAT,5.09817141_
     *FLOAT,6.12149689_FLOAT,5.54083734_FLOAT,5.06696878_FLOAT,4.8700510
     *8_FLOAT,4.59089647_FLOAT,4.31176304_FLOAT,9.55461698_FLOAT,8.67396
     *077_FLOAT,7.97210197_FLOAT,7.43439917_FLOAT,6.58711862_FLOAT,6.195
     *36215_FLOAT,6.01517290_FLOAT,5.81623410_FLOAT,5.65710424_FLOAT,5.5
     *2640661_FLOAT,5.44263305_FLOAT,5.58285373_FLOAT,7.02081898_FLOAT,6
     *.46815523_FLOAT,5.98089120_FLOAT,5.81686657_FLOAT,5.53321815_FLOAT
     *,5.25477007_FLOAT,11.02204549_FLOAT,10.15679528_FLOAT,9.35167836_F
     *LOAT,9.06926079_FLOAT,8.97241155_FLOAT,8.90092807_FLOAT,8.85984840
     *_FLOAT,8.81736827_FLOAT,8.79317710_FLOAT,7.89969626_FLOAT,8.805884
     *54_FLOAT,8.42439218_FLOAT,8.54289262_FLOAT,8.47583370_FLOAT,8.4509
     *0888_FLOAT,8.47339339_FLOAT,7.83525634_FLOAT,8.20702843_FLOAT/
      DATA (R2R4(I),I=73,94)/7.70559063_FLOAT,7.32755997_FLOAT,7.0388738
     *1_FLOAT,6.68978720_FLOAT,6.05450052_FLOAT,5.88752022_FLOAT,5.70661
     *499_FLOAT,5.78450695_FLOAT,7.79780729_FLOAT,7.26443867_FLOAT,6.781
     *51984_FLOAT,6.67883169_FLOAT,6.39024318_FLOAT,6.09527958_FLOAT,11.
     *79156076_FLOAT,11.10997644_FLOAT,9.51377795_FLOAT,8.67197068_FLOAT
     *,8.77140725_FLOAT,8.65402716_FLOAT,8.53923501_FLOAT,8.85024712_FLO
     *AT/
      DATA (RCOV(I),I=1,72)/0.80628308_FLOAT,1.15903197_FLOAT,3.02356173
     *_FLOAT,2.36845659_FLOAT,1.94011865_FLOAT,1.88972601_FLOAT,1.788940
     *56_FLOAT,1.58736983_FLOAT,1.61256616_FLOAT,1.68815527_FLOAT,3.5274
     *8848_FLOAT,3.14954334_FLOAT,2.84718717_FLOAT,2.62041997_FLOAT,2.77
     *159820_FLOAT,2.57002732_FLOAT,2.49443835_FLOAT,2.41884923_FLOAT,4.
     *43455700_FLOAT,3.88023730_FLOAT,3.35111422_FLOAT,3.07395437_FLOAT,
     *3.04875805_FLOAT,2.77159820_FLOAT,2.69600923_FLOAT,2.62041997_FLOA
     *T,2.51963467_FLOAT,2.49443835_FLOAT,2.54483100_FLOAT,2.74640188_FL
     *OAT,2.82199085_FLOAT,2.74640188_FLOAT,2.89757982_FLOAT,2.77159820_
     *FLOAT,2.87238349_FLOAT,2.94797246_FLOAT,4.76210950_FLOAT,4.2077898
     *0_FLOAT,3.70386304_FLOAT,3.50229216_FLOAT,3.32591790_FLOAT,3.12434
     *702_FLOAT,2.89757982_FLOAT,2.84718717_FLOAT,2.84718717_FLOAT,2.721
     *20556_FLOAT,2.89757982_FLOAT,3.09915070_FLOAT,3.22513231_FLOAT,3.1
     *7473967_FLOAT,3.17473967_FLOAT,3.09915070_FLOAT,3.32591790_FLOAT,3
     *.30072128_FLOAT,5.26603625_FLOAT,4.43455700_FLOAT,4.08180818_FLOAT
     *,3.70386304_FLOAT,3.98102289_FLOAT,3.95582657_FLOAT,3.93062995_FLO
     *AT,3.90543362_FLOAT,3.80464833_FLOAT,3.82984466_FLOAT,3.80464833_F
     *LOAT,3.77945201_FLOAT,3.75425569_FLOAT,3.75425569_FLOAT,3.72905937
     *_FLOAT,3.85504098_FLOAT,3.67866672_FLOAT,3.45189952_FLOAT/
      DATA (RCOV(I),I=73,94)/3.30072128_FLOAT,3.09915070_FLOAT,2.9731687
     *8_FLOAT,2.92277614_FLOAT,2.79679452_FLOAT,2.82199085_FLOAT,2.84718
     *717_FLOAT,3.32591790_FLOAT,3.27552496_FLOAT,3.27552496_FLOAT,3.426
     *70319_FLOAT,3.30072128_FLOAT,3.47709584_FLOAT,3.57788113_FLOAT,5.0
     *6446567_FLOAT,4.56053862_FLOAT,4.20778980_FLOAT,3.98102289_FLOAT,3
     *.82984466_FLOAT,3.85504098_FLOAT,3.88023730_FLOAT,3.90543362_FLOAT
     */
      if(echo) then
         
      endif


      AUTOEV=PAR(4)
      AUTOANG=PAR(32)
      CALL CRYALLOC(R0AB,MAX_ELEM,MAX_ELEM,'PBCDFTD3','R0AB')
      CALL CRYALLOC(C6AB,MAX_ELEM,MAX_ELEM,MAXC,MAXC,3,'PBCDFTD3','C6AB'
     *)
      CALL CRYALLOC(MXC,MAX_ELEM,'PBCDFTD3','MXC')
      call setr0ab(max_elem,r0ab)

      if(nversion.eq.2)then
          if(echo)write(IOUT,'(''loading DFT-D2 parameters ...'')')
          call loadoldpar(max_elem,maxc,c6ab,r0ab)
          mxc=1
          c6ab=c6ab*17.3452656005951_FLOAT
      else
        call copyc6(maxc,max_elem,c6ab,mxc)
      endif

      call set_criteria(rthr,lat,dum_vec)
      rep_vdw=NINT(dum_vec)+1
      call set_criteria(rthr2,lat,dum_vec)
      rep_cn=NINT(dum_vec)+1
      if (present(s6x).and.present(rs6x)
     *    .and.present(s8x).and.present(rs8x)) then
         s6=s6x
         rs6=rs6x
         s8=s8x
         rs8=rs8x
      else
           call setfuncpar(method,nversion,s6,rs6,s8,rs8,alp6)
      endif
      if (nversion.eq.4) then
         a1=rs6
         a2=rs8
      endif
      alp8=alp6+2
      
      if(lprint(127).eq.1) then
      if(.false.)then
         ida=0
         do i=1,n
            ida(iz(i))=ida(iz(i))+1
         end do
        write(IOUT,'('' C6 coefficients used:'')')
         do i=1,94
           if (ida(i).gt.0)then
           write(IOUT,*) mxc(i),' C6 for element ',i
            do j=1,maxc
            if (c6ab(i,i,j,j,1).gt.0)then
            write(IOUT,'('' Z='',i3,'' CN='',F6.3,5x,''C6(AA)='',F8.2)')
     .  i,c6ab(i,i,j,j,2),c6ab(i,i,j,j,1)
            end if
            end do
           end if
        end do
      endif
      if (echo) then
              write(IOUT,"(/'# Z ',
     .              ' R0(AA)[Ang.]',3x,
     .              'CN',7x,
     .              'C6(AA)     C8(AA)   C10(AA) [au] ') ")
      x=0
      call pbcncoord(n,rcov,iz,cn,lat,rep_cn,rthr2)

      do i=1,n
         z=iz(i)
         call getc6(maxc,max_elem,c6ab,mxc,iz(i),iz(i),cn(i),cn(i),c6)
         do j=1,n
          call getc6(maxc,max_elem,c6ab,mxc,iz(i),iz(j),cn(i),cn(j),dum)
          x=x+dum
         end do
         c8 =r2r4(iz(i))**2*3.0_FLOAT*c6
         c10=(49.0_FLOAT/40.0_FLOAT)*c8**2/c6
         dum=0.5_FLOAT*autoang*r0ab(z,z)
         if (version.eq.4)dum=rs6*0.5_FLOAT*autoang*sqrt(c8/c6)
         write(IOUT,'(i4,F7.3,4x,F7.3,3F12.1,L2)')
     .        iz(i),
     .        dum,cn(i),
     .        c6,c8,c10
      end do
      write(IOUT,'(/'' molecular/unit cell C6(AA) [au] = '',F12.2)')x
      endif
      endif

      call pbcedisp(max_elem,maxc,n,iz,c6ab,mxc,r2r4,r0ab,
     .           rcov,rs6,rs8,alp6,alp8,nversion,noabc,
     .           e6,e8,e6abc,lat,rthr,rep_vdw,rthr2,rep_cn,rthr3)
      
      e6   = e6   *s6
      e6abc= e6abc*s6
      e8   = e8   *s8
      
      edisp =-e6-e8-e6abc

      
      if(grad)THEN
      call pbcgdisp(max_elem,maxc,n,iz,c6ab,mxc,r2r4,r0ab,
     *rcov,s6,s8,rs6,rs8,alp6,alp8,noabc,nversion,gdisp,
     *edisp_check,gnorm,sdisp,lat,rep_vdw,rep_cn,rthr,rthr2,rthr3)     
     
      EDISP=EDISP_CHECK
      ENDIF
      if (echo) then
      if(nversion.lt.4)then
      write(IOUT,'(/10x,'' DFT-D V'',i1)') nversion
      else
      write(IOUT,'(/10x,'' DFT-D V3(BJ)'')')
      endif
      write(IOUT,'('' DF '',a50)') method
      write(IOUT,'('' parameters'')')
      if(nversion.eq.2)then
         write(IOUT,'('' s6       :'',f10.4)') s6
         write(IOUT,'('' alpha6   :'',f10.4)') alp6
      endif
      if(nversion.eq.3)then
         write(IOUT,'('' s6       :'',f10.4)') s6
         write(IOUT,'('' s8       :'',f10.4)') s8
         write(IOUT,'('' rs6      :'',f10.4)') rs6
         write(IOUT,'('' rs18     :'',f10.4)') rs8
         write(IOUT,'('' alpha6   :'',f10.4)') alp6
         write(IOUT,'('' alpha8   :'',f10.4)') alp8
         write(IOUT,'('' k1 k3    :'',2f10.4)')PARK1,PARK3
      endif
      if(nversion.eq.4)then
         write(IOUT,'('' s6       :'',f10.4)') s6
         write(IOUT,'('' s8       :'',f10.4)') s8
         write(IOUT,'('' a1       :'',f10.4)') rs6
         write(IOUT,'('' a2       :'',f10.4)') rs8
         write(IOUT,'('' k1 k3    :'',2f10.4)')PARK1,PARK3
      endif
      write(IOUT,'('' Cutoff   :'',f10.4)') sqrt(rthr)*autoang
      write(IOUT,'('' CN-Cutoff:'',f10.4)') sqrt(rthr2)*autoang
      write(IOUT,'('' ABC-Cutoff:'',f10.4)') sqrt(rthr3)*autoang
      write(IOUT,*)
      write(IOUT,'('' Edisp /kcal,au:'',f11.4,f12.8)')edisp*autokcal,
     *edisp
      write(IOUT,'(/'' E6    /kcal :'',f11.4)')-e6*autokcal
      if(nversion.gt.2)then
        write(IOUT,'('' E8    /kcal :'',f11.4)')-e8*autokcal
        if(.not.noabc)
     .  write(IOUT,'('' E6(ABC) "   :'',2f11.6,F16.12)')-e6abc*autokcal
      endif !version.gt.2      
      write(IOUT,'('' '')')
      endif             !echo
      CALL CRYDEALLOC(MXC,'PBCDFTD3','MXC')
      CALL CRYDEALLOC(C6AB,'PBCDFTD3','C6AB')
      CALL CRYDEALLOC(R0AB,'PBCDFTD3','R0AB')
      end subroutine pbcdftd3

      subroutine pbcedisp(max_elem,maxc,n,iz,c6ab,mxc,r2r4,r0ab,
     .           rcov,rs6,rs8,alp6,alp8,nversion,noabc,
     .           e6,e8,e63,lat,rthr,rep_vdw,cn_thr,rep_cn,rthr3)
      USE NUMBERS
      USE PARINF_MODULE
      USE BASATO_MODULE,ONLY:XA
      USE MEMORY_USE
      IMPLICIT REAL(FLOAT)(A-H,O-Z)
      REAL(FLOAT),DIMENSION(:),ALLOCATABLE :: CN,R2AB,CC6AB
      INTEGER,DIMENSION(:),ALLOCATABLE :: ICOMP
      DIMENSION R0AB(MAX_ELEM,MAX_ELEM),R2R4(MAX_ELEM),RCOV(MAX_ELEM)
      DIMENSION C6AB(MAX_ELEM,MAX_ELEM,MAXC,MAXC,3),MXC(MAX_ELEM),IZ(N)
      REAL(FLOAT),DIMENSION(3,3) :: lat
      REAL(FLOAT),DIMENSION(3) :: rxyz,dxyz,d2,tau
      integer rep_vdw(3),rep_cn(3)
      integer taux,tauy,tauz,counter
      logical noabc
      e6 =0._FLOAT
      e8 =0._FLOAT
      e63=0._FLOAT
      tau=(/0._FLOAT,0._FLOAT,0._FLOAT/)
      counter=0
      crit_cn=cn_thr
      a1=rs6
      a2=rs8      
      I=N*N
      CALL CRYALLOC(CN,N,'PBCEDISP','CN')
      CALL CRYALLOC(R2AB,I,'PBCEDISP','R2AB')
      CALL CRYALLOC(CC6AB,I,'PBCEDISP','CC6AB')
      CALL CRYALLOC(ICOMP,I,'PBCEDISP','ICOMP')
      if(nversion.eq.2)then

      do iat=1,n-1
         do jat=iat+1,n
           c6=c6ab(iz(jat),iz(iat),1,1,1)
           do taux=-rep_vdw(1),rep_vdw(1)
           do tauy=-rep_vdw(2),rep_vdw(2)
           do tauz=-rep_vdw(3),rep_vdw(3)
            tau=taux*lat(:,1)+tauy*lat(:,2)+tauz*lat(:,3)
            dx=XA(1,iat)-XA(1,jat)+tau(1)
            dy=XA(2,iat)-XA(2,jat)+tau(2)
            dz=XA(3,iat)-XA(3,jat)+tau(3)
            r2=dx*dx+dy*dy+dz*dz
           if(r2.gt.rthr) cycle
            r=sqrt(r2)
            damp6=1._FLOAT/(1._FLOAT+exp(-alp6*(r/(rs6*r0ab(iz(jat),
     *iz(iat)))-1._FLOAT)))
            r6=r2**3      
            e6 =e6+c6*damp6/r6
            counter=counter+1
           enddo !taux
           enddo !tauy
           enddo !tauz
         enddo
      enddo
      
      do iat=1,n
        jat=iat
        c6=c6ab(iz(jat),iz(iat),1,1,1)
        do taux=-rep_vdw(1),rep_vdw(1)
        do tauy=-rep_vdw(2),rep_vdw(2)
        do tauz=-rep_vdw(3),rep_vdw(3)
          if (taux.eq.0 .and. tauy.eq.0 .and. tauz.eq.0) cycle
          tau=taux*lat(:,1)+tauy*lat(:,2)+tauz*lat(:,3)
          dx=tau(1)
          dy=tau(2)
          dz=tau(3)
          r2=dx*dx+dy*dy+dz*dz
           if(r2.gt.rthr) cycle
          r=sqrt(r2)
          damp6=1._FLOAT/(1._FLOAT+exp(-alp6*(r/(rs6*r0ab(iz(jat),
     *iz(iat)))-1._FLOAT)))
          r6=r2**3      
          e6 =e6+c6*damp6/r6*0.5_FLOAT
          counter=counter+1
        enddo
        enddo
        enddo
      enddo !iat
!      write(*,*)'counter: ',counter
      
      

      else if (nversion.eq.3) then

      call pbcncoord(n,rcov,iz,cn,lat,rep_cn,crit_cn)

      icomp=0
      do iat=1,n-1
         do jat=iat+1,n
          call getc6(maxc,max_elem,c6ab,mxc,iz(iat),iz(jat),
     .                                  cn(iat),cn(jat),c6)

           do taux=-rep_vdw(1),rep_vdw(1)
           do tauy=-rep_vdw(2),rep_vdw(2)
           do tauz=-rep_vdw(3),rep_vdw(3)
            tau=taux*lat(:,1)+tauy*lat(:,2)+tauz*lat(:,3)

            dx=XA(1,iat)-XA(1,jat)+tau(1)
            dy=XA(2,iat)-XA(2,jat)+tau(2)
            dz=XA(3,iat)-XA(3,jat)+tau(3)
            r2=dx*dx+dy*dy+dz*dz

           if(r2.gt.rthr) cycle
            r =sqrt(r2)
            rr=r0ab(iz(jat),iz(iat))/r
            tmp=rs6*rr   
            damp6 =1._FLOAT/(1._FLOAT+6._FLOAT*tmp**alp6)
            tmp=rs8*rr     
            damp8 =1._FLOAT/(1._FLOAT+6._FLOAT*tmp**alp8)

            if(.not.noabc)then
              ij=lin(jat,iat)
              icomp(ij)=1
              cc6ab(ij)=sqrt(c6)
            endif

            r6=r2**3      
            e6 =e6+c6*damp6/r6
            c8 =3._FLOAT*c6*r2r4(iz(iat))*r2r4(iz(jat))
            r8 =r6*r2

            e8 =e8+c8*damp8/r8

            counter=counter+1

           enddo !tauz
           enddo !tauy
           enddo !taux
         enddo !jat
      enddo !iat
      
      do iat=1,n
        jat=iat
        call getc6(maxc,max_elem,c6ab,mxc,iz(iat),iz(jat),
     .                                  cn(iat),cn(jat),c6)
         
        do taux=-rep_vdw(1),rep_vdw(1)
         do tauy=-rep_vdw(2),rep_vdw(2)
          do tauz=-rep_vdw(3),rep_vdw(3)
            if (taux.eq.0 .and. tauy.eq.0 .and. tauz.eq.0) cycle
            tau=taux*lat(:,1)+tauy*lat(:,2)+tauz*lat(:,3)

            dx=tau(1)
            dy=tau(2)
            dz=tau(3)
            r2=dx*dx+dy*dy+dz*dz
           if(r2.gt.rthr) cycle
            r =sqrt(r2)
            rr=r0ab(iz(jat),iz(iat))/r
            tmp=rs6*rr
            damp6 =1._FLOAT/(1._FLOAT+6._FLOAT*tmp**alp6)
            tmp=rs8*rr
            damp8 =1._FLOAT/(1._FLOAT+6._FLOAT*tmp**alp8)

            if(.not.noabc)then
              ij=lin(jat,iat)
              icomp(ij)=1
              cc6ab(ij)=sqrt(c6)
            endif

            r6=r2**3

            e6 =e6+c6*damp6/r6*0.5_FLOAT

            c8 =3._FLOAT*c6*r2r4(iz(iat))*r2r4(iz(jat))
            r8 =r6*r2

            e8 =e8+c8*damp8/r8*0.5_FLOAT
            counter=counter+1

         enddo !tauz
        enddo !tauy
       enddo !taux
      enddo !iat
      else if (nversion.eq.4) then


      call pbcncoord(n,rcov,iz,cn,lat,rep_cn,crit_cn)
      icomp=0
      do iat=1,n
         do jat=iat+1,n
           call getc6(maxc,max_elem,c6ab,mxc,iz(iat),iz(jat),
     .                                  cn(iat),cn(jat),c6)

           rxyz=XA(:,iat)-XA(:,jat)
           r42=r2r4(iz(iat))*r2r4(iz(jat))
           bj_dmp6=(a1*sqrt(3.0_FLOAT*r42)+a2)**6
           bj_dmp8=(a1*sqrt(3.0_FLOAT*r42)+a2)**8

           do taux=-rep_vdw(1),rep_vdw(1)
           do tauy=-rep_vdw(2),rep_vdw(2)
           do tauz=-rep_vdw(3),rep_vdw(3)
            tau=taux*lat(:,1)+tauy*lat(:,2)+tauz*lat(:,3)
            
            dxyz=rxyz+tau

            r2=sum(dxyz*dxyz)
           if(r2.gt.rthr) cycle
            r =sqrt(r2)
            rr=r0ab(iz(jat),iz(iat))/r

            if(.not.noabc)then
              ij=lin(jat,iat)
              icomp(ij)=1
              cc6ab(ij)=sqrt(c6)
            endif

            r6=r2**3      

            e6 =e6+c6/(r6+bj_dmp6)
            c8 =3._FLOAT*c6*r42
            r8 =r6*r2

            e8 =e8+c8/(r8+bj_dmp8)

            counter=counter+1

           enddo !tauz
           enddo !tauy
           enddo !taux
         enddo !jat

        jat=iat
        call getc6(maxc,max_elem,c6ab,mxc,iz(iat),iz(jat),
     .                                  cn(iat),cn(jat),c6)
        r42=r2r4(iz(iat))*r2r4(iz(iat))
        bj_dmp6=(a1*sqrt(3.0_FLOAT*r42)+a2)**6
        bj_dmp8=(a1*sqrt(3.0_FLOAT*r42)+a2)**8
         
        do taux=-rep_vdw(1),rep_vdw(1)
         do tauy=-rep_vdw(2),rep_vdw(2)
          do tauz=-rep_vdw(3),rep_vdw(3)
            if (taux.eq.0 .and. tauy.eq.0 .and. tauz.eq.0) cycle
            tau=taux*lat(:,1)+tauy*lat(:,2)+tauz*lat(:,3)

            r2=sum(tau*tau)
           if(r2.gt.rthr) cycle
            r =sqrt(r2)
            rr=r0ab(iz(jat),iz(iat))/r

            if(.not.noabc)then
              ij=lin(jat,iat)
              icomp(ij)=1
              cc6ab(ij)=sqrt(c6)
            endif

            r6=r2**3      

            e6 =e6+c6/(r6+bj_dmp6)*0.5_FLOAT

            c8 =3._FLOAT*c6*r42
            r8 =r6*r2

            e8 =e8+c8/(r8+bj_dmp8)*0.5_FLOAT
            counter=counter+1


         enddo !tauz
        enddo !tauy
       enddo !taux
      enddo !iat
      endif !nversion
      if(noabc)return
      call pbcthreebody(max_elem,lat,n,iz,rep_cn,rthr3,cc6ab,r0ab,e63)
      CALL CRYDEALLOC(ICOMP,'PBCEDISP','ICOMP')
      CALL CRYDEALLOC(CC6AB,'PBCEDISP','CC6AB')
      CALL CRYDEALLOC(R2AB,'PBCEDISP','R2AB')
      CALL CRYDEALLOC(CN,'PBCEDISP','CN')
      end subroutine pbcedisp
      SUBROUTINE pbcthreebody(max_elem,lat,n,iz,repv,abcthr,cc6ab,r0ab,
     *eabc)
      USE NUMBERS
      USE PARINF_MODULE
      USE BASATO_MODULE,ONLY:XA
      USE PARAL1_MODULE
      IMPLICIT REAL(FLOAT)(A-H,O-Z)
      PARAMETER(SR9=.75_FLOAT,ALP9=-16._FLOAT) !reciprocal radii scaling parameter for damping function (s_r=4/3)
      PARAMETER(THIRD=1._FLOAT/3._FLOAT,SIXTH=THIRD*.5_FLOAT)
      REAL(FLOAT),INTENT(OUT)::eabc
      REAL(FLOAT),DIMENSION(3,3),INTENT(IN)::lat
      REAL(FLOAT),DIMENSION(3):: jtau,ktau,jxyz,kxyz,ijvec,ikvec,jkvec
      REAL(FLOAT),DIMENSION(3):: TUX0,TUY0,TUZ0,TUX,TUY
      REAL(FLOAT),DIMENSION(3):: TUY0K,TUZ0K,TUXK,TUYK
      REAL(FLOAT),INTENT(IN) :: cc6ab(n*n),R0AB(MAX_ELEM,MAX_ELEM)
      integer,intent(in) :: IZ(N)
      INTEGER,DIMENSION(3):: repv
      INTEGER,EXTERNAL :: lin
      eabc=0._FLOAT
      IREPV1=REPV(1)
      IREPV2=REPV(2)
      IREPV3=REPV(3)
      MREPV1=-IREPV1
      MREPV2=-IREPV2
      MREPV3=-IREPV3
      TUX0(1:3)=LAT(1:3,1)*REAL(MREPV1,FLOAT)
      TUY0(1:3)=LAT(1:3,2)*REAL(MREPV2,FLOAT)
      TUZ0(1:3)=LAT(1:3,3)*REAL(MREPV3,FLOAT)
      CALL IGPVAL
      do iat=3,n
      IZIAT=IZ(IAT)
        do jat=2,iat-1
      ITASK=ITASK+1
      IF(MPOINT.NE.ITASK)CYCLE
         IZJAT=IZ(JAT)
         ijvec(1:3)=XA(1:3,jat)-XA(1:3,iat)
         ij=lin(iat,jat)
         R0IJ=1._FLOAT/R0AB(IZIAT,IZJAT)
          do kat=1,jat-1
            ik=lin(iat,kat)
            jk=lin(jat,kat)
            ikvec(1:3)=XA(1:3,kat)-XA(1:3,iat)
            jkvec(1:3)=XA(1:3,kat)-XA(1:3,jat)
            c9=-1.0_FLOAT*(cc6ab(ij)*cc6ab(ik)*cc6ab(jk))
            IZKAT=IZ(KAT)
            R0IK=1._FLOAT/R0AB(IZIAT,IZKAT)
            R0JK=1._FLOAT/R0AB(IZJAT,IZKAT)
            TUX(1:3)=TUX0(1:3)
            do jtaux=MREPV1,IREPV1
            IREPMIN1=MAX(MREPV1,JTAUX+MREPV1)
            IREPMAX1=MIN(IREPV1,JTAUX+IREPV1)
            TUY(1:3)=TUY0(1:3)+TUX(1:3)
            do jtauy=MREPV2,IREPV2
            IREPMIN2=MAX(MREPV2,JTAUY+MREPV2)
            IREPMAX2=MIN(IREPV2,JTAUY+IREPV2)
            JTAU(1:3)=TUZ0(1:3)+TUY(1:3)
            do jtauz=MREPV3,IREPV3
              IREPMIN3=MAX(MREPV3,JTAUZ+MREPV3)
              IREPMAX3=MIN(IREPV3,JTAUZ+IREPV3)
              rij2=SUM((ijvec(1:3)+JTAU(1:3))**2)
              IF(RIJ2.GT.ABCTHR)THEN
              JTAU(1:3)=LAT(1:3,3)+JTAU(1:3)
              CYCLE
              ENDIF
              RR0IJ=SQRT(rij2)*R0IJ
              TUXK(1:3)=LAT(1:3,1)*REAL(IREPMIN1,FLOAT)
              TUY0K(1:3)=LAT(1:3,2)*REAL(IREPMIN2,FLOAT)
              TUZ0K(1:3)=LAT(1:3,3)*REAL(IREPMIN3,FLOAT)
              do ktaux=IREPMIN1,IREPMAX1
              TUYK(1:3)=TUY0K(1:3)+TUXK(1:3)
              do ktauy=IREPMIN2,IREPMAX2
              KTAU(1:3)=TUZ0K(1:3)+TUYK(1:3)
              do ktauz=IREPMIN3,IREPMAX3
                rik2=SUM((IKVEC(1:3)+KTAU(1:3))**2)
              IF(RIK2.GT.ABCTHR)THEN
              KTAU(1:3)=LAT(1:3,3)+KTAU(1:3)
              CYCLE
              ENDIF
                rjk2=SUM((JKVEC(1:3)+KTAU(1:3)-JTAU(1:3))**2)
              IF(RJK2.GT.ABCTHR)THEN
              KTAU(1:3)=LAT(1:3,3)+KTAU(1:3)
              CYCLE
              ENDIF
              geomean=(SQRT(RIK2*RJK2)*R0IK*R0JK*rr0ij)**THIRD
                fdamp=1._FLOAT/(1._FLOAT+6._FLOAT*(sr9*geomean)**alp9)
                TMP=1._FLOAT/(RIJ2*RJK2*RIK2)
      ANG=((RIJ2+rjk2-RIK2)*(RIJ2+RIK2-RJK2)*(RIK2+RJK2-RIJ2)*TMP*
     *0.375_FLOAT+1._FLOAT)*TMP**1.5_FLOAT
                EABC=fdamp*C9*ANG+EABC
              KTAU(1:3)=LAT(1:3,3)+KTAU(1:3)
              ENDDO !ktauz
              TUYK(1:3)=LAT(1:3,2)+TUYK(1:3)
              ENDDO !ktauy
              TUXK(1:3)=LAT(1:3,1)+TUXK(1:3)
              ENDDO !ktaux
            JTAU(1:3)=LAT(1:3,3)+JTAU(1:3)
            ENDDO !jtauz
            TUY(1:3)=LAT(1:3,2)+TUY(1:3)
            ENDDO !jtauy
            TUX(1:3)=LAT(1:3,1)+TUX(1:3)
            ENDDO !jtaux
          ENDDO !kat
      CALL IGPVAL
        ENDDO !jat
      ENDDO !iat
      DO iat=2,n
      ij=lin(IAT,IAT)
      IZIAT=IZ(IAT)
      R0IJ=1._FLOAT/R0AB(IZIAT,IZIAT)
        DO kat=1,iat-1
      ITASK=ITASK+1
      IF(MPOINT.NE.ITASK)CYCLE
          ik=lin(iat,kat)
          ikvec(1:3)=XA(1:3,kat)-XA(1:3,iat)
          c9=-(cc6ab(ij)*cc6ab(ik)*cc6ab(IK))
      IZKAT=IZ(KAT)
      R0IK=1._FLOAT/R0AB(IZIAT,IZKAT)
      TUX(1:3)=TUX0(1:3)
          do jtaux=MREPV1,IREPV1
          IREPMIN1=MAX(MREPV1,JTAUX+MREPV1)
          IREPMAX1=MIN(IREPV1,JTAUX+IREPV1)
          TUY(1:3)=TUY0(1:3)+TUX(1:3)
          do jtauy=MREPV2,IREPV2
          IREPMIN2=MAX(MREPV2,JTAUY+MREPV2)
          IREPMAX2=MIN(IREPV2,JTAUY+IREPV2)
          JTAU(1:3)=TUZ0(1:3)+TUY(1:3)
          do jtauz=MREPV3,IREPV3
          IREPMIN3=MAX(MREPV3,JTAUZ+MREPV3)
          IREPMAX3=MIN(IREPV3,JTAUZ+IREPV3)
          IF(ABS(jtaux)+ABS(jtauy)+ABS(jtauz).eq.0)THEN
          JTAU(1:3)=LAT(1:3,3)+JTAU(1:3)
          CYCLE
          ENDIF
          rij2=SUM(JTAU(1:3)**2)
            IF(RIJ2.GT.ABCTHR)THEN
            JTAU(1:3)=LAT(1:3,3)+JTAU(1:3)
            CYCLE
            ENDIF
            RR0IJ=SQRT(rij2)*R0IJ
            TUXK(1:3)=LAT(1:3,1)*REAL(IREPMIN1,FLOAT)
            TUY0K(1:3)=LAT(1:3,2)*REAL(IREPMIN2,FLOAT)
            TUZ0K(1:3)=LAT(1:3,3)*REAL(IREPMIN3,FLOAT)
            do ktaux=IREPMIN1,IREPMAX1
            TUYK(1:3)=TUY0K(1:3)+TUXK(1:3)
            do ktauy=IREPMIN2,IREPMAX2
            KTAU(1:3)=TUZ0K(1:3)+TUYK(1:3)
            do ktauz=IREPMIN3,IREPMAX3
            rik2=SUM((IKVEC(1:3)+KTAU(1:3))**2)
            IF(RIK2.GT.ABCTHR)THEN
            KTAU(1:3)=LAT(1:3,3)+KTAU(1:3)
            CYCLE
            ENDIF
            rjk2=SUM((IKVEC(1:3)+KTAU(1:3)-JTAU(1:3))**2)
            IF(RJK2.GT.ABCTHR)THEN
            KTAU(1:3)=LAT(1:3,3)+KTAU(1:3)
            CYCLE
            ENDIF
            geomean=(SQRT(RIK2*rjk2)*R0IK*R0IK*rr0ij)**THIRD
            fdamp=1._FLOAT/(1._FLOAT+6._FLOAT*(sr9*geomean)**alp9)
            TMP=1._FLOAT/(RIJ2*RJK2*RIK2)
      ANG=((RIJ2+rjk2-RIK2)*(RIJ2+RIK2-RJK2)*(RIK2+RJK2-RIJ2)*TMP*
     *0.375_FLOAT+1._FLOAT)*TMP**1.5_FLOAT
                EABC=fdamp*C9*ANG*0.5_FLOAT+EABC
            KTAU(1:3)=LAT(1:3,3)+KTAU(1:3)
            ENDDO !ktauz
            TUYK(1:3)=LAT(1:3,2)+TUYK(1:3)
            ENDDO !ktauy
            TUXK(1:3)=LAT(1:3,1)+TUXK(1:3)
            ENDDO !ktaux
          JTAU(1:3)=LAT(1:3,3)+JTAU(1:3)
          ENDDO !jtauz
          TUY(1:3)=LAT(1:3,2)+TUY(1:3)
          ENDDO !jtauy
          TUX(1:3)=LAT(1:3,1)+TUX(1:3)
          ENDDO !jtaux
      CALL IGPVAL
        ENDDO !kat
      ENDDO !iat
      DO iat=2,n
      IZIAT=IZ(IAT)
        DO JAT=1,iat-1
      ITASK=ITASK+1
      IF(MPOINT.NE.ITASK)CYCLE
          ij=lin(iat,JAT)
          jk=lin(JAT,JAT)
          IJVEC(1:3)=XA(1:3,JAT)-XA(1:3,iat)
          c9=-(cc6ab(ij)*cc6ab(IJ)*cc6ab(jk))
      IZJAT=IZ(JAT)
      R0IJ=1._FLOAT/R0AB(IZIAT,IZJAT)
      R0JK=1._FLOAT/R0AB(IZJAT,IZJAT)
      TUX(1:3)=TUX0(1:3)
          do jtaux=MREPV1,IREPV1
          IREPMIN1=MAX(MREPV1,JTAUX+MREPV1)
          IREPMAX1=MIN(IREPV1,JTAUX+IREPV1)
          TUY(1:3)=TUY0(1:3)+TUX(1:3)
          do jtauy=MREPV2,IREPV2
          IREPMIN2=MAX(MREPV2,JTAUY+MREPV2)
          IREPMAX2=MIN(IREPV2,JTAUY+IREPV2)
          JTAU(1:3)=TUZ0(1:3)+TUY(1:3)
          do jtauz=MREPV3,IREPV3
          IREPMIN3=MAX(MREPV3,JTAUZ+MREPV3)
          IREPMAX3=MIN(IREPV3,JTAUZ+IREPV3)
            rij2=SUM((IJVEC(1:3)+JTAU(1:3))**2)
            IF(RIJ2.GT.ABCTHR)THEN
            JTAU(1:3)=LAT(1:3,3)+JTAU(1:3)
            CYCLE
            ENDIF
            RR0IJ=SQRT(rij2)*R0IJ
            TUXK(1:3)=LAT(1:3,1)*REAL(IREPMIN1,FLOAT)
            TUY0K(1:3)=LAT(1:3,2)*REAL(IREPMIN2,FLOAT)
            TUZ0K(1:3)=LAT(1:3,3)*REAL(IREPMIN3,FLOAT)
            do ktaux=IREPMIN1,IREPMAX1
            TUYK(1:3)=TUY0K(1:3)+TUXK(1:3)
            do ktauy=IREPMIN2,IREPMAX2
            KTAU(1:3)=TUZ0K(1:3)+TUYK(1:3)
            do ktauz=IREPMIN3,IREPMAX3
      IF(ABS(jtaux-KTAUX)+ABS(jtauy-KTAUY)+ABS(jtauz-KTAUZ).EQ.0)THEN
            KTAU(1:3)=LAT(1:3,3)+KTAU(1:3)
            CYCLE
            ENDIF
            rik2=SUM((IJVEC(1:3)+KTAU(1:3))**2)
            IF(RIK2.GT.ABCTHR)THEN
            KTAU(1:3)=LAT(1:3,3)+KTAU(1:3)
            CYCLE
            ENDIF
            rjk2=SUM((KTAU(1:3)-JTAU(1:3))**2)
            IF(RJK2.GT.ABCTHR)THEN 
            KTAU(1:3)=LAT(1:3,3)+KTAU(1:3)
            CYCLE
            ENDIF
            geomean=(SQRT(RIK2*rjk2)*R0IJ*R0JK*rr0ij)**THIRD
            fdamp=1._FLOAT/(1._FLOAT+6._FLOAT*(sr9*geomean)**alp9)
            TMP=1._FLOAT/(RIJ2*RJK2*RIK2)
      ANG=((RIJ2+rjk2-RIK2)*(RIJ2+RIK2-RJK2)*(RIK2+RJK2-RIJ2)*TMP*
     *0.375_FLOAT+1._FLOAT)*TMP**1.5_FLOAT
              EABC=fdamp*C9*ANG*0.5_FLOAT+EABC
            KTAU(1:3)=LAT(1:3,3)+KTAU(1:3)
            ENDDO !ktauz
            TUYK(1:3)=LAT(1:3,2)+TUYK(1:3)
            ENDDO !ktauy
            TUXK(1:3)=LAT(1:3,1)+TUXK(1:3)
            ENDDO !ktaux
          JTAU(1:3)=LAT(1:3,3)+JTAU(1:3)
          ENDDO !jtauz
          TUY(1:3)=LAT(1:3,2)+TUY(1:3)
          ENDDO !jtauy
          TUX(1:3)=LAT(1:3,1)+TUX(1:3)
          ENDDO !jtaux
      CALL IGPVAL
        ENDDO !kat
      ENDDO !iat
      DO iat=1,n
      ITASK=ITASK+1
      IF(MPOINT.NE.ITASK)CYCLE
      ij=lin(iat,IAT)
      c9=-cc6ab(ij)**3
      IZIAT=IZ(IAT)
      R0IJ=1._FLOAT/R0AB(IZIAT,IZIAT)
      TUX(1:3)=TUX0(1:3)
        do jtaux=MREPV1,IREPV1
        IREPMIN1=MAX(MREPV1,JTAUX+MREPV1)
        IREPMAX1=MIN(IREPV1,JTAUX+IREPV1)
        TUY(1:3)=TUY0(1:3)+TUX(1:3)
        do jtauy=MREPV2,IREPV2
        IREPMIN2=MAX(MREPV2,JTAUY+MREPV2)
        IREPMAX2=MIN(IREPV2,JTAUY+IREPV2)
        JTAU(1:3)=TUZ0(1:3)+TUY(1:3)
        do jtauz=MREPV3,IREPV3
        IREPMIN3=MAX(MREPV3,JTAUZ+MREPV3)
        IREPMAX3=MIN(IREPV3,JTAUZ+IREPV3)
        IF(ABS(jtaux)+ABS(jtauy)+ABS(jtauz).eq.0)THEN
        JTAU(1:3)=LAT(1:3,3)+JTAU(1:3)
        CYCLE
        ENDIF
        rij2=SUM(JTAU(1:3)**2)
        IF(RIJ2.GT.ABCTHR)THEN
        JTAU(1:3)=LAT(1:3,3)+JTAU(1:3)
        CYCLE
        ENDIF
        RR0IJ=SQRT(rij2)*R0IJ
        TUXK(1:3)=LAT(1:3,1)*REAL(IREPMIN1,FLOAT)
        TUY0K(1:3)=LAT(1:3,2)*REAL(IREPMIN2,FLOAT)
        TUZ0K(1:3)=LAT(1:3,3)*REAL(IREPMIN3,FLOAT)
          do ktaux=IREPMIN1,IREPMAX1
          TUYK(1:3)=TUY0K(1:3)+TUXK(1:3)
          do ktauy=IREPMIN2,IREPMAX2
          KTAU(1:3)=TUZ0K(1:3)+TUYK(1:3)
          do ktauz=IREPMIN3,IREPMAX3
          IF(ABS(ktaux)+ABS(ktauy)+ABS(ktauz).eq.0.OR.ABS(ktaux-jtaux)+
     *ABS(ktauy-jtauy)+ABS(ktauz-jtauz).EQ.0)THEN
          KTAU(1:3)=LAT(1:3,3)+KTAU(1:3)
          CYCLE
          ENDIF
          rik2=SUM(KTAU(1:3)**2)
          IF(RIK2.GT.ABCTHR)THEN
          KTAU(1:3)=LAT(1:3,3)+KTAU(1:3)
          CYCLE
          ENDIF
          rjk2=SUM((KTAU(1:3)-JTAU(1:3))**2)
          IF(RJK2.GT.ABCTHR)THEN
          KTAU(1:3)=LAT(1:3,3)+KTAU(1:3)
          CYCLE
          ENDIF
          geomean=(SQRT(RIK2*RJK2)*R0IJ*R0IJ*rr0ij)**THIRD
          fdamp=1._FLOAT/(1._FLOAT+6._FLOAT*(sr9*geomean)**alp9)
          TMP=1._FLOAT/(RIJ2*RJK2*RIK2)
      ANG=((RIJ2+rjk2-RIK2)*(RIJ2+RIK2-RJK2)*(RIK2+RJK2-RIJ2)*TMP*
     *0.375_FLOAT+1._FLOAT)*TMP**1.5_FLOAT
          EABC=fdamp*C9*ANG*SIXTH+EABC
          KTAU(1:3)=LAT(1:3,3)+KTAU(1:3)
          ENDDO !ktauz
          TUYK(1:3)=LAT(1:3,2)+TUYK(1:3)
          ENDDO !ktauy
          TUXK(1:3)=LAT(1:3,1)+TUXK(1:3)
          ENDDO !ktaux
        JTAU(1:3)=LAT(1:3,3)+JTAU(1:3)
        ENDDO !jtauz
        TUY(1:3)=LAT(1:3,2)+TUY(1:3)
        ENDDO !jtauy
        TUX(1:3)=LAT(1:3,1)+TUX(1:3)
        ENDDO !jtaux
      CALL IGPVAL
      ENDDO !iat
      CALL IGPRST
      CALL GSUM(EABC,1)
      END SUBROUTINE pbcthreebody
      
      subroutine pbcgdisp(max_elem,maxc,n,iz,c6ab,mxc,r2r4,r0ab,
     *            rcov,s6,s18,rs6,rs8,alp6,alp8,noabc,
     *                 nversion,g,disp,gnorm,stress,lat,rep_v,rep_cn,
     *                 crit_vdw,crit_cn,crit_abc)
      USE NUMBERS
      USE PARINF_MODULE
      USE BASATO_MODULE,ONLY:XA
      USE MEMORY_USE
      USE PARAL1_MODULE
      IMPLICIT REAL(FLOAT)(A-H,O-Z)
      PARAMETER(SR9=.75_FLOAT,ALP9=-16._FLOAT,AUTOKCAL=627.509541_FLOAT)
      PARAMETER(PARK1=-16._FLOAT,THIRD=1._FLOAT/3._FLOAT,HALF=.5_FLOAT)
      PARAMETER(SIXTH=THIRD*HALF,ALP92=ALP9+ALP9,C38TH=-.375_FLOAT)
      REAL(FLOAT),DIMENSION(:,:),ALLOCATABLE :: DC6IJ !dC6(iat,jat)/dCN(iat) in dc6ij(i,j)
      REAL(FLOAT),DIMENSION(:),ALLOCATABLE :: CN,C6SAVE,DC6I
      DIMENSION R0AB(MAX_ELEM,MAX_ELEM),R2R4(MAX_ELEM),RCOV(MAX_ELEM),
     *G(3,N)
      DIMENSION C6AB(MAX_ELEM,MAX_ELEM,MAXC,MAXC,3),MXC(MAX_ELEM),IZ(N)
      logical noabc
      REAL(FLOAT),DIMENSION(3,3) :: lat,stress,sigma,lat_1
      REAL(FLOAT),DIMENSION(3) :: tau,vec,dxyz,rij,rik,rjk,TUX,TUY,JTAU
      REAL(FLOAT),DIMENSION(3) :: TUX0,TUY0,TUZ0,IJVEC,IKVEC,JKVEC,KTAU
      REAL(FLOAT),DIMENSION(3) :: KTX,KTY
      integer,DIMENSION(3) :: rep_v,rep_cn
      integer :: a,b !czw
      AUTOEV=PAR(4)
      AUTOANG=PAR(32)
      sigma=0._FLOAT
      g(1:3,1:n)=0._FLOAT
      DISP=0._FLOAT
      IREPV1=REP_V(1)
      IREPV2=REP_V(2)
      IREPV3=REP_V(3)
      REPV1=REAL(-IREPV1,FLOAT)
      REPV2=REAL(-IREPV2,FLOAT)
      REPV3=REAL(-(IREPV3+1),FLOAT)
      if(nversion.eq.2)then
      do iat=1,n-1
       IZIAT=IZ(IAT)
         do jat=iat+1,n
           R0=1._FLOAT/(r0ab(iz(jat),IZIAT)*rs6)
           c6=c6ab(iz(jat),IZIAT,1,1,1)*s6
           do Jtaux=-IREPV1,IREPV1
           do Jtauy=-IREPV2,IREPV2
           do Jtauz=-IREPV3,IREPV3
            tau=REAL(Jtaux,FLOAT)*lat(:,1)+REAL(Jtauy,FLOAT)*lat(:,2)+
     *REAL(Jtauz,FLOAT)*lat(:,3)
              dxyz=XA(:,iat)-XA(:,jat)+tau            
            r2  =sum(dxyz*dxyz)
           if(r2.gt.CRIT_VDW) cycle
            r235=r2**(-3.5_FLOAT)
            r=sqrt(r2)
            damp6=exp(-alp6*(r*R0-1._FLOAT))
            R=1._FLOAT/R
            damp1=1._FLOAT/(damp6+1._FLOAT)
            R235=R235*DAMP1
            tmp1=damp6*R0*R235*DAMP1
            tmp2=6._FLOAT*r*R235
      term=(alp6*damp6*R0*damp1-6._FLOAT*R)*R235*C6
              g(:,iat)=g(:,iat)-term*dxyz
              g(:,jat)=g(:,jat)+term*dxyz
            disp=disp+c6*damp1*r2**(-3)
            do ny=1,3
              R20=DXYZ(NY)*TERM
            do my=1,3
              sigma(my,ny)=sigma(my,ny)+dxyz(my)*R20
            enddo !my
            enddo !ny
           enddo !tauz
           enddo !tauy
           enddo !taux
         enddo !jat
      enddo !iat
      do iat=1,n
       IZIAT=IZ(IAT)
           R0=1._FLOAT/(r0ab(IZIAT,IZIAT)*rs6)
           c6=c6ab(IZIAT,IZIAT,1,1,1)*s6
           C6H=C6*HALF
           do Jtaux=-IREPV1,IREPV1
           do Jtauy=-IREPV2,IREPV2
           do Jtauz=-IREPV3,IREPV3
            if(ABS(Jtaux)+ABS(Jtauy)+ABS(Jtauz).eq.0)cycle
            tau=REAL(Jtaux,FLOAT)*lat(:,1)+REAL(Jtauy,FLOAT)*lat(:,2)+
     *REAL(Jtauz,FLOAT)*lat(:,3)
            dxyz=tau
            r2=sum(dxyz*dxyz)
            if(r2.gt.CRIT_VDW) cycle
            r235=r2**(-3.5_FLOAT)
            r=sqrt(r2)
            damp6=exp(-alp6*(r*R0-1._FLOAT))
            R=1._FLOAT/R
            damp1=1._FLOAT/(damp6+1._FLOAT)
            R235=R235*DAMP1
            tmp1=damp6*R0*damp1*r235
            tmp2=6._FLOAT*r*r235
            disp=disp+c6*damp1*r2**(-3)*HALF
            term=(alp6*tmp1-tmp2)*C6H
            do ny=1,3
            TMP1=DXYZ(NY)*TERM
            do my=1,3
           sigma(my,ny)=sigma(my,ny)+term*dxyz(my)*TMP1
            enddo !my
            enddo !ny
           enddo !tauz
           enddo !tauy
           enddo !taux
      enddo !iat
      disp=-disp
      ELSE
      IF(NVERSION.GT.2)THEN
      IAT0=((N+1)*N)/2
      CALL CRYALLOC(DC6IJ,N,N,'PBCGDISP','DC6IJ')
      CALL CRYALLOC(C6SAVE,IAT0,'PBCGDISP','C6SAVE')
      CALL CRYALLOC(CN,N,'PBCGDISP','CN')
      CALL CRYALLOC(DC6I,N,'PBCGDISP','DC6I')
      DC6I(1:N)=0._FLOAT
      DC6IJ(1:N,1:N)=0._FLOAT
      C6SAVE(1:IAT0)=0._FLOAT
      DO I=1,3
      TUX0(I)=LAT(I,1)*REPV1
      TUY0(I)=LAT(I,2)*REPV2
      TUZ0(I)=LAT(I,3)*REPV3
      ENDDO
      IAT0=0
      INDI=0
      CALL IGPVAL
      call pbcncoord(n,rcov,iz,cn,lat,rep_cn,crit_cn)
      if(nversion.eq.3)then
      do iat=1,n
      INDI=INDI+IAT0
      IAT0=IAT0+1
      ITASK=ITASK+1
      IF(MPOINT.NE.ITASK)CYCLE
      IZIAT=IZ(IAT)
      call get_dC6_dCNij(maxc,max_elem,c6ab,mxc(IZIAT),mxc(IZIAT),cn(iat
     *),cn(iat),IZIAT,IZIAT,iat,iat,c6,dc6ij(iat,iat),fdum)
      C6SAVE(INDI+IAT)=C6
      DC6IIJ=DC6IJ(IAT,IAT)*2._FLOAT
        r0=r0ab(IZIAT,IZIAT)
        RR0=1._FLOAT/(RS6*R0)
        R20=1._FLOAT/(RS8*R0)
        R2R4I=R2R4(IZIAT)
        r42=R2R4I*R2R4I
        RCOVI=RCOV(IZIAT)
        rcovij=RCOVI+RCOVI
        TUX(1:3)=TUX0(1:3)
        do Jtaux=-IREPV1,IREPV1
        TUY(1:3)=TUY0(1:3)+TUX(1:3)
        do Jtauy=-IREPV2,IREPV2
        TAU(1:3)=TUZ0(1:3)+TUY(1:3)
        do Jtauz=-IREPV3,IREPV3
          TAU(1:3)=LAT(1:3,3)+TAU(1:3)
          R2=TAU(1)**2+TAU(2)**2+TAU(3)**2
          
          if(R2.GE.CRIT_VDW.OR.R2.LE.0.1_FLOAT)CYCLE
          r=sqrt(r2)
          r6=r2*r2*r2
          r8=r6*r2
          r9=r8*r2
          R8=1._FLOAT/R8
          R9=1._FLOAT/R9
          t6=(r*RR0)**(-alp6)
          damp6 =1._FLOAT/(t6*6._FLOAT+1._FLOAT)
          t8=(r*R20)**(-alp8)
          damp8=1._FLOAT/(t8*6._FLOAT+1._FLOAT)
      VEC(1:3)=((alp6*t6*damp6-1._FLOAT)*s6*R8*damp6*3._FLOAT+(alp8*t8*
     *damp8*9._FLOAT-12._FLOAT)*S18*R9*R42*DAMP8)*C6*TAU(1:3)
      DO I=1,3
      SIGMA(1:3,I)=VEC(1:3)*TAU(I)+SIGMA(1:3,I)
      enddo
      DC6_REST=(S18*r42*R8*damp8*3._FLOAT+s6/r6*damp6)*HALF
      disp=disp-DC6_REST*c6     ! calculate E_disp for sanity check
      DC6I(IAT)=DC6_REST*DC6IIJ+DC6I(IAT)
        ENDDO !tauz
        TUY(1:3)=LAT(1:3,2)+TUY(1:3)
        ENDDO !tauy
        TUX(1:3)=LAT(1:3,1)+TUX(1:3)
        ENDDO !taux
        RIJ(1:3)=-XA(1:3,IAT)
        do jat=1,iat-1
        IZJAT=IZ(JAT)
      call get_dC6_dCNij(maxc,max_elem,c6ab,mxc(IZIAT),mxc(IZJAT),cn(iat
     *),cn(jat),IZIAT,IZJAT,iat,jat,c6,dc6ij(iat,jat),dc6ij(jat,iat))
      C6SAVE(INDI+JAT)=C6
      DC6IIJ=DC6IJ(IAT,JAT)
      DC6IJI=DC6IJ(JAT,IAT)
      RIK(1:3)=XA(1:3,JAT)+RIJ(1:3)
          r0=r0ab(IZJAT,IZIAT)
        RR0=1._FLOAT/(RS6*R0)
        R20=1._FLOAT/(RS8*R0)
          r42=r2r4(IZJAT)*R2R4I
          rcovij=rcov(IZJAT)+RCOVI
        TUX(1:3)=TUX0(1:3)+RIK(1:3)
            do Jtaux=-IREPV1,IREPV1
        TUY(1:3)=TUY0(1:3)+TUX(1:3)
            do Jtauy=-IREPV2,IREPV2
        TAU(1:3)=TUZ0(1:3)+TUY(1:3)
            do Jtauz=-IREPV3,IREPV3
          TAU(1:3)=LAT(1:3,3)+TAU(1:3)
          R2=TAU(1)**2+TAU(2)**2+TAU(3)**2
            if(r2.gt.CRIT_VDW)CYCLE
            r=sqrt(r2)
            r6=r2*r2*r2
            r8=r6*r2
            r9=r8*r2
          R8=1._FLOAT/R8
          R9=1._FLOAT/R9
            t6=(r*RR0)**(-alp6)
            damp6=1._FLOAT/(1._FLOAT+6._FLOAT*t6)
            t8=(r*R20)**(-alp8)
            damp8=1._FLOAT/(1._FLOAT+6._FLOAT*t8)
      VEC(1:3)=((alp6*t6*damp6-1._FLOAT)*s6*R8*damp6*6._FLOAT+(alp8*t8*
     *damp8*18._FLOAT-24._FLOAT)*S18*R9*R42*DAMP8)*C6*TAU(1:3)
      DO I=1,3
      SIGMA(1:3,I)=VEC(1:3)*TAU(I)+SIGMA(1:3,I)
      G(I,IAT)=G(I,IAT)+VEC(I)
      G(I,JAT)=G(I,JAT)-VEC(I)
      enddo
      DC6_REST=S18*r42*R8*damp8*3._FLOAT+s6/r6*damp6
      disp=disp-DC6_REST*c6  ! calculate E_disp for sanity check
      DC6I(IAT)=DC6_REST*DC6IIJ+DC6I(IAT)
      DC6I(JAT)=DC6_REST*DC6IJI+DC6I(JAT)
          enddo !tauz
          TUY(1:3)=LAT(1:3,2)+TUY(1:3)
          enddo !tauy
        TUX(1:3)=LAT(1:3,1)+TUX(1:3)
          enddo !taux
        enddo !jat
      CALL IGPVAL
      enddo !iat
      elseif(nversion.eq.4)then
      do iat=1,n
      INDI=INDI+IAT0
      IAT0=IAT0+1
      ITASK=ITASK+1
      IF(MPOINT.NE.ITASK)CYCLE
      IZIAT=IZ(IAT)
      call get_dC6_dCNij(maxc,max_elem,c6ab,mxc(IZIAT),mxc(IZIAT),
     *cn(iat),cn(iat),IZIAT,IZIAT,iat,iat,c6,dc6ij(iat,iat),fdum)
      C6SAVE(INDI+IAT)=C6
      R2R4I=R2R4(IZIAT)
      r42=R2R4I*R2R4I
        RCOVI=RCOV(IZIAT)
        rcovij=RCOVI+RCOVI
        FDUM=R42*3._FLOAT
        R0=SQRT(FDUM)*RS6+RS8
        FDUO=FDUM*S18
        FDUM=-FDUO*4._FLOAT
        FDUN=-S6*3._FLOAT
        R20=R0**6
        R0=R20*R0*R0
        TUX(1:3)=TUX0(1:3)
        DC6IAT=0._FLOAT
        do Jtaux=-IREPV1,IREPV1
        TUY(1:3)=TUY0(1:3)+TUX(1:3)
        do Jtauy=-IREPV2,IREPV2
        TAU(1:3)=TUZ0(1:3)+TUY(1:3)
        do Jtauz=-IREPV3,IREPV3
          TAU(1:3)=LAT(1:3,3)+TAU(1:3)
          R2=SUM(TAU(1:3)**2)
          if(R2.GE.CRIT_VDW.OR.R2.LE.0.1_FLOAT)CYCLE
          r4=r2*r2
          r6=r4*r2
          r8=r6*r2
          t6=1._FLOAT/(r6+R20)
          t8=1._FLOAT/(r8+R0)
      vec(1:3)=(FDUN*r4*t6*t6+FDUM*r6*t8*t8)*C6*tau(1:3)
          DO I=1,3
          SIGMA(1:3,I)=VEC(1:3)*TAU(I)+SIGMA(1:3,I)
          enddo !i
          DC6_REST=(s6*t6+FDUO*t8)*HALF
          disp=disp-DC6_REST*c6  ! calculate E_disp for sanity check
        DC6IAT=DC6IAT+DC6_REST
        ENDDO !tauz
        TUY(1:3)=LAT(1:3,2)+TUY(1:3)
        ENDDO !tauy
        TUX(1:3)=LAT(1:3,1)+TUX(1:3)
        ENDDO !taux
      DC6I(IAT)=DC6IJ(IAT,IAT)*DC6IAT*2._FLOAT+DC6I(IAT)
      DC6IAT=0._FLOAT
        RIJ(1:3)=-XA(1:3,IAT)
        do jat=1,iat-1
        IZJAT=IZ(JAT)
      call get_dC6_dCNij(maxc,max_elem,c6ab,mxc(IZIAT),mxc(IZJAT),cn(iat
     *),cn(jat),IZIAT,IZJAT,iat,jat,c6,dc6ij(iat,jat),dc6ij(jat,iat))
      C6SAVE(INDI+JAT)=C6
      RIK(1:3)=XA(1:3,JAT)+RIJ(1:3)
        r42=r2r4(IZJAT)*R2R4I
        rcovij=rcov(IZJAT)+RCOVI
        FDUM=R42*3._FLOAT
        R0=SQRT(FDUM)*RS6+RS8
        FDUO=FDUM*S18
        FDUM=-FDUO*8._FLOAT
        FDUN=-S6*6._FLOAT
        R20=R0**6
        R0=R20*R0*R0
      DC6JAT=0._FLOAT
        TUX(1:3)=TUX0(1:3)+RIK(1:3)
            do Jtaux=-IREPV1,IREPV1
        TUY(1:3)=TUY0(1:3)+TUX(1:3)
            do Jtauy=-IREPV2,IREPV2
        TAU(1:3)=TUZ0(1:3)+TUY(1:3)
            do Jtauz=-IREPV3,IREPV3
          TAU(1:3)=LAT(1:3,3)+TAU(1:3)
          R2=SUM(TAU(1:3)**2)
            if(r2.gt.CRIT_VDW)CYCLE
            r4=r2*r2
            r6=r4*r2
            r8=r6*r2
            t6=1._FLOAT/(r6+R20)
            t8=1._FLOAT/(r8+R0)
      vec(1:3)=(FDUN*r4*t6*t6+FDUM*r6*t8*t8)*C6*tau(1:3)
      do i=1,3
      SIGMA(1:3,I)=VEC(1:3)*TAU(I)+SIGMA(1:3,I)
      G(I,IAT)=G(I,IAT)+VEC(I)
      G(I,JAT)=G(I,JAT)-VEC(I)
      enddo !i
      DC6_REST=s6*t6+FDUO*t8
      disp=disp-DC6_REST*c6
      DC6JAT=DC6_REST+DC6JAT
          enddo !tauz
          TUY(1:3)=LAT(1:3,2)+TUY(1:3)
          enddo !tauy
        TUX(1:3)=LAT(1:3,1)+TUX(1:3)
          enddo !taux
      DC6IAT=DC6IJ(IAT,JAT)*DC6JAT+DC6IAT
      DC6I(JAT)=DC6IJ(JAT,IAT)*DC6JAT+DC6I(JAT)
        enddo !jat
      CALL IGPVAL
      DC6I(IAT)=DC6I(IAT)+DC6IAT
      enddo !iat
      ENDIF
      CALL GSUM(DC6IJ,N*N)
      CALL GSUM(C6SAVE,((N+1)*N)/2)
      IREPV1=REP_CN(1)
      IREPV2=REP_CN(2)
      IREPV3=REP_CN(3)
      MREPV1=-IREPV1
      MREPV2=-IREPV2
      MREPV3=-IREPV3
      REPV1=REAL(MREPV1,FLOAT)
      REPV2=REAL(MREPV2,FLOAT)
      REPV3=REAL(MREPV3-1,FLOAT)
      if(.not.noabc)then
      eabc=0._FLOAT
      DO I=1,3
      TUX0(I)=LAT(I,1)*REPV1
      TUY0(I)=LAT(I,2)*REPV2
      TUZ0(I)=LAT(I,3)*REPV3
      ENDDO
      IAT0=2
      INDI=1
      do iat=3,n
      IZIAT=IZ(IAT)
      DC6IAT=0._FLOAT
      INDI=INDI+IAT0
      IAT0=IAT0+1
      INDII=0
      JAT0=1
      do jat=2,iat-1
      INDII=INDII+JAT0
      JAT0=JAT0+1
      ITASK=ITASK+1
      IF(MPOINT.NE.ITASK)CYCLE
      DC6JAT=0._FLOAT
      c6ijD=c6save(INDI+JAT)
      C6IJ=1._FLOAT/C6IJD
      ijvec(1:3)=XA(1:3,jat)-XA(1:3,iat)
      IZJAT=IZ(JAT)
      R0IJ=1._FLOAT/R0AB(IZIAT,IZJAT)
      do kat=1,jat-1
      ikvec(1:3)=XA(1:3,kat)-XA(1:3,iat)
      jkvec(1:3)=XA(1:3,kat)-XA(1:3,jat)
      DC6KAT=0._FLOAT
      c6ik=c6save(INDI+KAT)
      c6jk=c6save(INDII+KAT)
      c9=-SQRT(c6ijD*c6ik*c6jk)
      C96=C9*HALF
      FDUM=0._FLOAT
      C6IK=1._FLOAT/C6IK
      C6JK=1._FLOAT/C6JK
      IZKAT=IZ(KAT)
      R0IK=R0IJ/(R0AB(IZIAT,IZKAT)*R0AB(IZJAT,IZKAT))
      TUX(1:3)=TUX0(1:3)
      do jtaux=MREPV1,IREPV1
      IREPMIN1=MAX(MREPV1,JTAUX+MREPV1)
      IREPMAX1=MIN(IREPV1,JTAUX+IREPV1)
      TUY(1:3)=TUY0(1:3)+TUX(1:3)
      do jtauy=MREPV2,IREPV2
      IREPMIN2=MAX(MREPV2,JTAUY+MREPV2)
      IREPMAX2=MIN(IREPV2,JTAUY+IREPV2)
      JTAU(1:3)=TUZ0(1:3)+TUY(1:3)
      do jtauz=MREPV3,IREPV3
      IREPMIN3=MAX(MREPV3,JTAUZ+MREPV3)
      IREPMAX3=MIN(IREPV3,JTAUZ+IREPV3)
      JTAU(1:3)=LAT(1:3,3)+JTAU(1:3)
      RIJ(1:3)=IJVEC(1:3)+JTAU(1:3)
      rij2=SUM(RIJ(1:3)**2)
      if(rij2.gt.CRIT_ABC)cycle
      RIJ2I=C9/RIJ2
      rr0ij=SQRT(rij2)*R0IK
      RIJ2Q=RIJ2*RIJ2
      RIJ2C=RIJ2Q*RIJ2
      KTX(1:3)=LAT(1:3,1)*REAL(IREPMIN1,FLOAT)
      DXYZ(1:3)=LAT(1:3,2)*REAL(IREPMIN2,FLOAT)
      TAU(1:3)=LAT(1:3,3)*REAL(IREPMIN3-1,FLOAT)
      do ktaux=IREPMIN1,IREPMAX1
      KTY(1:3)=DXYZ(1:3)+KTX(1:3)
      do ktauy=IREPMIN2,IREPMAX2
      KTAU(1:3)=TAU(1:3)+KTY(1:3)
      do ktauz=IREPMIN3,IREPMAX3
      KTAU(1:3)=LAT(1:3,3)+KTAU(1:3)
      RIK(1:3)=KTAU(1:3)+IKVEC(1:3)
      rik2=SUM(RIK(1:3)**2)
      if(rik2.gt.CRIT_ABC)cycle
      RJK(1:3)=KTAU(1:3)-JTAU(1:3)+JKVEC(1:3)
      RJK2=SUM(RJK(1:3)**2)
      if(rjk2.gt.CRIT_ABC)cycle
      RIK2Q=RIK2*RIK2
      RJK2Q=RJK2*RJK2
      geomean2=1._FLOAT/(rij2*rjk2*rik2)
      r0av=(SQRT(rik2*RJK2)*rr0ij)**THIRD
      damp9=1._FLOAT/(6._FLOAT*(sr9*r0av)**alp9+1._FLOAT)  !alp9 is already saved with "-"
      geomean3=SQRT(GEOMEAN2)*geomean2
      geomean5=GEOMEAN3*geomean2
      ang=(rij2+rjk2-rik2)*(rij2-rjk2+rik2)*(-rij2+rjk2+rik2)*GEOMEAN5*
     *.375_FLOAT+GEOMEAN3
      dc6_rest=ang*damp9
      eabc=dc6_rest*c9+eabc
      dfdmp=(0.75_FLOAT*r0av)**alp9*damp9*damp9*ALP92
      dang=(RIJ2C+RIJ2Q*(rjk2+rik2)+((RJK2Q+RIK2Q)*3._FLOAT+(RJK2+RJK2)*
     *rik2)*RIJ2+(RIK2-RJK2)*(RJK2Q-RIK2Q)*5._FLOAT)*geomean5*C38TH
      VEC(1:3)=(dang*damp9-dfdmp*ang)*RIJ2I*RIJ(1:3)
      DO I=1,3
      SIGMA(1:3,I)=RIJ(I)*VEC(1:3)+SIGMA(1:3,I)
      G(I,IAT)=G(I,IAT)+VEC(I)
      G(I,JAT)=G(I,JAT)-VEC(I)
      ENDDO
      dang=(RIK2Q*RIK2+RIK2Q*(rjk2+rij2)+((RJK2Q+RIJ2Q)*3._FLOAT+(RJK2+R
     *JK2)*rij2)*RIK2+(RIJ2-RJK2)*(RJK2Q-RIJ2Q)*5._FLOAT)*geomean5*C38TH
      VEC(1:3)=(dang*damp9-dfdmp*ang)/RIK2*C9*RIK(1:3)
      DO I=1,3
      SIGMA(1:3,I)=RIK(I)*VEC(1:3)+SIGMA(1:3,I)
      G(I,IAT)=G(I,IAT)+VEC(I)
      G(I,KAT)=G(I,KAT)-VEC(I)
      ENDDO
      dang=(RJK2Q*RJK2+RJK2Q*(rik2+rij2)+((RIK2Q+RIJ2Q)*3._FLOAT+(RIK2+R
     *IK2)*rij2)*RJK2+(RIJ2-RIK2)*(RIK2Q-RIJ2Q)*5._FLOAT)*geomean5*C38TH
      VEC(1:3)=(dang*damp9-dfdmp*ang)/RJK2*C9*RJK(1:3)
      DO I=1,3
      SIGMA(1:3,I)=RJK(I)*VEC(1:3)+SIGMA(1:3,I)
      G(I,JAT)=G(I,JAT)+VEC(I)
      G(I,KAT)=G(I,KAT)-VEC(I)
      ENDDO
      FDUM=FDUM+dc6_rest
      enddo !ktauz
      KTY(1:3)=LAT(1:3,2)+KTY(1:3)
      enddo !ktauy
      KTX(1:3)=LAT(1:3,1)+KTX(1:3)
      enddo !ktaux
      enddo !jtauz
      TUY(1:3)=LAT(1:3,2)+TUY(1:3)
      enddo !jtauy
      TUX(1:3)=LAT(1:3,1)+TUX(1:3)
      enddo !jtaux
      FDUM=FDUM*C96
      DC6IAT=(dc6ij(iat,jat)*C6IJ+dc6ij(iat,kat)*C6IK)*FDUM+DC6IAT
      DC6JAT=(dc6ij(jat,iat)*C6IJ+dc6ij(jat,kat)*C6JK)*FDUM+DC6JAT
      DC6I(KAT)=(dc6ij(kat,iat)*C6IK+dc6ij(kat,jat)*C6JK)*FDUM+DC6I(KAT)
      enddo !kat
      CALL IGPVAL
      DC6I(JAT)=DC6I(JAT)+DC6JAT
      enddo !jat
      DC6I(IAT)=DC6I(IAT)+DC6IAT
      enddo !iat
      IAT0=1
      INDI=0
      DO I=1,3
      TUX0(I)=LAT(I,1)*REPV1
      TUY0(I)=LAT(I,2)*REPV2
      TUZ0(I)=LAT(I,3)*REPV3
      ENDDO
      DO iat=2,n
      INDI=INDI+IAT0
      IAT0=IAT0+1
      DC6IAT=0._FLOAT
      c6ijD=c6save(INDI+IAT)
      C6IJ=1._FLOAT/C6IJD
      DC6IJI=dc6ij(iat,IAT)
      IZIAT=IZ(IAT)
      R0IJ=1._FLOAT/r0ab(IZIAT,IZIAT)
      DO kat=1,iat-1
      ITASK=ITASK+1
      IF(MPOINT.NE.ITASK)CYCLE
      c6ik=c6save(INDI+KAT)
      IZKAT=IZ(KAT)
      R0IK=R0IJ/(r0ab(IZIAT,IZKAT)**2)
      ikvec(1:3)=XA(1:3,kat)-XA(1:3,iat)
      c9=-sqrt(c6ijD*c6ik*c6Ik)
      C96=C9*HALF
      C6IK=1._FLOAT/C6IK
      FDUM=0._FLOAT
      TUX(1:3)=TUX0(1:3)
      do jtaux=MREPV1,IREPV1
      IREPMIN1=MAX(MREPV1,JTAUX+MREPV1)
      IREPMAX1=MIN(IREPV1,JTAUX+IREPV1)
      TUY(1:3)=TUY0(1:3)+TUX(1:3)
      do jtauy=MREPV2,IREPV2
      IREPMIN2=MAX(MREPV2,JTAUY+MREPV2)
      IREPMAX2=MIN(IREPV2,JTAUY+IREPV2)
      JTAU(1:3)=TUZ0(1:3)+TUY(1:3)
      do jtauz=MREPV3,IREPV3
      IREPMIN3=MAX(MREPV3,JTAUZ+MREPV3)
      IREPMAX3=MIN(IREPV3,JTAUZ+IREPV3)
      JTAU(1:3)=LAT(1:3,3)+JTAU(1:3)
      IF (ABS(jtaux)+ABS(jtauy)+ABS(jtauz).eq.0)cycle
      rij2=SUM(JTAU(1:3)**2)
      if(rij2.gt.CRIT_ABC)cycle
      RIJ2I=C96/RIJ2
      rr0ij=SQRT(rij2)*R0IK
      RIJ2Q=RIJ2*RIJ2
      RIJ2C=RIJ2Q*RIJ2
      KTX(1:3)=LAT(1:3,1)*REAL(IREPMIN1,FLOAT)
      DXYZ(1:3)=LAT(1:3,2)*REAL(IREPMIN2,FLOAT)
      TAU(1:3)=LAT(1:3,3)*REAL(IREPMIN3-1,FLOAT)
      do ktaux=IREPMIN1,IREPMAX1
      KTY(1:3)=DXYZ(1:3)+KTX(1:3)
      do ktauy=IREPMIN2,IREPMAX2
      KTAU(1:3)=TAU(1:3)+KTY(1:3)
      do ktauz=IREPMIN3,IREPMAX3
      KTAU(1:3)=LAT(1:3,3)+KTAU(1:3)
      RIK(1:3)=KTAU(1:3)+IKVEC(1:3)
      rik2=SUM(RIK(1:3)**2)
      if(rik2.gt.CRIT_ABC)cycle
      RJK(1:3)=KTAU(1:3)-JTAU(1:3)+IKVEC(1:3)
      RJK2=SUM(RJK(1:3)**2)
      if(rjk2.gt.CRIT_ABC)cycle
      RIK2Q=RIK2*RIK2
      RJK2Q=RJK2*RJK2
      geomean2=1._FLOAT/(rij2*rjk2*rik2)
      r0av=(SQRT(rik2*RJK2)*rr0ij)**THIRD
      damp9=1._FLOAT/(1._FLOAT+6._FLOAT*(sr9*r0av)**alp9)  !alp9 is already saved with "-"
      GEOMEAN3=SQRT(GEOMEAN2)*GEOMEAN2
      GEOMEAN5=GEOMEAN3*GEOMEAN2
      ang=(rij2+rjk2-rik2)*(rij2-rjk2+rik2)*(-rij2+rjk2+rik2)*GEOMEAN5*
     *.375_FLOAT+GEOMEAN3
      dc6_rest=ang*damp9*HALF   !factor 1/2 for doublecounting
      eabc=dc6_rest*c9+eabc
      dfdmp=(0.75_FLOAT*r0av)**(alp9)*damp9*damp9*alp92
      dang=(RIJ2C+RIJ2Q*(rjk2+rik2)+((RJK2Q+RIK2Q)*3._FLOAT+(RJK2+RJK2)*
     *rik2)*RIJ2+(RIK2-RJK2)*(RJK2Q-RIK2Q)*5._FLOAT)*GEOMEAN5*C38TH
      VEC(1:3)=(dang*damp9-dfdmp*ang)*RIJ2I*JTAU(1:3)
      DO I=1,3
      SIGMA(1:3,I)=JTAU(I)*VEC(1:3)+SIGMA(1:3,I)
      ENDDO
      dang=(RIK2Q*RIK2+RIK2Q*(rjk2+rij2)+((RJK2Q+RIJ2Q)*3._FLOAT+(RJK2+R
     *JK2)*rij2)*RIK2+(RIJ2-RJK2)*(RJK2Q-RIJ2Q)*5._FLOAT)*GEOMEAN5*C38TH
      VEC(1:3)=(dang*damp9-dfdmp*ang)/RIK2*C96*RIK(1:3)
      DO I=1,3
      SIGMA(1:3,I)=RIK(I)*VEC(1:3)+SIGMA(1:3,I)
      G(I,IAT)=G(I,IAT)+VEC(I)
      G(I,KAT)=G(I,KAT)-VEC(I)
      ENDDO
      dang=(RJK2Q*RJK2+RJK2Q*(rik2+rij2)+((RIK2Q+RIJ2Q)*3._FLOAT+(RIK2+R
     *IK2)*rij2)*RJK2+(RIJ2-RIK2)*(RIK2Q-RIJ2Q)*5._FLOAT)*GEOMEAN5*C38TH
      VEC(1:3)=(dang*damp9-dfdmp*ang)/RJK2*C96*RJK(1:3)
      DO I=1,3
      SIGMA(1:3,I)=RJK(I)*VEC(1:3)+SIGMA(1:3,I)
      G(I,IAT)=G(I,IAT)+VEC(I)
      G(I,KAT)=G(I,KAT)-VEC(I)
      ENDDO
      FDUM=FDUM+dc6_rest
      ENDDO !ktauz
      KTY(1:3)=LAT(1:3,2)+KTY(1:3)
      ENDDO !ktauy
      KTX(1:3)=LAT(1:3,1)+KTX(1:3)
      ENDDO !ktaux
      ENDDO !jtauz
      TUY(1:3)=LAT(1:3,2)+TUY(1:3)
      ENDDO !jtauy
      TUX(1:3)=LAT(1:3,1)+TUX(1:3)
      ENDDO !jtaux
      FDUM=FDUM*C9
      DC6IAT=(DC6IJI*c6ij+dc6ij(iat,kat)*c6ik)*FDUM+DC6IAT
      DC6I(KAT)=dc6ij(kat,iat)*c6ik*FDUM+DC6I(KAT)
      CALL IGPVAL
      ENDDO !kat
      DC6I(IAT)=DC6I(IAT)+DC6IAT
      ENDDO !iat
      INDI=0
      IAT0=1
      DO I=1,3
      TUX0(I)=LAT(I,1)*REPV1
      TUY0(I)=LAT(I,2)*REPV2
      TUZ0(I)=LAT(I,3)*REPV3
      ENDDO
      DO iat=2,n
      IZIAT=IZ(IAT)
      DC6IAT=0._FLOAT
      INDI=INDI+IAT0
      IAT0=IAT0+1
      INDII=0
      JAT0=0
      DO jat=1,iat-1
      ITASK=ITASK+1
      INDII=INDII+JAT0
      JAT0=JAT0+1
      IF(MPOINT.NE.ITASK)CYCLE
      c6ij=c6save(INDI+JAT)
      c6jk=c6save(INDII+JAT)
      IJVEC(1:3)=XA(1:3,JAT)-XA(1:3,IAT)
      c9=-SQRT(C6IJ*C6IJ*C6JK)
      C96=C9*HALF
      C6IJ=1._FLOAT/C6IJ
      C6JK=1._FLOAT/C6JK
      IZJAT=IZ(JAT)
      R0IJ=1._FLOAT/(R0AB(IZIAT,IZJAT)**2*R0AB(IZJAT,IZJAT))
      FDUM=0._FLOAT
      TUX(1:3)=TUX0(1:3)
      do jtaux=MREPV1,IREPV1
      IREPMIN1=MAX(MREPV1,JTAUX+MREPV1)
      IREPMAX1=MIN(IREPV1,JTAUX+IREPV1)
      TUY(1:3)=TUY0(1:3)+TUX(1:3)
      do jtauy=MREPV2,IREPV2
      IREPMIN2=MAX(MREPV2,JTAUY+MREPV2)
      IREPMAX2=MIN(IREPV2,JTAUY+IREPV2)
      JTAU(1:3)=TUZ0(1:3)+TUY(1:3)
      do jtauz=MREPV3,IREPV3
      IREPMIN3=MAX(MREPV3,JTAUZ+MREPV3)
      IREPMAX3=MIN(IREPV3,JTAUZ+IREPV3)
      JTAU(1:3)=LAT(1:3,3)+JTAU(1:3)
      RIJ(1:3)=IJVEC(1:3)+JTAU(1:3)
      rij2=SUM(RIJ(1:3)**2)
      if(rij2.gt.CRIT_ABC)cycle
      RIJ2I=C96/RIJ2
      rr0ij=SQRT(rij2)*R0IJ
      RIJ2Q=RIJ2*RIJ2
      RIJ2C=RIJ2Q*RIJ2
      KTX(1:3)=LAT(1:3,1)*REAL(IREPMIN1,FLOAT)
      DXYZ(1:3)=LAT(1:3,2)*REAL(IREPMIN2,FLOAT)
      TAU(1:3)=LAT(1:3,3)*REAL(IREPMIN3-1,FLOAT)
      do ktaux=IREPMIN1,IREPMAX1
      KTY(1:3)=DXYZ(1:3)+KTX(1:3)
      do ktauy=IREPMIN2,IREPMAX2
      KTAU(1:3)=TAU(1:3)+KTY(1:3)
      do ktauz=IREPMIN3,IREPMAX3
      KTAU(1:3)=LAT(1:3,3)+KTAU(1:3)
      IF(ABS(jtaux-ktaux)+ABS(jtauy-ktauy)+ABS(jtauz-ktauz).EQ.0)cycle
      RIK(1:3)=KTAU(1:3)+IJVEC(1:3)
      rik2=SUM(RIK(1:3)**2)
      if(rik2.gt.CRIT_ABC)cycle
      RJK(1:3)=KTAU(1:3)-JTAU(1:3)
      rjk2=SUM(RJK(1:3)**2)
      if(rjk2.gt.CRIT_ABC)cycle
      RIK2Q=RIK2*RIK2
      RJK2Q=RJK2*RJK2
      geomean2=1._FLOAT/(rij2*rjk2*rik2)
      r0av=(SQRT(rik2*rjk2)*rr0ij)**THIRD
      damp9=1._FLOAT/(1._FLOAT+6._FLOAT*(sr9*r0av)**alp9)  !alp9 is already saved with "-"
      GEOMEAN3=SQRT(GEOMEAN2)*GEOMEAN2
      GEOMEAN5=GEOMEAN3*GEOMEAN2
      ang=(rij2+rjk2-rik2)*(rij2-rjk2+rik2)*(-rij2+rjk2+rik2)*GEOMEAN5*
     *0.375_FLOAT+GEOMEAN3
      dc6_rest=ang*damp9*HALF   !factor 1/2 for doublecounting
      eabc=eabc+dc6_rest*c9
      dfdmp=(R0AV*0.75_FLOAT)**(alp9)*damp9*damp9*ALP92
      dang=(RIJ2C+RIJ2Q*(rjk2+rik2)+((RJK2Q+RIK2Q)*3._FLOAT+(RJK2+RJK2)*
     *rik2)*RIJ2+(RIK2-RJK2)*(RJK2Q-RIK2Q)*5._FLOAT)*GEOMEAN5*C38TH
      VEC(1:3)=(dang*damp9-dfdmp*ang)*RIJ2I*RIJ(1:3)
      DO I=1,3
      SIGMA(1:3,I)=RIJ(I)*VEC(1:3)+SIGMA(1:3,I)
      G(I,IAT)=G(I,IAT)+VEC(I)
      G(I,JAT)=G(I,JAT)-VEC(I)
      ENDDO
      dang=(RIK2Q*RIK2+RIK2Q*(rjk2+rij2)+((RJK2Q+RIJ2Q)*3._FLOAT+(RJK2+R
     *JK2)*rij2)*RIK2+(RIJ2-RJK2)*(RJK2Q-RIJ2Q)*5._FLOAT)*GEOMEAN5*C38TH
      VEC(1:3)=(dang*damp9-dfdmp*ang)/RIK2*C96*RIK(1:3)
      DO I=1,3
      SIGMA(1:3,I)=RIK(I)*VEC(1:3)+SIGMA(1:3,I)
      G(I,IAT)=G(I,IAT)+VEC(I)
      G(I,JAT)=G(I,JAT)-VEC(I)
      ENDDO
      dang=(RJK2Q*RJK2+RJK2Q*(rik2+rij2)+((RIK2Q+RIJ2Q)*3._FLOAT+(RIK2+R
     *IK2)*rij2)*RJK2+(RIJ2-RIK2)*(RIK2Q-RIJ2Q)*5._FLOAT)*GEOMEAN5*C38TH
      VEC(1:3)=(dang*damp9-dfdmp*ang)/RJK2*C96*RJK(1:3)
      DO I=1,3
      SIGMA(1:3,I)=RJK(I)*VEC(1:3)+SIGMA(1:3,I)
      ENDDO
      FDUM=FDUM+dc6_rest
      ENDDO !ktauz
      KTY(1:3)=LAT(1:3,2)+KTY(1:3)
      ENDDO !ktauy
      KTX(1:3)=LAT(1:3,1)+KTX(1:3)
      ENDDO !ktaux
      ENDDO !jtauz
      TUY(1:3)=LAT(1:3,2)+TUY(1:3)
      ENDDO !jtauy
      TUX(1:3)=LAT(1:3,1)+TUX(1:3)
      ENDDO !jtaux
      CALL IGPVAL
      FDUM=FDUM*C9
      DC6IAT=dc6ij(iat,jat)*c6ij*FDUM+DC6IAT
      dc6i(jat)=(dc6ij(jat,iat)*c6ij+dc6ij(jat,JAT)*c6jk)*FDUM+dc6i(jat)
      ENDDO !JAT
      DC6I(IAT)=DC6I(IAT)+DC6IAT
      ENDDO !iat
      IAT0=0
      INDI=0
      DO I=1,3
      TUX0(I)=LAT(I,1)*REPV1
      TUY0(I)=LAT(I,2)*REPV2
      TUZ0(I)=LAT(I,3)*REPV3
      ENDDO
      DO iat=1,n
      INDI=INDI+IAT0+1
      IAT0=IAT0+1
      ITASK=ITASK+1
      IF(MPOINT.NE.ITASK)CYCLE
      c6ij=c6save(INDI)
      c9=-SQRT(c6ij**3)
      C96=C9*SIXTH
      C6IJ=1._FLOAT/C6IJ
      IZIAT=IZ(IAT)
      R0IJ=1._FLOAT/(R0AB(IZIAT,IZIAT)**3)
      FDUM=0._FLOAT
      TUX(1:3)=TUX0(1:3)
      do jtaux=MREPV1,IREPV1
      IREPMIN1=MAX(MREPV1,JTAUX+MREPV1)
      IREPMAX1=MIN(IREPV1,JTAUX+IREPV1)
      TUY(1:3)=TUY0(1:3)+TUX(1:3)
      do jtauy=MREPV2,IREPV2
      IREPMIN2=MAX(MREPV2,JTAUY+MREPV2)
      IREPMAX2=MIN(IREPV2,JTAUY+IREPV2)
      JTAU(1:3)=TUZ0(1:3)+TUY(1:3)
      do jtauz=MREPV3,IREPV3
      IREPMIN3=MAX(MREPV3,JTAUZ+MREPV3)
      IREPMAX3=MIN(IREPV3,JTAUZ+IREPV3)
      JTAU(1:3)=LAT(1:3,3)+JTAU(1:3)
      IF (ABS(jtaux)+ABS(jtauy)+ABS(jtauz).eq.0)cycle !IF iat and jat are the same then cycle
      rij2=SUM(JTAU(1:3)**2)
      if(rij2.gt.CRIT_ABC)cycle
      RIJ2I=C96/RIJ2
      rr0ij=SQRT(RIJ2)*R0IJ
      RIJ2Q=RIJ2*RIJ2
      RIJ2C=RIJ2Q*RIJ2
      KTX(1:3)=LAT(1:3,1)*REAL(IREPMIN1,FLOAT)
      DXYZ(1:3)=LAT(1:3,2)*REAL(IREPMIN2,FLOAT)
      TAU(1:3)=LAT(1:3,3)*REAL(IREPMIN3-1,FLOAT)
      do ktaux=IREPMIN1,IREPMAX1
      KTY(1:3)=DXYZ(1:3)+KTX(1:3)
      do ktauy=IREPMIN2,IREPMAX2
      KTAU(1:3)=TAU(1:3)+KTY(1:3)
      do ktauz=IREPMIN3,IREPMAX3
      KTAU(1:3)=LAT(1:3,3)+KTAU(1:3)
      IF(ABS(ktaux)+ABS(ktauy)+ABS(ktauz).eq.0.OR.ABS(ktaux-jtaux)+
     *ABS(ktauy-jtauy)+ABS(ktauz-jtauz).EQ.0)cycle      !If kat and jat are the same then cycle
      rik2=SUM(KTAU(1:3)**2)
      if(rik2.gt.CRIT_ABC)cycle
      RJK(1:3)=KTAU(1:3)-JTAU(1:3)
      rjk2=SUM(RJK(1:3)**2)
      if(rjk2.gt.CRIT_ABC)cycle
      RIK2Q=RIK2*RIK2
      RJK2Q=RJK2*RJK2
      GEOMEAN2=1._FLOAT/(rij2*rjk2*rik2)
      r0av=(SQRT(rik2*rjk2)*rr0ij)**THIRD
      damp9=1._FLOAT/(1._FLOAT+6._FLOAT*(sr9*r0av)**alp9)  !alp9 is already saved with "-"
      GEOMEAN3=SQRT(GEOMEAN2)*GEOMEAN2
      GEOMEAN5=GEOMEAN3*GEOMEAN2
      ang=(rij2+rjk2-rik2)*(rij2-rjk2+rik2)*(-rij2+rjk2+rik2)*GEOMEAN5*
     *0.375_FLOAT+GEOMEAN3
      dc6_rest=ang*damp9*SIXTH
      eabc=eabc+c9*dc6_rest
      dfdmp=(0.75_FLOAT*r0av)**(alp9)*damp9*damp9*ALP92
      dang=(RIJ2C+RIJ2Q*(rjk2+rik2)+((RJK2Q+RIK2Q)*3._FLOAT+(RJK2+RJK2)*
     *rik2)*RIJ2+(RIK2-RJK2)*(RJK2Q-RIK2Q)*5._FLOAT)*GEOMEAN5*C38TH
      VEC(1:3)=(dang*damp9-dfdmp*ang)*RIJ2I*JTAU(1:3)
      DO I=1,3
      SIGMA(1:3,I)=JTAU(I)*VEC(1:3)+SIGMA(1:3,I)
      ENDDO
      dang=(RIK2Q*RIK2+RIK2Q*(rjk2+rij2)+((RJK2Q+RIJ2Q)*3._FLOAT+(RIJ2+R
     *IJ2)*rjk2)*rik2+(RIJ2-RJK2)*(RJK2Q-RIJ2Q)*5._FLOAT)*GEOMEAN5*C38TH
      VEC(1:3)=(dang*damp9-dfdmp*ang)/RIK2*C96*KTAU(1:3)
      DO I=1,3
      SIGMA(1:3,I)=KTAU(I)*VEC(1:3)+SIGMA(1:3,I)
      ENDDO
      dang=(RJK2Q*RJK2+RJK2Q*(rik2+rij2)+((RIK2Q+RIJ2Q)*3._FLOAT+(RIK2+R
     *IK2)*rij2)*rjk2+(RIJ2-RIK2)*(RIK2Q-RIJ2Q)*5._FLOAT)*GEOMEAN5*C38TH
      VEC(1:3)=(dang*damp9-dfdmp*ang)/RJK2*C96*RJK(1:3)
      DO I=1,3
      SIGMA(1:3,I)=RJK(I)*VEC(1:3)+SIGMA(1:3,I)
      ENDDO
      FDUM=FDUM+dc6_rest
      ENDDO !ktauz
      KTY(1:3)=LAT(1:3,2)+KTY(1:3)
      ENDDO !ktauy
      KTX(1:3)=LAT(1:3,1)+KTX(1:3)
      ENDDO !ktaux
      ENDDO !jtauz
      TUY(1:3)=LAT(1:3,2)+TUY(1:3)
      ENDDO !jtauy
      TUX(1:3)=LAT(1:3,1)+TUX(1:3)
      ENDDO !jtaux
      CALL IGPVAL
      dc6i(IAT)=dc6ij(IAT,IAT)*c6ij*C9*FDUM*3._FLOAT+dc6i(IAT)
      ENDDO !iat
      disp=disp-eabc
      endif !.not.noabc
      CALL GSUM(DC6I,N)
      DO IAT=2,N
      RCOVI=RCOV(IZ(IAT))
      DC6IIJ=DC6I(IAT)
      DO JAT=1,IAT-1
      ITASK=ITASK+1
      IF(MPOINT.NE.ITASK)CYCLE
      DO I=1,3
      TUX(I)=LAT(I,1)*REPV1
      TUY0(I)=LAT(I,2)*REPV2
      TUZ0(I)=LAT(I,3)*REPV3
      ENDDO
      TUX(1)=XA(1,JAT)-XA(1,IAT)+TUX(1)
      TUY0(2)=XA(2,JAT)-XA(2,IAT)+TUY0(2)
      TUZ0(3)=XA(3,JAT)-XA(3,IAT)+TUZ0(3)
       ahah=REAL(MREPV1,FLOAT)
      DC6IJI=DC6I(JAT)+DC6IIJ
      RCOVIJ=(RCOV(IZ(JAT))+RCOVI)*PARK1
      do Jtaux=MREPV1,IREPV1
      TUY(1:3)=TUY0(1:3)+TUX(1:3)
      do Jtauy=MREPV2,IREPV2
      TAU(1:3)=TUZ0(1:3)+TUY(1:3)
      do Jtauz=MREPV3,IREPV3
      TAU(1:3)=LAT(1:3,3)+TAU(1:3)
      R20=SUM(TAU(1:3)**2)
      IF(R20.LE.HALF.OR.R20.GE.CRIT_ABC)CYCLE
      R=1._FLOAT/SQRT(R20)
      EXPTERM=EXP(RCOVIJ*R-PARK1)
      DCNII=RCOVIJ*EXPTERM/(R20*(EXPTERM+1._FLOAT)*(EXPTERM+1._FLOAT))
      VEC(1:3)=DC6IJI*DCNII*R*TAU(1:3)
      G(1:3,IAT)=G(1:3,IAT)+VEC(1:3)
      G(1:3,JAT)=G(1:3,JAT)-VEC(1:3)
      DO I=1,3
      SIGMA(1:3,I)=VEC(1:3)*TAU(I)+SIGMA(1:3,I)
      ENDDO
      ENDDO
      TUY(1:3)=LAT(1:3,2)+TUY(1:3)
      ENDDO
      TUX(1:3)=LAT(1:3,1)+TUX(1:3)
      ENDDO
      CALL IGPVAL
      ENDDO
      ENDDO
      DO IAT=1,N
      ITASK=ITASK+1
      IF(MPOINT.NE.ITASK)CYCLE
      RCOVIJ=RCOV(IZ(IAT))*(PARK1+PARK1)
      DC6IIJ=DC6I(IAT)
      DO I=1,3
      TUX(I)=LAT(I,1)*REPV1
      TUY0(I)=LAT(I,2)*REPV2
      TUZ0(I)=LAT(I,3)*REPV3
      ENDDO
      do Jtaux=MREPV1,IREPV1
      TUY(1:3)=TUY0(1:3)+TUX(1:3)
      do Jtauy=MREPV2,IREPV2
      TAU(1:3)=TUZ0(1:3)+TUY(1:3)
      do Jtauz=MREPV3,IREPV3
      TAU(1:3)=LAT(1:3,3)+TAU(1:3)
      IF(ABS(JTAUX)+ABS(JTAUY)+ABS(JTAUZ).EQ.0)CYCLE
      R20=SUM(TAU(1:3)**2)
      IF(R20.GE.CRIT_ABC)CYCLE
      R=1._FLOAT/SQRT(R20)
      EXPTERM=EXP(RCOVIJ*R-PARK1)
      DCNII=RCOVIJ*EXPTERM/(R20*(EXPTERM+1._FLOAT)*(EXPTERM+1._FLOAT))
      VEC(1:3)=DC6IIJ*DCNII*R*TAU(1:3)
      DO I=1,3
      SIGMA(1:3,I)=VEC(1:3)*TAU(I)+SIGMA(1:3,I)
      ENDDO
      ENDDO
      TUY(1:3)=LAT(1:3,2)+TUY(1:3)
      ENDDO
      TUX(1:3)=LAT(1:3,1)+TUX(1:3)
      ENDDO
      CALL IGPVAL
      ENDDO
      CALL IGPRST
      CALL GSUM(G,N*3)
      CALL GSUM(SIGMA,9)
      CALL GSUM(DISP,1)
      IF(.NOT.noabc) THEN
      CALL GSUM(EABC,1)
      ENDIF
      call MINV3(lat,lat_1,STRESS)
      stress=0._FLOAT
      CALL MXMBN(SIGMA,1,3,LAT_1,3,1,STRESS,1,3,3,3,3)
      CALL CRYDEALLOC(DC6I,'PBCGDISP','DC6I')
      CALL CRYDEALLOC(CN,'PBCGDISP','CN')
      CALL CRYDEALLOC(C6SAVE,'PBCGDISP','C6SAVE')
      CALL CRYDEALLOC(DC6IJ,'PBCGDISP','DC6IJ')
      endif ! nversion
      ENDIF
      gnorm=sum(abs(g(1:3,1:n)))
      RETURN
      end subroutine pbcgdisp

      subroutine pbcncoord(natoms,rcov,iz,cn,lat,rep_cn,crit_cn)
      USE NUMBERS
      USE PARINF_MODULE
      USE BASATO_MODULE,ONLY:XA
      IMPLICIT REAL(FLOAT)(A-H,O-Z)
      PARAMETER(PARK1=-16._FLOAT)
      PARAMETER(MAX_ELEM=94)
      REAL(FLOAT),intent(in)  :: rcov(MAX_ELEM),crit_cn
      REAL(FLOAT) :: cn(*),lat(3,3)
      REAL(FLOAT),DIMENSION(3) :: TAU
      integer,intent(in) :: natoms,iz(*)
      integer rep_cn(3),taux,tauy,tauz
      do i=1,natoms
      xn=0._FLOAT
      do iat=1,natoms
        do taux=-rep_cn(1),rep_cn(1)
         do tauy=-rep_cn(2),rep_cn(2)
          do tauz=-rep_cn(3),rep_cn(3)
            if(iat.eq.i .and. taux.eq.0 .and. tauy.eq.0 .and. 
     .       tauz.eq.0)        cycle
            tau=taux*lat(:,1)+tauy*lat(:,2)+tauz*lat(:,3)
            dx=XA(1,iat)-XA(1,i)+tau(1)
            dy=XA(2,iat)-XA(2,i)+tau(2)
            dz=XA(3,iat)-XA(3,i)+tau(3)
            r=(dx*dx+dy*dy+dz*dz)
            if (r.gt.crit_cn) cycle
            r=sqrt(r)
            rco=rcov(iz(i))+rcov(iz(iat))
            rr=rco/r
            damp=1._FLOAT/(1._FLOAT+exp(PARK1*(rr-1._FLOAT)))
            xn=xn+damp

          enddo !tauz
         enddo !tauy
        enddo !taux
      enddo !iat
      cn(i)=xn  
      enddo !i
      end subroutine pbcncoord

      integer function lin(i1,i2)
      integer i1,i2,idum1,idum2
      idum1=max(i1,i2)
      idum2=min(i1,i2)
      lin=idum2+idum1*(idum1-1)/2
      return
      end function lin


      subroutine getc6(maxc,max_elem,c6ab,mxc,iat,jat,nci,ncj,c6)
      USE NUMBERS
      IMPLICIT REAL(FLOAT)(A-H,O-Z)
      PARAMETER(PARK3=-4._FLOAT)
      real(FLOAT) ::  c6ab(max_elem,max_elem,maxc,maxc,3),nci,ncj
      integer mxc(max_elem)
      c6mem=-1.E99_FLOAT
      rsum=0._FLOAT
      csum=0._FLOAT
      c6  =0._FLOAT
      r_save=1.E90_FLOAT
      do i=1,mxc(iat)
      do j=1,mxc(jat)
         c6=c6ab(iat,jat,i,j,1)
         if(c6.gt.0._FLOAT)then
            cn1=c6ab(iat,jat,i,j,2)
            cn2=c6ab(iat,jat,i,j,3)
            r=(cn1-nci)**2+(cn2-ncj)**2
            if (r.lt.r_save) then
               r_save=r
               c6mem=c6
            endif
            tmp1=exp(PARK3*r)
            rsum=rsum+tmp1     
            csum=csum+tmp1*c6
         endif
      enddo
      enddo

      if(rsum.gt.1.E-99_FLOAT)then
         c6=csum/rsum
      else
         c6=c6mem
      endif

      end subroutine getc6


      subroutine get_dC6_dCNij(maxc,max_elem,c6ab,mxci,mxcj,cni,cnj,
     .           izi,izj,iat,jat,c6check,dc6i,dc6j)

      USE NUMBERS
      IMPLICIT REAL(FLOAT)(A-H,O-Z)
      PARAMETER(PARK3=-4._FLOAT)
      real(FLOAT) :: c6ab(max_elem,max_elem,maxc,maxc,3),nenner
      !czw
      integer :: a,b

      c6mem=-1.E99_FLOAT
      r_save=9999._FLOAT
      zaehler=0._FLOAT
      nenner=0._FLOAT

      dzaehler_i=0._FLOAT
      dnenner_i=0._FLOAT
      dzaehler_j=0._FLOAT
      dnenner_j=0._FLOAT


      DO a=1,mxci
        DO b=1,mxcj
          c6ref=c6ab(izi,izj,a,b,1)
          if(c6ref.gt.0._FLOAT)then
            cn_refi=c6ab(izi,izj,a,b,2)
            cn_refj=c6ab(izi,izj,a,b,3)
            r=(cn_refi-cni)*(cn_refi-cni)+(cn_refj-cnj)*(cn_refj-cnj)
            if (r.lt.r_save) then
               r_save=r
               c6mem=c6ref
            endif
            expterm=exp(PARK3*r)
            zaehler=zaehler+c6ref*expterm
            nenner=nenner+expterm
            dzaehler_i=dzaehler_i+c6ref*expterm*
     .             2._FLOAT*PARK3*(cni-cn_refi)
            dnenner_i=dnenner_i+expterm*
     .             2._FLOAT*PARK3*(cni-cn_refi)

            dzaehler_j=dzaehler_j+c6ref*expterm*
     .             2._FLOAT*PARK3*(cnj-cn_refj)
            dnenner_j=dnenner_j+expterm*
     .             2._FLOAT*PARK3*(cnj-cn_refj)
          endif
        ENDDO !b
      ENDDO !a

      if(nenner.gt.1.E-99_FLOAT)then
        c6check=zaehler/nenner
        dc6i=((dzaehler_i*nenner)-(dnenner_i*zaehler))
     .    /(nenner*nenner)
        dc6j=((dzaehler_j*nenner)-(dnenner_j*zaehler))
     .    /(nenner*nenner)
      else
        c6check=c6mem
        dc6i=0._FLOAT
        dc6j=0._FLOAT
      endif
      end subroutine get_dC6_dCNij

      subroutine setfuncpar(func,nversion,s6,rs6,s18,rs18,alp)
      USE NUMBERS
      IMPLICIT REAL(FLOAT)(A-H,O-Z)
      character*(*) func     
      if(nversion.eq.4)then
      s6=1._FLOAT
      alp =14._FLOAT
      select case (func)
         case ("b-p")
              rs6 =0.3946_FLOAT
              s18 =3.2822_FLOAT
              rs18=4.8516_FLOAT
         case ("b-lyp")
              rs6 =0.4298_FLOAT
              s18 =2.6996_FLOAT
              rs18=4.2359_FLOAT
         case ("revpbe")
              rs6 =0.5238_FLOAT
              s18 =2.3550_FLOAT
              rs18=3.5016_FLOAT
         case ("rpbe")
              rs6 =0.1820_FLOAT
              s18 =0.8318_FLOAT
              rs18=4.0094_FLOAT
         case ("b97-d")
              rs6 =0.5545_FLOAT
              s18 =2.2609_FLOAT
              rs18=3.2297_FLOAT
         case ("pbe")
              rs6 =0.4289_FLOAT
              s18 =0.7875_FLOAT
              rs18=4.4407_FLOAT
         case ("rpw86-pbe")
              rs6 =0.4613_FLOAT
              s18 =1.3845_FLOAT
              rs18=4.5062_FLOAT
         case ("b3-lyp")
              rs6 =0.3981_FLOAT
              s18 =1.9889_FLOAT
              rs18=4.4211_FLOAT
         case ("tpss")
              rs6 =0.4535_FLOAT
              s18 =1.9435_FLOAT
              rs18=4.4752_FLOAT
         case ("hf")
              rs6 =0.3385_FLOAT
              s18 =0.9171_FLOAT
              rs18=2.8830_FLOAT
         case ("tpss0")
              rs6 =0.3768_FLOAT
              s18 =1.2576_FLOAT
              rs18=4.5865_FLOAT
         case ("pbe0")
              rs6 =0.4145_FLOAT
              s18 =1.2177_FLOAT
              rs18=4.8593_FLOAT
         case ("revpbe38")
              rs6 =0.4309_FLOAT
              s18 =1.4760_FLOAT
              rs18=3.9446_FLOAT
         case ("pw6b95")
              rs6 =0.2076_FLOAT
              s18 =0.7257_FLOAT
              rs18=6.3750_FLOAT
         case ("b2-plyp")
              rs6 =0.3065_FLOAT
              s18 =0.9147_FLOAT
              rs18=5.0570_FLOAT
              s6=0.64_FLOAT
         case ("dsd-blyp")
              rs6 =0._FLOAT
              s18 =0.2130_FLOAT
              rs18=6.0519_FLOAT
              s6=0.5_FLOAT
         case ("dsd-blyp-fc")
              rs6 =0.0009_FLOAT
              s18 =0.2112_FLOAT
              rs18=5.9807_FLOAT
              s6=0.5_FLOAT
         case ("bop")
              rs6 =0.4870_FLOAT
              s18 =3.2950_FLOAT
              rs18=3.5043_FLOAT
         case ("mpwlyp")
              rs6 =0.4831_FLOAT
              s18 =2.0077_FLOAT
              rs18=4.5323_FLOAT
         case ("o-lyp")
              rs6 =0.5299_FLOAT
              s18 =2.6205_FLOAT
              rs18=2.8065_FLOAT
         case ("pbesol")
              rs6 =0.4466_FLOAT
              s18 =2.9491_FLOAT
              rs18=6.1742_FLOAT
         case ("bpbe")
              rs6 =0.4567_FLOAT
              s18 =4.0728_FLOAT
              rs18=4.3908_FLOAT
         case ("opbe")
              rs6 =0.5512_FLOAT
              s18 =3.3816_FLOAT
              rs18=2.9444_FLOAT
         case ("ssb")
              rs6 =-0.0952_FLOAT
              s18 =-0.1744_FLOAT
              rs18=5.2170_FLOAT
         case ("revssb")
              rs6 =0.4720_FLOAT
              s18 =0.4389_FLOAT
              rs18=4.0986_FLOAT
         case ("otpss")
              rs6 =0.4634_FLOAT
              s18 =2.7495_FLOAT
              rs18=4.3153_FLOAT
        case ("b3pw91")
              rs6 =0.4312_FLOAT
              s18 =2.8524_FLOAT
              rs18=4.4693_FLOAT
         case ("bh-lyp")
              rs6 =0.2793_FLOAT
              s18 =1.0354_FLOAT
              rs18=4.9615_FLOAT
         case ("revpbe0")
              rs6 =0.4679_FLOAT
              s18 =1.7588_FLOAT
              rs18=3.7619_FLOAT
         case ("tpssh")
              rs6 =0.4529_FLOAT
              s18 =2.2382_FLOAT
              rs18=4.6550_FLOAT
         case ("mpw1b95")
              rs6 =0.1955_FLOAT
              s18 =1.0508_FLOAT
              rs18=6.4177_FLOAT
         case ("pwb6k")
              rs6 =0.1805_FLOAT
              s18 =0.9383_FLOAT
              rs18=7.7627_FLOAT
         case ("b1b95")
              rs6 =0.2092_FLOAT
              s18 =1.4507_FLOAT
              rs18=5.5545_FLOAT
         case ("mpwb1k")
              rs6 =0.1474_FLOAT
              s18 =0.9499_FLOAT
              rs18=6.6223_FLOAT
         case ("bmk")
              rs6 =0.1940_FLOAT
              s18 =2.0860_FLOAT
              rs18=5.9197_FLOAT
         case ("cam-b3lyp")
              rs6 =0.3708_FLOAT
              s18 =2.0674_FLOAT
              rs18=5.4743_FLOAT
         case ("lc-wpbe")
              rs6 =0.3919_FLOAT
              s18 =1.8541_FLOAT
              rs18=5.0897_FLOAT
         case ("b2gp-plyp")
              rs6 =0._FLOAT
              s18 =0.2597_FLOAT
              rs18=6.3332_FLOAT
                s6=0.56_FLOAT
         case ("ptpss")
              rs6 =0._FLOAT
              s18 =0.2804_FLOAT
              rs18=6.5745_FLOAT
                s6=0.750_FLOAT
         case ("pwpb95")
              rs6 =0._FLOAT
              s18 =0.2904_FLOAT
              rs18=7.3141_FLOAT
                s6=0.82_FLOAT
      case ("pw1pw")
              rs6 =0.3807_FLOAT
              s18 =2.3363_FLOAT
              rs18=5.8844_FLOAT
       case ("pwgga")
              rs6 =0.2211_FLOAT
              s18 =2.6910_FLOAT
              rs18=6.7278_FLOAT
       case ("hse06")
              rs6 =0.3820_FLOAT
              s18 =2.3100_FLOAT
              rs18=5.6850_FLOAT
       case ("hsesol")
              rs6 =0.4650_FLOAT
              s18 =2.9215_FLOAT
              rs18=6.2003_FLOAT
       case ("scan")
              rs6 =0.5380_FLOAT
              s18 =0._FLOAT
              rs18=5.4200_FLOAT
       case ("rscan")
              rs6 =0.4702_FLOAT
              s18 =1.0886_FLOAT
              rs18=5.7341_FLOAT
       case ("r2scan")
              rs6 =0.4948_FLOAT
              s18 =0.7898_FLOAT
              rs18=5.7308_FLOAT
       case ("r2scan0")
              rs6 =0.4534_FLOAT
              s18 =1.1846_FLOAT
              rs18=5.8972_FLOAT
       case ("r2scanh")
              rs6 =0.4709_FLOAT
              s18 =1.1236_FLOAT
              rs18=5.9157_FLOAT
       case ("r2scan50")
              rs6 =0.4311_FLOAT
              s18 =1.3294_FLOAT
              rs18=5.9240_FLOAT
       case ("mn15")
              rs6 =2.0971_FLOAT
              s18 =0.7862_FLOAT
              rs18=7.5923_FLOAT
         case ("hf/mixed")
              rs6 =0.5607_FLOAT
              s18 =3.9027_FLOAT
              rs18=4.5622_FLOAT
         case ("hf/sv")
              rs6 =0.4249_FLOAT
              s18 =2.1849_FLOAT
              rs18=4.2783_FLOAT
         case ("hf/minis")
              rs6 =0.1702_FLOAT
              s18 =0.9841_FLOAT
              rs18=3.8506_FLOAT
         case ("b3-lyp/6-31gd")
              rs6 =0.5014_FLOAT
              s18 =4.0672_FLOAT
              rs18=4.8409_FLOAT
         case ("hcth120")
              rs6=0.3563_FLOAT
              s18=1.0821_FLOAT
              rs18=4.3359_FLOAT
         case ("hf3c")
              rs6=0.4171_FLOAT
              s18=0.8777_FLOAT
              rs18=2.9149_FLOAT
         case ("hf3cv")
              rs6=0.3063_FLOAT
              s18=0.5022_FLOAT
              rs18=3.9856_FLOAT
         case ("hfsol3c")
              rs6=0.4171_FLOAT
              s18=0.236979_FLOAT
              rs18=2.9149_FLOAT
         case("dftb")
            rs6=0.5719_FLOAT
            s18=0.5883_FLOAT
            rs18=3.6017_FLOAT
         case("pbeh2c","pbeh3c")
            rs6=0.4860_FLOAT
            s18=0.0000_FLOAT
            rs18=4.5_FLOAT
         case("hse3c")
            rs6=0.4411_FLOAT
            s18=0.0000_FLOAT
            rs18=4.5182_FLOAT
         case("b973c")            
            rs6=0.3700_FLOAT
            s18=1.5000_FLOAT
            rs18=4.1000_FLOAT
         case("pbesol03c")
          rs6=0.53633_FLOAT
          s18=0.0000_FLOAT
          rs18=4.64499_FLOAT
         case("hsesol3c")
          rs6=0.51958_FLOAT
          s18=0.0000_FLOAT
          rs18=4.93854_FLOAT
         case("r2scansol3c")
          rs6=0.63937_FLOAT
          s18=0.0000_FLOAT
          rs18=5.94659_FLOAT
         case("r2scan0sol3c")
          rs6=0.57478_FLOAT
          s18=0.0000_FLOAT
          rs18=5.74152_FLOAT
         case("r2scanpob3c")
          rs6=0.91196_FLOAT
          s18=0.0000_FLOAT
          rs18=7.62991_FLOAT
         case("r2scan0pob3c")
          rs6=0.77939_FLOAT
          s18=0.0000_FLOAT
          rs18=6.93380_FLOAT
         case DEFAULT
              call ERRVRS(0,'SETFUNCPAR','functional name unknown')
      end select
      endif

      if(nversion.eq.3)then
      s6  =1._FLOAT
      alp =14._FLOAT
      rs18=1._FLOAT
      select case (func)
         case ("slater-dirac-exchange")
              rs6 =0.999_FLOAT
              s18 =-1.957_FLOAT
              rs18=0.697_FLOAT
         case ("b-lyp")
              rs6=1.094_FLOAT
              s18=1.682_FLOAT
         case ("b-p")
              rs6=1.139_FLOAT
              s18=1.683_FLOAT
         case ("b97-d")
              rs6=0.892_FLOAT
              s18=0.909_FLOAT
         case ("revpbe")
              rs6=0.923_FLOAT
              s18=1.010_FLOAT
         case ("pbe")
              rs6=1.217_FLOAT
              s18=0.722_FLOAT
         case ("pbesol")
              rs6=1.345_FLOAT
              s18=0.612_FLOAT
         case ("rpw86-pbe")
              rs6=1.224_FLOAT
              s18=0.901_FLOAT
         case ("rpbe")
              rs6=0.872_FLOAT
              s18=0.514_FLOAT
         case ("tpss")
              rs6=1.166_FLOAT
              s18=1.105_FLOAT
         case ("b3-lyp")
              rs6=1.261_FLOAT
              s18=1.703_FLOAT
         case ("pbe0")
              rs6=1.287_FLOAT
              s18=0.928_FLOAT
         case ("revpbe38")
              rs6=1.021_FLOAT
              s18=0.862_FLOAT
         case ("pw6b95")
              rs6=1.532_FLOAT
              s18=0.862_FLOAT
         case ("tpss0")
              rs6=1.252_FLOAT
              s18=1.242_FLOAT
         case ("b2-plyp")
              rs6=1.427_FLOAT
              s18=1.022_FLOAT
              s6=0.64_FLOAT
         case ("pwpb95")
              rs6=1.557_FLOAT
              s18=0.705_FLOAT
              s6=0.82_FLOAT
         case ("b2gp-plyp")
              rs6=1.586_FLOAT
              s18=0.760_FLOAT
              s6=0.56_FLOAT
         case ("ptpss")
              rs6=1.541_FLOAT
              s18=0.879_FLOAT
              s6=0.75_FLOAT
         case ("hf")
              rs6=1.158_FLOAT
              s18=1.746_FLOAT
         case ("mpwlyp")
              rs6=1.239_FLOAT
              s18=1.098_FLOAT
         case ("bpbe")
              rs6=1.087_FLOAT
              s18=2.033_FLOAT
         case ("bh-lyp")
              rs6=1.370_FLOAT
              s18=1.442_FLOAT
         case ("tpssh")
              rs6=1.223_FLOAT
              s18=1.219_FLOAT
         case ("pwb6k")
              rs6=1.660_FLOAT
              s18=0.550_FLOAT
         case ("b1b95")
              rs6=1.613_FLOAT
              s18=1.868_FLOAT
         case ("bop")
              rs6=0.929_FLOAT
              s18=1.975_FLOAT
         case ("o-lyp")
              rs6=0.806_FLOAT
              s18=1.764_FLOAT
         case ("o-pbe")
              rs6=0.837_FLOAT
              s18=2.055_FLOAT
         case ("ssb")
              rs6=1.215_FLOAT
              s18=0.663_FLOAT
         case ("revssb")
              rs6=1.221_FLOAT
              s18=0.560_FLOAT
         case ("otpss")
              rs6=1.128_FLOAT
              s18=1.494_FLOAT
         case ("b3pw91")
              rs6=1.176_FLOAT
              s18=1.775_FLOAT
         case ("revpbe0")
              rs6=0.949_FLOAT
              s18=0.792_FLOAT
         case ("pbe38")
              rs6=1.333_FLOAT
              s18=0.998_FLOAT
         case ("mpw1b95")
              rs6=1.605_FLOAT
              s18=1.118_FLOAT
         case ("mpwb1k")
              rs6=1.671_FLOAT
              s18=1.061_FLOAT
         case ("bmk")
              rs6=1.931_FLOAT
              s18=2.168_FLOAT
         case ("cam-b3lyp")
              rs6=1.378_FLOAT
              s18=1.217_FLOAT
         case ("lc-wpbe")
              rs6=1.355_FLOAT
              s18=1.279_FLOAT
         case ("m05")
              rs6=1.373_FLOAT
              s18=0.595_FLOAT
         case ("m052x")
              rs6=1.417_FLOAT
              s18=0._FLOAT
         case ("m06l")
              rs6=1.581_FLOAT
              s18=0._FLOAT
         case ("m06")
              rs6=1.325_FLOAT
              s18=0._FLOAT
         case ("m062x")
              rs6=1.619_FLOAT
              s18=0._FLOAT
         case ("m06hf")
              rs6=1.446_FLOAT
              s18=0._FLOAT
         case ("dftb")
              rs6=1.699_FLOAT
              s18=1.504_FLOAT
         case ("hcth120")
              rs6=1.221_FLOAT
              s18=1.206_FLOAT
      case ("pw1pw")
              rs6 =1.1324_FLOAT
              s18 =0._FLOAT
      case ("pwgga")
              rs6 =1.1271_FLOAT
              s18 =0._FLOAT
      case ("hse06")
              rs6 =1.1290_FLOAT
              s18 =0.1090_FLOAT
      case ("hsesol")
              rs6 =1.2146_FLOAT
              s18 =0._FLOAT
      case ("scan")
              rs6 =1.3240_FLOAT
              s18 =0._FLOAT
      case ("mn15l")
              rs6 =3.3388_FLOAT
              s18 =0._FLOAT
      case ("hse3c")
              rs6 =1.1269_FLOAT
              s18 =0._FLOAT
      case ("pbeh3c")
              rs6 =1.16_FLOAT
              s18 =0._FLOAT
      case ("hf3c")
              rs6= 0.9500_FLOAT
              s18=0.8877_FLOAT
      case ("b973c")
              rs6 =1.0600_FLOAT
              s18 =1.500_FLOAT
      case ("pbesol03c")
              rs6 =1.1593_FLOAT
              s18 =0._FLOAT
      case ("hsesol3c")
              rs6 =1.0396_FLOAT
              s18 =0._FLOAT
      case ("wb97xd3")
               rs6=1.2810_FLOAT
               s18=1.0000_FLOAT
               rs18=1.0940_FLOAT
         case DEFAULT
              call ERRVRS(0,'SETFUNCPAR','functional name unknown')
      end select
      endif !nversion =3
      if(nversion.eq.2)then
      rs6=1.1_FLOAT
      s18=0._FLOAT
      alp=20._FLOAT
      select case (func)
         case ("b-lyp")
              s6=1.2_FLOAT  
         case ("b-p")
              s6=1.05_FLOAT
         case ("b97-d")
              s6=1.25_FLOAT
         case ("revpbe")
              s6=1.25_FLOAT
         case ("pbe")
              s6=0.75_FLOAT
         case ("tpss")
              s6=1._FLOAT
         case ("b3-lyp")
              s6=1.05_FLOAT
         case ("pbe0")
              s6=0.6_FLOAT
         case ("pw6b95")
              s6=0.5_FLOAT
         case ("tpss0")
              s6=0.85_FLOAT
         case ("b2-plyp")
              s6=0.55_FLOAT
         case ("b2gp-plyp")
              s6=0.4_FLOAT
         case ("dsd-blyp")
              s6=0.41_FLOAT
              alp=60._FLOAT
         case DEFAULT
              call ERRVRS(0,'SETFUNCPAR','functional name unknown')
      end select

      endif
      end subroutine setfuncpar

      SUBROUTINE SET_CRITERIA(rthr,lat,tau_max)
      USE NUMBERS
      IMPLICIT REAL(FLOAT)(A-H,O-Z)
        REAL(FLOAT) :: lat(3,3),tau_max(3),norm1(3),norm2(3),norm3(3)
        real(FLOAT),external :: vectorsize

        r_cutoff=sqrt(rthr)
        call kreuzprodukt(lat(:,2),lat(:,3),norm1)
        call kreuzprodukt(lat(:,3),lat(:,1),norm2)
        call kreuzprodukt(lat(:,1),lat(:,2),norm3)
        norm1=norm1/VECTORSIZE(norm1)
        norm2=norm2/VECTORSIZE(norm2)
        norm3=norm3/VECTORSIZE(norm3)
        cos10=SUM(norm1*lat(:,1))
        cos21=SUM(norm2*lat(:,2))
        cos32=SUM(norm3*lat(:,3))
        tau_max(1)=abs(r_cutoff/cos10)
        tau_max(2)=abs(r_cutoff/cos21)
        tau_max(3)=abs(r_cutoff/cos32)
      END SUBROUTINE SET_CRITERIA

      SUBROUTINE kreuzprodukt(A,B,C)
      USE NUMBERS
      IMPLICIT REAL(FLOAT)(A-H,O-Z)
        REAL(FLOAT),DIMENSION(3) :: A,B,C
        X=A(2)*B(3)-B(2)*A(3)
        Y=A(3)*B(1)-B(3)*A(1)
        Z=A(1)*B(2)-B(1)*A(2)
        C=(/X,Y,Z/)
      END SUBROUTINE kreuzprodukt

       FUNCTION VECTORSIZE(VECT)
      USE NUMBERS
      IMPLICIT REAL(FLOAT)(A-H,O-Z)
         REAL(FLOAT),DIMENSION(3) :: VECT
         VECTORSIZE=SQRT(SUM(VECT(1:3)**2))
       END FUNCTION VECTORSIZE

      subroutine setr0ab(max_elem,r)
      USE NUMBERS
      USE PARINF_MODULE
      IMPLICIT REAL(FLOAT)(A-H,O-Z)
      real(FLOAT) :: r(max_elem,max_elem),T0AB(4465)
      DATA T0AB(1:70)/2.1823_FLOAT,1.8547_FLOAT,1.7347_FLOAT,
     *2.9086_FLOAT,2.5732_FLOAT,3.4956_FLOAT,2.3550_FLOAT,2.5095_FLOAT,
     *2.9802_FLOAT,3.0982_FLOAT,2.5141_FLOAT,2.3917_FLOAT,2.9977_FLOAT,
     *2.9484_FLOAT,3.2160_FLOAT,2.4492_FLOAT,2.2527_FLOAT,3.1933_FLOAT,
     *3.0214_FLOAT,2.9531_FLOAT,2.9103_FLOAT,2.3667_FLOAT,2.1328_FLOAT,
     *2.8784_FLOAT,2.7660_FLOAT,2.7776_FLOAT,2.7063_FLOAT,2.6225_FLOAT,
     *2.1768_FLOAT,2.0625_FLOAT,2.6395_FLOAT,2.6648_FLOAT,2.6482_FLOAT,
     *2.5697_FLOAT,2.4846_FLOAT,2.4817_FLOAT,2.0646_FLOAT,1.9891_FLOAT,
     *2.5086_FLOAT,2.6908_FLOAT,2.6233_FLOAT,2.4770_FLOAT,2.3885_FLOAT,
     *2.3511_FLOAT,2.2996_FLOAT,1.9892_FLOAT,1.9251_FLOAT,2.4190_FLOAT,
     *2.5473_FLOAT,2.4994_FLOAT,2.4091_FLOAT,2.3176_FLOAT,2.2571_FLOAT,
     *2.1946_FLOAT,2.1374_FLOAT,2.9898_FLOAT,2.6397_FLOAT,3.6031_FLOAT,
     *3.1219_FLOAT,3.7620_FLOAT,3.2485_FLOAT,2.9357_FLOAT,2.7093_FLOAT,
     *2.5781_FLOAT,2.4839_FLOAT,3.7082_FLOAT,2.5129_FLOAT,2.7321_FLOAT,
     *3.1052_FLOAT,3.2962_FLOAT/
      DATA T0AB(71:140)/3.1331_FLOAT,3.2000_FLOAT,2.9586_FLOAT,
     *3.0822_FLOAT,2.8582_FLOAT,2.7120_FLOAT,3.2570_FLOAT,3.4839_FLOAT,
     *2.8766_FLOAT,2.7427_FLOAT,3.2776_FLOAT,3.2363_FLOAT,3.5929_FLOAT,
     *3.2826_FLOAT,3.0911_FLOAT,2.9369_FLOAT,2.9030_FLOAT,2.7789_FLOAT,
     *3.3921_FLOAT,3.3970_FLOAT,4.0106_FLOAT,2.8884_FLOAT,2.6605_FLOAT,
     *3.7513_FLOAT,3.1613_FLOAT,3.3605_FLOAT,3.3325_FLOAT,3.0991_FLOAT,
     *2.9297_FLOAT,2.8674_FLOAT,2.7571_FLOAT,3.8129_FLOAT,3.3266_FLOAT,
     *3.7105_FLOAT,3.7917_FLOAT,2.8304_FLOAT,2.5538_FLOAT,3.3932_FLOAT,
     *3.1193_FLOAT,3.1866_FLOAT,3.1245_FLOAT,3.0465_FLOAT,2.8727_FLOAT,
     *2.7664_FLOAT,2.6926_FLOAT,3.4608_FLOAT,3.2984_FLOAT,3.5142_FLOAT,
     *3.5418_FLOAT,3.5017_FLOAT,2.6190_FLOAT,2.4797_FLOAT,3.1331_FLOAT,
     *3.0540_FLOAT,3.0651_FLOAT,2.9879_FLOAT,2.9054_FLOAT,2.8805_FLOAT,
     *2.7330_FLOAT,2.6331_FLOAT,3.2096_FLOAT,3.5668_FLOAT,3.3684_FLOAT,
     *3.3686_FLOAT,3.3180_FLOAT,3.3107_FLOAT,2.4757_FLOAT,2.4019_FLOAT,
     *2.9789_FLOAT,3.1468_FLOAT/
      DATA T0AB(141:210)/2.9768_FLOAT,2.8848_FLOAT,2.7952_FLOAT,
     *2.7457_FLOAT,2.6881_FLOAT,2.5728_FLOAT,3.0574_FLOAT,3.3264_FLOAT,
     *3.3562_FLOAT,3.2529_FLOAT,3.1916_FLOAT,3.1523_FLOAT,3.1046_FLOAT,
     *2.3725_FLOAT,2.3289_FLOAT,2.8760_FLOAT,2.9804_FLOAT,2.9093_FLOAT,
     *2.8040_FLOAT,2.7071_FLOAT,2.6386_FLOAT,2.5720_FLOAT,2.5139_FLOAT,
     *2.9517_FLOAT,3.1606_FLOAT,3.2085_FLOAT,3.1692_FLOAT,3.0982_FLOAT,
     *3.0352_FLOAT,2.9730_FLOAT,2.9148_FLOAT,3.2147_FLOAT,2.8315_FLOAT,
     *3.8724_FLOAT,3.4621_FLOAT,3.8823_FLOAT,3.3760_FLOAT,3.0746_FLOAT,
     *2.8817_FLOAT,2.7552_FLOAT,2.6605_FLOAT,3.9740_FLOAT,3.6192_FLOAT,
     *3.6569_FLOAT,3.9586_FLOAT,3.6188_FLOAT,3.3917_FLOAT,3.2479_FLOAT,
     *3.1434_FLOAT,4.2411_FLOAT,2.7597_FLOAT,3.0588_FLOAT,3.3474_FLOAT,
     *3.6214_FLOAT,3.4353_FLOAT,3.4729_FLOAT,3.2487_FLOAT,3.3200_FLOAT,
     *3.0914_FLOAT,2.9403_FLOAT,3.4972_FLOAT,3.7993_FLOAT,3.6773_FLOAT,
     *3.8678_FLOAT,3.5808_FLOAT,3.8243_FLOAT,3.5826_FLOAT,3.4156_FLOAT,
     *3.8765_FLOAT,4.1035_FLOAT/
      DATA T0AB(211:280)/2.7361_FLOAT,2.9765_FLOAT,3.2475_FLOAT,
     *3.5004_FLOAT,3.4185_FLOAT,3.4378_FLOAT,3.2084_FLOAT,3.2787_FLOAT,
     *3.0604_FLOAT,2.9187_FLOAT,3.4037_FLOAT,3.6759_FLOAT,3.6586_FLOAT,
     *3.8327_FLOAT,3.5372_FLOAT,3.7665_FLOAT,3.5310_FLOAT,3.3700_FLOAT,
     *3.7788_FLOAT,3.9804_FLOAT,3.8903_FLOAT,2.6832_FLOAT,2.9060_FLOAT,
     *3.2613_FLOAT,3.4359_FLOAT,3.3538_FLOAT,3.3860_FLOAT,3.1550_FLOAT,
     *3.2300_FLOAT,3.0133_FLOAT,2.8736_FLOAT,3.4024_FLOAT,3.6142_FLOAT,
     *3.5979_FLOAT,3.5295_FLOAT,3.4834_FLOAT,3.7140_FLOAT,3.4782_FLOAT,
     *3.3170_FLOAT,3.7434_FLOAT,3.9623_FLOAT,3.8181_FLOAT,3.7642_FLOAT,
     *2.6379_FLOAT,2.8494_FLOAT,3.1840_FLOAT,3.4225_FLOAT,3.2771_FLOAT,
     *3.3401_FLOAT,3.1072_FLOAT,3.1885_FLOAT,2.9714_FLOAT,2.8319_FLOAT,
     *3.3315_FLOAT,3.5979_FLOAT,3.5256_FLOAT,3.4980_FLOAT,3.4376_FLOAT,
     *3.6714_FLOAT,3.4346_FLOAT,3.2723_FLOAT,3.6859_FLOAT,3.8985_FLOAT,
     *3.7918_FLOAT,3.7372_FLOAT,3.7211_FLOAT,2.9230_FLOAT,2.6223_FLOAT,
     *3.4161_FLOAT,2.8999_FLOAT/
      DATA T0AB(281:350)/3.0557_FLOAT,3.3308_FLOAT,3.0555_FLOAT,
     *2.8508_FLOAT,2.7385_FLOAT,2.6640_FLOAT,3.5263_FLOAT,3.0277_FLOAT,
     *3.2990_FLOAT,3.7721_FLOAT,3.5017_FLOAT,3.2751_FLOAT,3.1368_FLOAT,
     *3.0435_FLOAT,3.7873_FLOAT,3.2858_FLOAT,3.2140_FLOAT,3.1727_FLOAT,
     *3.2178_FLOAT,3.4414_FLOAT,2.5490_FLOAT,2.7623_FLOAT,3.0991_FLOAT,
     *3.3252_FLOAT,3.1836_FLOAT,3.2428_FLOAT,3.0259_FLOAT,3.1225_FLOAT,
     *2.9032_FLOAT,2.7621_FLOAT,3.2490_FLOAT,3.5110_FLOAT,3.4429_FLOAT,
     *3.3845_FLOAT,3.3574_FLOAT,3.6045_FLOAT,3.3658_FLOAT,3.2013_FLOAT,
     *3.6110_FLOAT,3.8241_FLOAT,3.7090_FLOAT,3.6496_FLOAT,3.6333_FLOAT,
     *3.0896_FLOAT,3.5462_FLOAT,2.4926_FLOAT,2.7136_FLOAT,3.0693_FLOAT,
     *3.2699_FLOAT,3.1272_FLOAT,3.1893_FLOAT,2.9658_FLOAT,3.0972_FLOAT,
     *2.8778_FLOAT,2.7358_FLOAT,3.2206_FLOAT,3.4566_FLOAT,3.3896_FLOAT,
     *3.3257_FLOAT,3.2946_FLOAT,3.5693_FLOAT,3.3312_FLOAT,3.1670_FLOAT,
     *3.5805_FLOAT,3.7711_FLOAT,3.6536_FLOAT,3.5927_FLOAT,3.5775_FLOAT,
     *3.0411_FLOAT,3.4885_FLOAT/
      DATA T0AB(351:420)/3.4421_FLOAT,2.4667_FLOAT,2.6709_FLOAT,
     *3.0575_FLOAT,3.2357_FLOAT,3.0908_FLOAT,3.1537_FLOAT,2.9235_FLOAT,
     *3.0669_FLOAT,2.8476_FLOAT,2.7054_FLOAT,3.2064_FLOAT,3.4519_FLOAT,
     *3.3593_FLOAT,3.2921_FLOAT,3.2577_FLOAT,3.2161_FLOAT,3.2982_FLOAT,
     *3.1339_FLOAT,3.5606_FLOAT,3.7582_FLOAT,3.6432_FLOAT,3.5833_FLOAT,
     *3.5691_FLOAT,3.0161_FLOAT,3.4812_FLOAT,3.4339_FLOAT,3.4327_FLOAT,
     *2.4515_FLOAT,2.6338_FLOAT,3.0511_FLOAT,3.2229_FLOAT,3.0630_FLOAT,
     *3.1265_FLOAT,2.8909_FLOAT,3.0253_FLOAT,2.8184_FLOAT,2.6764_FLOAT,
     *3.1968_FLOAT,3.4114_FLOAT,3.3492_FLOAT,3.2691_FLOAT,3.2320_FLOAT,
     *3.1786_FLOAT,3.2680_FLOAT,3.1036_FLOAT,3.5453_FLOAT,3.7259_FLOAT,
     *3.6090_FLOAT,3.5473_FLOAT,3.5327_FLOAT,3.0018_FLOAT,3.4413_FLOAT,
     *3.3907_FLOAT,3.3593_FLOAT,3.3462_FLOAT,2.4413_FLOAT,2.6006_FLOAT,
     *3.0540_FLOAT,3.1987_FLOAT,3.0490_FLOAT,3.1058_FLOAT,2.8643_FLOAT,
     *2.9948_FLOAT,2.7908_FLOAT,2.6491_FLOAT,3.1950_FLOAT,3.3922_FLOAT,
     *3.3316_FLOAT,3.2585_FLOAT/
      DATA T0AB(421:490)/3.2136_FLOAT,3.1516_FLOAT,3.2364_FLOAT,
     *3.0752_FLOAT,3.5368_FLOAT,3.7117_FLOAT,3.5941_FLOAT,3.5313_FLOAT,
     *3.5164_FLOAT,2.9962_FLOAT,3.4225_FLOAT,3.3699_FLOAT,3.3370_FLOAT,
     *3.3234_FLOAT,3.3008_FLOAT,2.4318_FLOAT,2.5729_FLOAT,3.0416_FLOAT,
     *3.1639_FLOAT,3.0196_FLOAT,3.0843_FLOAT,2.8413_FLOAT,2.7436_FLOAT,
     *2.7608_FLOAT,2.6271_FLOAT,3.1811_FLOAT,3.3591_FLOAT,3.3045_FLOAT,
     *3.2349_FLOAT,3.1942_FLOAT,3.1291_FLOAT,3.2111_FLOAT,3.0534_FLOAT,
     *3.5189_FLOAT,3.6809_FLOAT,3.5635_FLOAT,3.5001_FLOAT,3.4854_FLOAT,
     *2.9857_FLOAT,3.3897_FLOAT,3.3363_FLOAT,3.3027_FLOAT,3.2890_FLOAT,
     *3.2655_FLOAT,3.2309_FLOAT,2.8502_FLOAT,2.6934_FLOAT,3.2467_FLOAT,
     *3.1921_FLOAT,3.5663_FLOAT,3.2541_FLOAT,3.0571_FLOAT,2.9048_FLOAT,
     *2.8657_FLOAT,2.7438_FLOAT,3.3547_FLOAT,3.3510_FLOAT,3.9837_FLOAT,
     *3.6871_FLOAT,3.4862_FLOAT,3.3389_FLOAT,3.2413_FLOAT,3.1708_FLOAT,
     *3.6096_FLOAT,3.6280_FLOAT,3.6860_FLOAT,3.5568_FLOAT,3.4836_FLOAT,
     *3.2868_FLOAT,3.3994_FLOAT/
      DATA T0AB(491:560)/3.3476_FLOAT,3.3170_FLOAT,3.2950_FLOAT,
     *3.2874_FLOAT,3.2606_FLOAT,3.9579_FLOAT,2.9226_FLOAT,2.6838_FLOAT,
     *3.7867_FLOAT,3.1732_FLOAT,3.3872_FLOAT,3.3643_FLOAT,3.1267_FLOAT,
     *2.9541_FLOAT,2.8505_FLOAT,2.7781_FLOAT,3.8475_FLOAT,3.3336_FLOAT,
     *3.7359_FLOAT,3.8266_FLOAT,3.5733_FLOAT,3.3959_FLOAT,3.2775_FLOAT,
     *3.1915_FLOAT,3.9878_FLOAT,3.8816_FLOAT,3.5810_FLOAT,3.5364_FLOAT,
     *3.5060_FLOAT,3.8097_FLOAT,3.3925_FLOAT,3.3348_FLOAT,3.3019_FLOAT,
     *3.2796_FLOAT,3.2662_FLOAT,3.2464_FLOAT,3.7136_FLOAT,3.8619_FLOAT,
     *2.9140_FLOAT,2.6271_FLOAT,3.4771_FLOAT,3.1774_FLOAT,3.2560_FLOAT,
     *3.1970_FLOAT,3.1207_FLOAT,2.9406_FLOAT,2.8322_FLOAT,2.7571_FLOAT,
     *3.5455_FLOAT,3.3514_FLOAT,3.5837_FLOAT,3.6177_FLOAT,3.5816_FLOAT,
     *3.3902_FLOAT,3.2604_FLOAT,3.1652_FLOAT,3.7037_FLOAT,3.6283_FLOAT,
     *3.5858_FLOAT,3.5330_FLOAT,3.4884_FLOAT,3.5789_FLOAT,3.4094_FLOAT,
     *3.3473_FLOAT,3.3118_FLOAT,3.2876_FLOAT,3.2707_FLOAT,3.2521_FLOAT,
     *3.5570_FLOAT,3.6496_FLOAT/
      DATA T0AB(561:630)/3.6625_FLOAT,2.7300_FLOAT,2.5870_FLOAT,
     *3.2471_FLOAT,3.1487_FLOAT,3.1667_FLOAT,3.0914_FLOAT,3.0107_FLOAT,
     *2.9812_FLOAT,2.8300_FLOAT,2.7284_FLOAT,3.3259_FLOAT,3.3182_FLOAT,
     *3.4707_FLOAT,3.4748_FLOAT,3.4279_FLOAT,3.4182_FLOAT,3.2547_FLOAT,
     *3.1353_FLOAT,3.5116_FLOAT,3.9432_FLOAT,3.8828_FLOAT,3.8303_FLOAT,
     *3.7880_FLOAT,3.3760_FLOAT,3.7218_FLOAT,3.3408_FLOAT,3.3059_FLOAT,
     *3.2698_FLOAT,3.2446_FLOAT,3.2229_FLOAT,3.4422_FLOAT,3.5023_FLOAT,
     *3.5009_FLOAT,3.5268_FLOAT,2.6026_FLOAT,2.5355_FLOAT,3.1129_FLOAT,
     *3.2863_FLOAT,3.1029_FLOAT,3.0108_FLOAT,2.9227_FLOAT,2.8694_FLOAT,
     *2.8109_FLOAT,2.6929_FLOAT,3.1958_FLOAT,3.4670_FLOAT,3.4018_FLOAT,
     *3.3805_FLOAT,3.3218_FLOAT,3.2815_FLOAT,3.2346_FLOAT,3.0994_FLOAT,
     *3.3937_FLOAT,3.7266_FLOAT,3.6697_FLOAT,3.6164_FLOAT,3.5730_FLOAT,
     *3.2522_FLOAT,3.5051_FLOAT,3.4686_FLOAT,3.4355_FLOAT,3.4084_FLOAT,
     *3.3748_FLOAT,3.3496_FLOAT,3.3692_FLOAT,3.4052_FLOAT,3.3910_FLOAT,
     *3.3849_FLOAT,3.3662_FLOAT/
      DATA T0AB(631:700)/2.5087_FLOAT,2.4814_FLOAT,3.0239_FLOAT,
     *3.1312_FLOAT,3.0535_FLOAT,2.9457_FLOAT,2.8496_FLOAT,2.7780_FLOAT,
     *2.7828_FLOAT,2.6532_FLOAT,3.1063_FLOAT,3.3143_FLOAT,3.3549_FLOAT,
     *3.3120_FLOAT,3.2421_FLOAT,3.1787_FLOAT,3.1176_FLOAT,3.0613_FLOAT,
     *3.3082_FLOAT,3.5755_FLOAT,3.5222_FLOAT,3.4678_FLOAT,3.4231_FLOAT,
     *3.1684_FLOAT,3.3528_FLOAT,3.3162_FLOAT,3.2827_FLOAT,3.2527_FLOAT,
     *3.2308_FLOAT,3.2029_FLOAT,3.3173_FLOAT,3.3343_FLOAT,3.3092_FLOAT,
     *3.2795_FLOAT,3.2452_FLOAT,3.2096_FLOAT,3.2893_FLOAT,2.8991_FLOAT,
     *4.0388_FLOAT,3.6100_FLOAT,3.9388_FLOAT,3.4475_FLOAT,3.1590_FLOAT,
     *2.9812_FLOAT,2.8586_FLOAT,2.7683_FLOAT,4.1428_FLOAT,3.7911_FLOAT,
     *3.8225_FLOAT,4.0372_FLOAT,3.7059_FLOAT,3.4935_FLOAT,3.3529_FLOAT,
     *3.2492_FLOAT,4.4352_FLOAT,4.0826_FLOAT,3.9733_FLOAT,3.9254_FLOAT,
     *3.8646_FLOAT,3.9315_FLOAT,3.7837_FLOAT,3.7465_FLOAT,3.7211_FLOAT,
     *3.7012_FLOAT,3.6893_FLOAT,3.6676_FLOAT,3.7736_FLOAT,4.0660_FLOAT,
     *3.7926_FLOAT,3.6158_FLOAT/
      DATA T0AB(701:770)/3.5017_FLOAT,3.4166_FLOAT,4.6176_FLOAT,
     *2.8786_FLOAT,3.1658_FLOAT,3.5823_FLOAT,3.7689_FLOAT,3.5762_FLOAT,
     *3.5789_FLOAT,3.3552_FLOAT,3.4004_FLOAT,3.1722_FLOAT,3.0212_FLOAT,
     *3.7241_FLOAT,3.9604_FLOAT,3.8500_FLOAT,3.9844_FLOAT,3.7035_FLOAT,
     *3.9161_FLOAT,3.6751_FLOAT,3.5075_FLOAT,4.1151_FLOAT,4.2877_FLOAT,
     *4.1579_FLOAT,4.1247_FLOAT,4.0617_FLOAT,3.4874_FLOAT,3.9848_FLOAT,
     *3.9280_FLOAT,3.9079_FLOAT,3.8751_FLOAT,3.8604_FLOAT,3.8277_FLOAT,
     *3.8002_FLOAT,3.9981_FLOAT,3.7544_FLOAT,4.0371_FLOAT,3.8225_FLOAT,
     *3.6718_FLOAT,4.3092_FLOAT,4.4764_FLOAT,2.8997_FLOAT,3.0953_FLOAT,
     *3.4524_FLOAT,3.6107_FLOAT,3.6062_FLOAT,3.5783_FLOAT,3.3463_FLOAT,
     *3.3855_FLOAT,3.1746_FLOAT,3.0381_FLOAT,3.6019_FLOAT,3.7938_FLOAT,
     *3.8697_FLOAT,3.9781_FLOAT,3.6877_FLOAT,3.8736_FLOAT,3.6451_FLOAT,
     *3.4890_FLOAT,3.9858_FLOAT,4.1179_FLOAT,4.0430_FLOAT,3.9563_FLOAT,
     *3.9182_FLOAT,3.4002_FLOAT,3.8310_FLOAT,3.7716_FLOAT,3.7543_FLOAT,
     *3.7203_FLOAT,3.7053_FLOAT/
      DATA T0AB(771:840)/3.6742_FLOAT,3.8318_FLOAT,3.7631_FLOAT,
     *3.7392_FLOAT,3.9892_FLOAT,3.7832_FLOAT,3.6406_FLOAT,4.1701_FLOAT,
     *4.3016_FLOAT,4.2196_FLOAT,2.8535_FLOAT,3.0167_FLOAT,3.3978_FLOAT,
     *3.5363_FLOAT,3.5393_FLOAT,3.5301_FLOAT,3.2960_FLOAT,3.3352_FLOAT,
     *3.1287_FLOAT,2.9967_FLOAT,3.6659_FLOAT,3.7239_FLOAT,3.8070_FLOAT,
     *3.7165_FLOAT,3.6368_FLOAT,3.8162_FLOAT,3.5885_FLOAT,3.4336_FLOAT,
     *3.9829_FLOAT,4.0529_FLOAT,3.9584_FLOAT,3.9025_FLOAT,3.8607_FLOAT,
     *3.3673_FLOAT,3.7658_FLOAT,3.7035_FLOAT,3.6866_FLOAT,3.6504_FLOAT,
     *3.6339_FLOAT,3.6024_FLOAT,3.7708_FLOAT,3.7283_FLOAT,3.6896_FLOAT,
     *3.9315_FLOAT,3.7250_FLOAT,3.5819_FLOAT,4.1457_FLOAT,4.2280_FLOAT,
     *4.1130_FLOAT,4.0597_FLOAT,3.0905_FLOAT,2.7998_FLOAT,3.6448_FLOAT,
     *3.0739_FLOAT,3.2996_FLOAT,3.5262_FLOAT,3.2559_FLOAT,3.0518_FLOAT,
     *2.9394_FLOAT,2.8658_FLOAT,3.7514_FLOAT,3.2295_FLOAT,3.5643_FLOAT,
     *3.7808_FLOAT,3.6931_FLOAT,3.4723_FLOAT,3.3357_FLOAT,3.2429_FLOAT,
     *4.0280_FLOAT,3.5589_FLOAT/
      DATA T0AB(841:910)/3.4636_FLOAT,3.4994_FLOAT,3.4309_FLOAT,
     *3.6177_FLOAT,3.2946_FLOAT,3.2376_FLOAT,3.2050_FLOAT,3.1847_FLOAT,
     *3.1715_FLOAT,3.1599_FLOAT,3.5555_FLOAT,3.8111_FLOAT,3.7693_FLOAT,
     *3.5718_FLOAT,3.4498_FLOAT,3.3662_FLOAT,4.1608_FLOAT,3.7417_FLOAT,
     *3.6536_FLOAT,3.6154_FLOAT,3.8596_FLOAT,3.0301_FLOAT,2.7312_FLOAT,
     *3.5821_FLOAT,3.0473_FLOAT,3.2137_FLOAT,3.4679_FLOAT,3.1975_FLOAT,
     *2.9969_FLOAT,2.8847_FLOAT,2.8110_FLOAT,3.6931_FLOAT,3.2076_FLOAT,
     *3.4943_FLOAT,3.5956_FLOAT,3.6379_FLOAT,3.4190_FLOAT,3.2808_FLOAT,
     *3.1860_FLOAT,3.9850_FLOAT,3.5105_FLOAT,3.4330_FLOAT,3.3797_FLOAT,
     *3.4155_FLOAT,3.6033_FLOAT,3.2737_FLOAT,3.2145_FLOAT,3.1807_FLOAT,
     *3.1596_FLOAT,3.1461_FLOAT,3.1337_FLOAT,3.4812_FLOAT,3.6251_FLOAT,
     *3.7152_FLOAT,3.5201_FLOAT,3.3966_FLOAT,3.3107_FLOAT,4.1128_FLOAT,
     *3.6899_FLOAT,3.6082_FLOAT,3.5604_FLOAT,3.7834_FLOAT,3.7543_FLOAT,
     *2.9189_FLOAT,2.6777_FLOAT,3.4925_FLOAT,2.9648_FLOAT,3.1216_FLOAT,
     *3.2940_FLOAT,3.0975_FLOAT/
      DATA T0AB(911:980)/2.9757_FLOAT,2.8493_FLOAT,2.7638_FLOAT,
     *3.6085_FLOAT,3.1214_FLOAT,3.4006_FLOAT,3.4793_FLOAT,3.5147_FLOAT,
     *3.3806_FLOAT,3.2356_FLOAT,3.1335_FLOAT,3.9144_FLOAT,3.4183_FLOAT,
     *3.3369_FLOAT,3.2803_FLOAT,3.2679_FLOAT,3.4871_FLOAT,3.1714_FLOAT,
     *3.1521_FLOAT,3.1101_FLOAT,3.0843_FLOAT,3.0670_FLOAT,3.0539_FLOAT,
     *3.3890_FLOAT,3.5086_FLOAT,3.5895_FLOAT,3.4783_FLOAT,3.3484_FLOAT,
     *3.2559_FLOAT,4.0422_FLOAT,3.5967_FLOAT,3.5113_FLOAT,3.4576_FLOAT,
     *3.6594_FLOAT,3.6313_FLOAT,3.5690_FLOAT,2.8578_FLOAT,2.6334_FLOAT,
     *3.4673_FLOAT,2.9245_FLOAT,3.0732_FLOAT,3.2435_FLOAT,3.0338_FLOAT,
     *2.9462_FLOAT,2.8143_FLOAT,2.7240_FLOAT,3.5832_FLOAT,3.0789_FLOAT,
     *3.3617_FLOAT,3.4246_FLOAT,3.4505_FLOAT,3.3443_FLOAT,3.1964_FLOAT,
     *3.0913_FLOAT,3.8921_FLOAT,3.3713_FLOAT,3.2873_FLOAT,3.2281_FLOAT,
     *3.2165_FLOAT,3.4386_FLOAT,3.1164_FLOAT,3.1220_FLOAT,3.0761_FLOAT,
     *3.0480_FLOAT,3.0295_FLOAT,3.0155_FLOAT,3.3495_FLOAT,3.4543_FLOAT,
     *3.5260_FLOAT,3.4413_FLOAT/
      DATA T0AB(981:1050)/3.3085_FLOAT,3.2134_FLOAT,4.0170_FLOAT,
     *3.5464_FLOAT,3.4587_FLOAT,3.4006_FLOAT,3.6027_FLOAT,3.5730_FLOAT,
     *3.4945_FLOAT,3.4623_FLOAT,2.8240_FLOAT,2.5960_FLOAT,3.4635_FLOAT,
     *2.9032_FLOAT,3.0431_FLOAT,3.2115_FLOAT,2.9892_FLOAT,2.9148_FLOAT,
     *2.7801_FLOAT,2.6873_FLOAT,3.5776_FLOAT,3.0568_FLOAT,3.3433_FLOAT,
     *3.3949_FLOAT,3.4132_FLOAT,3.3116_FLOAT,3.1616_FLOAT,3.0548_FLOAT,
     *3.8859_FLOAT,3.3719_FLOAT,3.2917_FLOAT,3.2345_FLOAT,3.2274_FLOAT,
     *3.4171_FLOAT,3.1293_FLOAT,3.0567_FLOAT,3.0565_FLOAT,3.0274_FLOAT,
     *3.0087_FLOAT,2.9939_FLOAT,3.3293_FLOAT,3.4249_FLOAT,3.4902_FLOAT,
     *3.4091_FLOAT,3.2744_FLOAT,3.1776_FLOAT,4.0078_FLOAT,3.5374_FLOAT,
     *3.4537_FLOAT,3.3956_FLOAT,3.5747_FLOAT,3.5430_FLOAT,3.4522_FLOAT,
     *3.4160_FLOAT,3.3975_FLOAT,2.8004_FLOAT,2.5621_FLOAT,3.4617_FLOAT,
     *2.9154_FLOAT,3.0203_FLOAT,3.1875_FLOAT,2.9548_FLOAT,2.8038_FLOAT,
     *2.7472_FLOAT,2.6530_FLOAT,3.5736_FLOAT,3.0584_FLOAT,3.3304_FLOAT,
     *3.3748_FLOAT,3.3871_FLOAT/
      DATA T0AB(1051:1120)/3.2028_FLOAT,3.1296_FLOAT,3.0214_FLOAT,
     *3.8796_FLOAT,3.3337_FLOAT,3.2492_FLOAT,3.1883_FLOAT,3.1802_FLOAT,
     *3.4050_FLOAT,3.0756_FLOAT,3.0478_FLOAT,3.0322_FLOAT,3.0323_FLOAT,
     *3.0163_FLOAT,3.0019_FLOAT,3.3145_FLOAT,3.4050_FLOAT,3.4656_FLOAT,
     *3.3021_FLOAT,3.2433_FLOAT,3.1453_FLOAT,3.9991_FLOAT,3.5017_FLOAT,
     *3.4141_FLOAT,3.3520_FLOAT,3.5583_FLOAT,3.5251_FLOAT,3.4243_FLOAT,
     *3.3851_FLOAT,3.3662_FLOAT,3.3525_FLOAT,2.7846_FLOAT,2.5324_FLOAT,
     *3.4652_FLOAT,2.8759_FLOAT,3.0051_FLOAT,3.1692_FLOAT,2.9273_FLOAT,
     *2.7615_FLOAT,2.7164_FLOAT,2.6212_FLOAT,3.5744_FLOAT,3.0275_FLOAT,
     *3.3249_FLOAT,3.3627_FLOAT,3.3686_FLOAT,3.1669_FLOAT,3.0584_FLOAT,
     *2.9915_FLOAT,3.8773_FLOAT,3.3099_FLOAT,3.2231_FLOAT,3.1600_FLOAT,
     *3.1520_FLOAT,3.4023_FLOAT,3.0426_FLOAT,3.0099_FLOAT,2.9920_FLOAT,
     *2.9809_FLOAT,2.9800_FLOAT,2.9646_FLOAT,3.3068_FLOAT,3.3930_FLOAT,
     *3.4486_FLOAT,3.2682_FLOAT,3.1729_FLOAT,3.1168_FLOAT,3.9952_FLOAT,
     *3.4796_FLOAT,3.3901_FLOAT/
      DATA T0AB(1121:1190)/3.3255_FLOAT,3.5530_FLOAT,3.5183_FLOAT,
     *3.4097_FLOAT,3.3683_FLOAT,3.3492_FLOAT,3.3360_FLOAT,3.3308_FLOAT,
     *2.5424_FLOAT,2.6601_FLOAT,3.2555_FLOAT,3.2807_FLOAT,3.1384_FLOAT,
     *3.1737_FLOAT,2.9397_FLOAT,2.8429_FLOAT,2.8492_FLOAT,2.7225_FLOAT,
     *3.3875_FLOAT,3.4910_FLOAT,3.4520_FLOAT,3.3608_FLOAT,3.3036_FLOAT,
     *3.2345_FLOAT,3.2999_FLOAT,3.1487_FLOAT,3.7409_FLOAT,3.8392_FLOAT,
     *3.7148_FLOAT,3.6439_FLOAT,3.6182_FLOAT,3.1753_FLOAT,3.5210_FLOAT,
     *3.4639_FLOAT,3.4265_FLOAT,3.4075_FLOAT,3.3828_FLOAT,3.3474_FLOAT,
     *3.4071_FLOAT,3.3754_FLOAT,3.3646_FLOAT,3.3308_FLOAT,3.4393_FLOAT,
     *3.2993_FLOAT,3.8768_FLOAT,3.9891_FLOAT,3.8310_FLOAT,3.7483_FLOAT,
     *3.3417_FLOAT,3.3019_FLOAT,3.2250_FLOAT,3.1832_FLOAT,3.1578_FLOAT,
     *3.1564_FLOAT,3.1224_FLOAT,3.4620_FLOAT,2.9743_FLOAT,2.8058_FLOAT,
     *3.4830_FLOAT,3.3474_FLOAT,3.6863_FLOAT,3.3617_FLOAT,3.1608_FLOAT,
     *3.0069_FLOAT,2.9640_FLOAT,2.8427_FLOAT,3.5885_FLOAT,3.5219_FLOAT,
     *4.1314_FLOAT,3.8120_FLOAT/
      DATA T0AB(1191:1260)/3.6015_FLOAT,3.4502_FLOAT,3.3498_FLOAT,
     *3.2777_FLOAT,3.8635_FLOAT,3.8232_FLOAT,3.8486_FLOAT,3.7215_FLOAT,
     *3.6487_FLOAT,3.4724_FLOAT,3.5627_FLOAT,3.5087_FLOAT,3.4757_FLOAT,
     *3.4517_FLOAT,3.4423_FLOAT,3.4139_FLOAT,4.1028_FLOAT,3.8388_FLOAT,
     *3.6745_FLOAT,3.5562_FLOAT,3.4806_FLOAT,3.4272_FLOAT,4.0182_FLOAT,
     *3.9991_FLOAT,4.0007_FLOAT,3.9282_FLOAT,3.7238_FLOAT,3.6498_FLOAT,
     *3.5605_FLOAT,3.5211_FLOAT,3.5009_FLOAT,3.4859_FLOAT,3.4785_FLOAT,
     *3.5621_FLOAT,4.2623_FLOAT,3.0775_FLOAT,2.8275_FLOAT,4.0181_FLOAT,
     *3.3385_FLOAT,3.5379_FLOAT,3.5036_FLOAT,3.2589_FLOAT,3.0804_FLOAT,
     *3.0094_FLOAT,2.9003_FLOAT,4.0869_FLOAT,3.5088_FLOAT,3.9105_FLOAT,
     *3.9833_FLOAT,3.7176_FLOAT,3.5323_FLOAT,3.4102_FLOAT,3.3227_FLOAT,
     *4.2702_FLOAT,4.0888_FLOAT,3.7560_FLOAT,3.7687_FLOAT,3.6681_FLOAT,
     *3.6405_FLOAT,3.5569_FLOAT,3.4990_FLOAT,3.4659_FLOAT,3.4433_FLOAT,
     *3.4330_FLOAT,3.4092_FLOAT,3.8867_FLOAT,4.0190_FLOAT,3.7961_FLOAT,
     *3.6412_FLOAT,3.5405_FLOAT/
      DATA T0AB(1261:1330)/3.4681_FLOAT,4.3538_FLOAT,4.2136_FLOAT,
     *3.9381_FLOAT,3.8912_FLOAT,3.9681_FLOAT,3.7909_FLOAT,3.6774_FLOAT,
     *3.6262_FLOAT,3.5999_FLOAT,3.5823_FLOAT,3.5727_FLOAT,3.5419_FLOAT,
     *4.0245_FLOAT,4.1874_FLOAT,3.0893_FLOAT,2.7917_FLOAT,3.7262_FLOAT,
     *3.3518_FLOAT,3.4241_FLOAT,3.5433_FLOAT,3.2773_FLOAT,3.0890_FLOAT,
     *2.9775_FLOAT,2.9010_FLOAT,3.8048_FLOAT,3.5362_FLOAT,3.7746_FLOAT,
     *3.7911_FLOAT,3.7511_FLOAT,3.5495_FLOAT,3.4149_FLOAT,3.3177_FLOAT,
     *4.0129_FLOAT,3.8370_FLOAT,3.7739_FLOAT,3.7125_FLOAT,3.7152_FLOAT,
     *3.7701_FLOAT,3.5813_FLOAT,3.5187_FLOAT,3.4835_FLOAT,3.4595_FLOAT,
     *3.4439_FLOAT,3.4242_FLOAT,3.7476_FLOAT,3.8239_FLOAT,3.8346_FLOAT,
     *3.6627_FLOAT,3.5479_FLOAT,3.4639_FLOAT,4.1026_FLOAT,3.9733_FLOAT,
     *3.9292_FLOAT,3.8667_FLOAT,3.9513_FLOAT,3.8959_FLOAT,3.7698_FLOAT,
     *3.7089_FLOAT,3.6765_FLOAT,3.6548_FLOAT,3.6409_FLOAT,3.5398_FLOAT,
     *3.8759_FLOAT,3.9804_FLOAT,4.0150_FLOAT,2.9091_FLOAT,2.7638_FLOAT,
     *3.5066_FLOAT,3.3377_FLOAT/
      DATA T0AB(1331:1400)/3.3481_FLOAT,3.2633_FLOAT,3.1810_FLOAT,
     *3.1428_FLOAT,2.9872_FLOAT,2.8837_FLOAT,3.5929_FLOAT,3.5183_FLOAT,
     *3.6729_FLOAT,3.6596_FLOAT,3.6082_FLOAT,3.5927_FLOAT,3.4224_FLOAT,
     *3.2997_FLOAT,3.8190_FLOAT,4.1865_FLOAT,4.1114_FLOAT,4.0540_FLOAT,
     *3.6325_FLOAT,3.5697_FLOAT,3.5561_FLOAT,3.5259_FLOAT,3.4901_FLOAT,
     *3.4552_FLOAT,3.4315_FLOAT,3.4091_FLOAT,3.6438_FLOAT,3.6879_FLOAT,
     *3.6832_FLOAT,3.7043_FLOAT,3.5557_FLOAT,3.4466_FLOAT,3.9203_FLOAT,
     *4.2919_FLOAT,4.2196_FLOAT,4.1542_FLOAT,3.7573_FLOAT,3.7039_FLOAT,
     *3.6546_FLOAT,3.6151_FLOAT,3.5293_FLOAT,3.4849_FLOAT,3.4552_FLOAT,
     *3.5192_FLOAT,3.7673_FLOAT,3.8359_FLOAT,3.8525_FLOAT,3.8901_FLOAT,
     *2.7806_FLOAT,2.7209_FLOAT,3.3812_FLOAT,3.4958_FLOAT,3.2913_FLOAT,
     *3.1888_FLOAT,3.0990_FLOAT,3.0394_FLOAT,2.9789_FLOAT,2.8582_FLOAT,
     *3.4716_FLOAT,3.6883_FLOAT,3.6105_FLOAT,3.5704_FLOAT,3.5059_FLOAT,
     *3.4619_FLOAT,3.4138_FLOAT,3.2742_FLOAT,3.7080_FLOAT,3.9773_FLOAT,
     *3.9010_FLOAT,3.8409_FLOAT/
      DATA T0AB(1401:1470)/3.7944_FLOAT,3.4465_FLOAT,3.7235_FLOAT,
     *3.6808_FLOAT,3.6453_FLOAT,3.6168_FLOAT,3.5844_FLOAT,3.5576_FLOAT,
     *3.5772_FLOAT,3.5959_FLOAT,3.5768_FLOAT,3.5678_FLOAT,3.5486_FLOAT,
     *3.4228_FLOAT,3.8107_FLOAT,4.0866_FLOAT,4.0169_FLOAT,3.9476_FLOAT,
     *3.6358_FLOAT,3.5800_FLOAT,3.5260_FLOAT,3.4838_FLOAT,3.4501_FLOAT,
     *3.4204_FLOAT,3.3553_FLOAT,3.6487_FLOAT,3.6973_FLOAT,3.7398_FLOAT,
     *3.7405_FLOAT,3.7459_FLOAT,3.7380_FLOAT,2.6848_FLOAT,2.6740_FLOAT,
     *3.2925_FLOAT,3.3386_FLOAT,3.2473_FLOAT,3.1284_FLOAT,3.0301_FLOAT,
     *2.9531_FLOAT,2.9602_FLOAT,2.8272_FLOAT,3.3830_FLOAT,3.5358_FLOAT,
     *3.5672_FLOAT,3.5049_FLOAT,3.4284_FLOAT,3.3621_FLOAT,3.3001_FLOAT,
     *3.2451_FLOAT,3.6209_FLOAT,3.8299_FLOAT,3.7543_FLOAT,3.6920_FLOAT,
     *3.6436_FLOAT,3.3598_FLOAT,3.5701_FLOAT,3.5266_FLOAT,3.4904_FLOAT,
     *3.4590_FLOAT,3.4364_FLOAT,3.4077_FLOAT,3.5287_FLOAT,3.5280_FLOAT,
     *3.4969_FLOAT,3.4650_FLOAT,3.4304_FLOAT,3.3963_FLOAT,3.7229_FLOAT,
     *3.9402_FLOAT,3.8753_FLOAT/
      DATA T0AB(1471:1540)/3.8035_FLOAT,3.5499_FLOAT,3.4913_FLOAT,
     *3.4319_FLOAT,3.3873_FLOAT,3.3520_FLOAT,3.3209_FLOAT,3.2948_FLOAT,
     *3.5052_FLOAT,3.6465_FLOAT,3.6696_FLOAT,3.6577_FLOAT,3.6388_FLOAT,
     *3.6142_FLOAT,3.5889_FLOAT,3.3968_FLOAT,3.0122_FLOAT,4.2241_FLOAT,
     *3.7887_FLOAT,4.0049_FLOAT,3.5384_FLOAT,3.2698_FLOAT,3.1083_FLOAT,
     *2.9917_FLOAT,2.9057_FLOAT,4.3340_FLOAT,3.9900_FLOAT,4.6588_FLOAT,
     *4.1278_FLOAT,3.8125_FLOAT,3.6189_FLOAT,3.4851_FLOAT,3.3859_FLOAT,
     *4.6531_FLOAT,4.3134_FLOAT,4.2258_FLOAT,4.1309_FLOAT,4.0692_FLOAT,
     *4.0944_FLOAT,3.9850_FLOAT,3.9416_FLOAT,3.9112_FLOAT,3.8873_FLOAT,
     *3.8736_FLOAT,3.8473_FLOAT,4.6027_FLOAT,4.1538_FLOAT,3.8994_FLOAT,
     *3.7419_FLOAT,3.6356_FLOAT,3.5548_FLOAT,4.8353_FLOAT,4.5413_FLOAT,
     *4.3891_FLOAT,4.3416_FLOAT,4.3243_FLOAT,4.2753_FLOAT,4.2053_FLOAT,
     *4.1790_FLOAT,4.1685_FLOAT,4.1585_FLOAT,4.1536_FLOAT,4.0579_FLOAT,
     *4.1980_FLOAT,4.4564_FLOAT,4.2192_FLOAT,4.0528_FLOAT,3.9489_FLOAT,
     *3.8642_FLOAT,5.0567_FLOAT/
      DATA T0AB(1541:1610)/3.0630_FLOAT,3.3271_FLOAT,4.0432_FLOAT,
     *4.0046_FLOAT,4.1555_FLOAT,3.7426_FLOAT,3.5130_FLOAT,3.5174_FLOAT,
     *3.2884_FLOAT,3.1378_FLOAT,4.1894_FLOAT,4.2321_FLOAT,4.1725_FLOAT,
     *4.1833_FLOAT,3.8929_FLOAT,4.0544_FLOAT,3.8118_FLOAT,3.6414_FLOAT,
     *4.6373_FLOAT,4.6268_FLOAT,4.4750_FLOAT,4.4134_FLOAT,4.3458_FLOAT,
     *3.8582_FLOAT,4.2583_FLOAT,4.1898_FLOAT,4.1562_FLOAT,4.1191_FLOAT,
     *4.1069_FLOAT,4.0639_FLOAT,4.1257_FLOAT,4.1974_FLOAT,3.9532_FLOAT,
     *4.1794_FLOAT,3.9660_FLOAT,3.8130_FLOAT,4.8160_FLOAT,4.8272_FLOAT,
     *4.6294_FLOAT,4.5840_FLOAT,4.0770_FLOAT,4.0088_FLOAT,3.9103_FLOAT,
     *3.8536_FLOAT,3.8324_FLOAT,3.7995_FLOAT,3.7826_FLOAT,4.2294_FLOAT,
     *4.3380_FLOAT,4.4352_FLOAT,4.1933_FLOAT,4.4580_FLOAT,4.2554_FLOAT,
     *4.1072_FLOAT,5.0454_FLOAT,5.1814_FLOAT,3.0632_FLOAT,3.2662_FLOAT,
     *3.6432_FLOAT,3.8088_FLOAT,3.7910_FLOAT,3.7381_FLOAT,3.5093_FLOAT,
     *3.5155_FLOAT,3.3047_FLOAT,3.1681_FLOAT,3.7871_FLOAT,3.9924_FLOAT,
     *4.0637_FLOAT,4.1382_FLOAT/
      DATA T0AB(1611:1680)/3.8591_FLOAT,4.0164_FLOAT,3.7878_FLOAT,
     *3.6316_FLOAT,4.1741_FLOAT,4.3166_FLOAT,4.2395_FLOAT,4.1831_FLOAT,
     *4.1107_FLOAT,3.5857_FLOAT,4.0270_FLOAT,3.9676_FLOAT,3.9463_FLOAT,
     *3.9150_FLOAT,3.9021_FLOAT,3.8708_FLOAT,4.0240_FLOAT,4.1551_FLOAT,
     *3.9108_FLOAT,4.1337_FLOAT,3.9289_FLOAT,3.7873_FLOAT,4.3666_FLOAT,
     *4.5080_FLOAT,4.4232_FLOAT,4.3155_FLOAT,3.8461_FLOAT,3.8007_FLOAT,
     *3.6991_FLOAT,3.6447_FLOAT,3.6308_FLOAT,3.5959_FLOAT,3.5749_FLOAT,
     *4.0359_FLOAT,4.3124_FLOAT,4.3539_FLOAT,4.1122_FLOAT,4.3772_FLOAT,
     *4.1785_FLOAT,4.0386_FLOAT,4.7004_FLOAT,4.8604_FLOAT,4.6261_FLOAT,
     *2.9455_FLOAT,3.2470_FLOAT,3.6108_FLOAT,3.8522_FLOAT,3.6625_FLOAT,
     *3.6598_FLOAT,3.4411_FLOAT,3.4660_FLOAT,3.2415_FLOAT,3.0944_FLOAT,
     *3.7514_FLOAT,4.0397_FLOAT,3.9231_FLOAT,4.0561_FLOAT,3.7860_FLOAT,
     *3.9845_FLOAT,3.7454_FLOAT,3.5802_FLOAT,4.1366_FLOAT,4.3581_FLOAT,
     *4.2351_FLOAT,4.2011_FLOAT,4.1402_FLOAT,3.5381_FLOAT,4.0653_FLOAT,
     *4.0093_FLOAT,3.9883_FLOAT/
      DATA T0AB(1681:1750)/3.9570_FLOAT,3.9429_FLOAT,3.9112_FLOAT,
     *3.8728_FLOAT,4.0682_FLOAT,3.8351_FLOAT,4.1054_FLOAT,3.8928_FLOAT,
     *3.7445_FLOAT,4.3415_FLOAT,4.5497_FLOAT,4.3833_FLOAT,4.3122_FLOAT,
     *3.8051_FLOAT,3.7583_FLOAT,3.6622_FLOAT,3.6108_FLOAT,3.5971_FLOAT,
     *3.5628_FLOAT,3.5408_FLOAT,4.0780_FLOAT,4.0727_FLOAT,4.2836_FLOAT,
     *4.0553_FLOAT,4.3647_FLOAT,4.1622_FLOAT,4.0178_FLOAT,4.5802_FLOAT,
     *4.9125_FLOAT,4.5861_FLOAT,4.6201_FLOAT,2.9244_FLOAT,3.2241_FLOAT,
     *3.5848_FLOAT,3.8293_FLOAT,3.6395_FLOAT,3.6400_FLOAT,3.4204_FLOAT,
     *3.4499_FLOAT,3.2253_FLOAT,3.0779_FLOAT,3.7257_FLOAT,4.0170_FLOAT,
     *3.9003_FLOAT,4.0372_FLOAT,3.7653_FLOAT,3.9672_FLOAT,3.7283_FLOAT,
     *3.5630_FLOAT,4.1092_FLOAT,4.3347_FLOAT,4.2117_FLOAT,4.1793_FLOAT,
     *4.1179_FLOAT,3.5139_FLOAT,4.0426_FLOAT,3.9867_FLOAT,3.9661_FLOAT,
     *3.9345_FLOAT,3.9200_FLOAT,3.8883_FLOAT,3.8498_FLOAT,4.0496_FLOAT,
     *3.8145_FLOAT,4.0881_FLOAT,3.8756_FLOAT,3.7271_FLOAT,4.3128_FLOAT,
     *4.5242_FLOAT,4.3578_FLOAT/
      DATA T0AB(1751:1820)/4.2870_FLOAT,3.7796_FLOAT,3.7318_FLOAT,
     *3.6364_FLOAT,3.5854_FLOAT,3.5726_FLOAT,3.5378_FLOAT,3.5155_FLOAT,
     *4.0527_FLOAT,4.0478_FLOAT,4.2630_FLOAT,4.0322_FLOAT,4.3449_FLOAT,
     *4.1421_FLOAT,3.9975_FLOAT,4.5499_FLOAT,4.8825_FLOAT,4.5601_FLOAT,
     *4.5950_FLOAT,4.5702_FLOAT,2.9046_FLOAT,3.2044_FLOAT,3.5621_FLOAT,
     *3.8078_FLOAT,3.6185_FLOAT,3.6220_FLOAT,3.4019_FLOAT,3.4359_FLOAT,
     *3.2110_FLOAT,3.0635_FLOAT,3.7037_FLOAT,3.9958_FLOAT,3.8792_FLOAT,
     *4.0194_FLOAT,3.7460_FLOAT,3.9517_FLOAT,3.7128_FLOAT,3.5474_FLOAT,
     *4.0872_FLOAT,4.3138_FLOAT,4.1906_FLOAT,4.1593_FLOAT,4.0973_FLOAT,
     *3.4919_FLOAT,4.0216_FLOAT,3.9657_FLOAT,3.9454_FLOAT,3.9134_FLOAT,
     *3.8986_FLOAT,3.8669_FLOAT,3.8289_FLOAT,4.0323_FLOAT,3.7954_FLOAT,
     *4.0725_FLOAT,3.8598_FLOAT,3.7113_FLOAT,4.2896_FLOAT,4.5021_FLOAT,
     *4.3325_FLOAT,4.2645_FLOAT,3.7571_FLOAT,3.7083_FLOAT,3.6136_FLOAT,
     *3.5628_FLOAT,3.5507_FLOAT,3.5155_FLOAT,3.4929_FLOAT,4.0297_FLOAT,
     *4.0234_FLOAT,4.2442_FLOAT/
      DATA T0AB(1821:1890)/4.0112_FLOAT,4.3274_FLOAT,4.1240_FLOAT,
     *3.9793_FLOAT,4.5257_FLOAT,4.8568_FLOAT,4.5353_FLOAT,4.5733_FLOAT,
     *4.5485_FLOAT,4.5271_FLOAT,2.8878_FLOAT,3.1890_FLOAT,3.5412_FLOAT,
     *3.7908_FLOAT,3.5974_FLOAT,3.6078_FLOAT,3.3871_FLOAT,3.4243_FLOAT,
     *3.1992_FLOAT,3.0513_FLOAT,3.6831_FLOAT,3.9784_FLOAT,3.8579_FLOAT,
     *4.0049_FLOAT,3.7304_FLOAT,3.9392_FLOAT,3.7002_FLOAT,3.5347_FLOAT,
     *4.0657_FLOAT,4.2955_FLOAT,4.1705_FLOAT,4.1424_FLOAT,4.0800_FLOAT,
     *3.4717_FLOAT,4.0043_FLOAT,3.9485_FLOAT,3.9286_FLOAT,3.8965_FLOAT,
     *3.8815_FLOAT,3.8500_FLOAT,3.8073_FLOAT,4.0180_FLOAT,3.7796_FLOAT,
     *4.0598_FLOAT,3.8470_FLOAT,3.6983_FLOAT,4.2678_FLOAT,4.4830_FLOAT,
     *4.3132_FLOAT,4.2444_FLOAT,3.7370_FLOAT,3.6876_FLOAT,3.5935_FLOAT,
     *3.5428_FLOAT,3.5314_FLOAT,3.4958_FLOAT,3.4730_FLOAT,4.0117_FLOAT,
     *4.0043_FLOAT,4.2287_FLOAT,3.9939_FLOAT,4.3134_FLOAT,4.1096_FLOAT,
     *3.9646_FLOAT,4.5032_FLOAT,4.8356_FLOAT,4.5156_FLOAT,4.5544_FLOAT,
     *4.5297_FLOAT,4.5083_FLOAT/
      DATA T0AB(1891:1960)/4.4896_FLOAT,2.8709_FLOAT,3.1737_FLOAT,
     *3.5199_FLOAT,3.7734_FLOAT,3.5802_FLOAT,3.5934_FLOAT,3.3724_FLOAT,
     *3.4128_FLOAT,3.1877_FLOAT,3.0396_FLOAT,3.6624_FLOAT,3.9608_FLOAT,
     *3.8397_FLOAT,3.9893_FLOAT,3.7145_FLOAT,3.9266_FLOAT,3.6877_FLOAT,
     *3.5222_FLOAT,4.0448_FLOAT,4.2771_FLOAT,4.1523_FLOAT,4.1247_FLOAT,
     *4.0626_FLOAT,3.4530_FLOAT,3.9866_FLOAT,3.9310_FLOAT,3.9115_FLOAT,
     *3.8792_FLOAT,3.8641_FLOAT,3.8326_FLOAT,3.7892_FLOAT,4.0025_FLOAT,
     *3.7636_FLOAT,4.0471_FLOAT,3.8343_FLOAT,3.6854_FLOAT,4.2464_FLOAT,
     *4.4635_FLOAT,4.2939_FLOAT,4.2252_FLOAT,3.7169_FLOAT,3.6675_FLOAT,
     *3.5739_FLOAT,3.5235_FLOAT,3.5126_FLOAT,3.4768_FLOAT,3.4537_FLOAT,
     *3.9932_FLOAT,3.9854_FLOAT,4.2123_FLOAT,3.9765_FLOAT,4.2992_FLOAT,
     *4.0951_FLOAT,3.9500_FLOAT,4.4811_FLOAT,4.8135_FLOAT,4.4959_FLOAT,
     *4.5351_FLOAT,4.5105_FLOAT,4.4891_FLOAT,4.4705_FLOAT,4.4515_FLOAT,
     *2.8568_FLOAT,3.1608_FLOAT,3.5050_FLOAT,3.7598_FLOAT,3.5665_FLOAT,
     *3.5803_FLOAT,3.3601_FLOAT/
      DATA T0AB(1961:2030)/3.4031_FLOAT,3.1779_FLOAT,3.0296_FLOAT,
     *3.6479_FLOAT,3.9471_FLOAT,3.8262_FLOAT,3.9773_FLOAT,3.7015_FLOAT,
     *3.9162_FLOAT,3.6771_FLOAT,3.5115_FLOAT,4.0306_FLOAT,4.2634_FLOAT,
     *4.1385_FLOAT,4.1116_FLOAT,4.0489_FLOAT,3.4366_FLOAT,3.9732_FLOAT,
     *3.9176_FLOAT,3.8983_FLOAT,3.8659_FLOAT,3.8507_FLOAT,3.8191_FLOAT,
     *3.7757_FLOAT,3.9907_FLOAT,3.7506_FLOAT,4.0365_FLOAT,3.8235_FLOAT,
     *3.6745_FLOAT,4.2314_FLOAT,4.4490_FLOAT,4.2792_FLOAT,4.2105_FLOAT,
     *3.7003_FLOAT,3.6510_FLOAT,3.5578_FLOAT,3.5075_FLOAT,3.4971_FLOAT,
     *3.4609_FLOAT,3.4377_FLOAT,3.9788_FLOAT,3.9712_FLOAT,4.1997_FLOAT,
     *3.9624_FLOAT,4.2877_FLOAT,4.0831_FLOAT,3.9378_FLOAT,4.4655_FLOAT,
     *4.7974_FLOAT,4.4813_FLOAT,4.5209_FLOAT,4.4964_FLOAT,4.4750_FLOAT,
     *4.4565_FLOAT,4.4375_FLOAT,4.4234_FLOAT,2.6798_FLOAT,3.0151_FLOAT,
     *3.2586_FLOAT,3.5292_FLOAT,3.5391_FLOAT,3.4902_FLOAT,3.2887_FLOAT,
     *3.3322_FLOAT,3.1228_FLOAT,2.9888_FLOAT,3.4012_FLOAT,3.7145_FLOAT,
     *3.7830_FLOAT,3.6665_FLOAT/
      DATA T0AB(2031:2100)/3.5898_FLOAT,3.8077_FLOAT,3.5810_FLOAT,
     *3.4265_FLOAT,3.7726_FLOAT,4.0307_FLOAT,3.9763_FLOAT,3.8890_FLOAT,
     *3.8489_FLOAT,3.2706_FLOAT,3.7595_FLOAT,3.6984_FLOAT,3.6772_FLOAT,
     *3.6428_FLOAT,3.6243_FLOAT,3.5951_FLOAT,3.7497_FLOAT,3.6775_FLOAT,
     *3.6364_FLOAT,3.9203_FLOAT,3.7157_FLOAT,3.5746_FLOAT,3.9494_FLOAT,
     *4.2076_FLOAT,4.1563_FLOAT,4.0508_FLOAT,3.5329_FLOAT,3.4780_FLOAT,
     *3.3731_FLOAT,3.3126_FLOAT,3.2846_FLOAT,3.2426_FLOAT,3.2135_FLOAT,
     *3.7491_FLOAT,3.9006_FLOAT,3.8332_FLOAT,3.8029_FLOAT,4.1436_FLOAT,
     *3.9407_FLOAT,3.7998_FLOAT,4.1663_FLOAT,4.5309_FLOAT,4.3481_FLOAT,
     *4.2911_FLOAT,4.2671_FLOAT,4.2415_FLOAT,4.2230_FLOAT,4.2047_FLOAT,
     *4.1908_FLOAT,4.1243_FLOAT,2.5189_FLOAT,2.9703_FLOAT,3.3063_FLOAT,
     *3.6235_FLOAT,3.4517_FLOAT,3.3989_FLOAT,3.2107_FLOAT,3.2434_FLOAT,
     *3.0094_FLOAT,2.8580_FLOAT,3.4253_FLOAT,3.8157_FLOAT,3.7258_FLOAT,
     *3.6132_FLOAT,3.5297_FLOAT,3.7566_FLOAT,3.5095_FLOAT,3.3368_FLOAT,
     *3.7890_FLOAT,4.1298_FLOAT/
      DATA T0AB(2101:2170)/4.0190_FLOAT,3.9573_FLOAT,3.9237_FLOAT,
     *3.2677_FLOAT,3.8480_FLOAT,3.8157_FLOAT,3.7656_FLOAT,3.7317_FLOAT,
     *3.7126_FLOAT,3.6814_FLOAT,3.6793_FLOAT,3.6218_FLOAT,3.5788_FLOAT,
     *3.8763_FLOAT,3.6572_FLOAT,3.5022_FLOAT,3.9737_FLOAT,4.3255_FLOAT,
     *4.1828_FLOAT,4.1158_FLOAT,3.5078_FLOAT,3.4595_FLOAT,3.3600_FLOAT,
     *3.3088_FLOAT,3.2575_FLOAT,3.2164_FLOAT,3.1856_FLOAT,3.8522_FLOAT,
     *3.8665_FLOAT,3.8075_FLOAT,3.7772_FLOAT,4.1391_FLOAT,3.9296_FLOAT,
     *3.7772_FLOAT,4.2134_FLOAT,4.7308_FLOAT,4.3787_FLOAT,4.3894_FLOAT,
     *4.3649_FLOAT,4.3441_FLOAT,4.3257_FLOAT,4.3073_FLOAT,4.2941_FLOAT,
     *4.1252_FLOAT,4.2427_FLOAT,3.0481_FLOAT,2.9584_FLOAT,3.6919_FLOAT,
     *3.5990_FLOAT,3.8881_FLOAT,3.4209_FLOAT,3.1606_FLOAT,3.1938_FLOAT,
     *2.9975_FLOAT,2.8646_FLOAT,3.8138_FLOAT,3.7935_FLOAT,3.7081_FLOAT,
     *3.9155_FLOAT,3.5910_FLOAT,3.4808_FLOAT,3.4886_FLOAT,3.3397_FLOAT,
     *4.1336_FLOAT,4.1122_FLOAT,3.9888_FLOAT,3.9543_FLOAT,3.8917_FLOAT,
     *3.5894_FLOAT,3.8131_FLOAT/
      DATA T0AB(2171:2240)/3.7635_FLOAT,3.7419_FLOAT,3.7071_FLOAT,
     *3.6880_FLOAT,3.6574_FLOAT,3.6546_FLOAT,3.9375_FLOAT,3.6579_FLOAT,
     *3.5870_FLOAT,3.6361_FLOAT,3.5039_FLOAT,4.3149_FLOAT,4.2978_FLOAT,
     *4.1321_FLOAT,4.1298_FLOAT,3.8164_FLOAT,3.7680_FLOAT,3.7154_FLOAT,
     *3.6858_FLOAT,3.6709_FLOAT,3.6666_FLOAT,3.6517_FLOAT,3.8174_FLOAT,
     *3.8608_FLOAT,4.1805_FLOAT,3.9102_FLOAT,3.8394_FLOAT,3.8968_FLOAT,
     *3.7673_FLOAT,4.5274_FLOAT,4.6682_FLOAT,4.3344_FLOAT,4.3639_FLOAT,
     *4.3384_FLOAT,4.3162_FLOAT,4.2972_FLOAT,4.2779_FLOAT,4.2636_FLOAT,
     *4.0253_FLOAT,4.1168_FLOAT,4.1541_FLOAT,2.8136_FLOAT,3.0951_FLOAT,
     *3.4635_FLOAT,3.6875_FLOAT,3.4987_FLOAT,3.5183_FLOAT,3.2937_FLOAT,
     *3.3580_FLOAT,3.1325_FLOAT,2.9832_FLOAT,3.6078_FLOAT,3.8757_FLOAT,
     *3.7616_FLOAT,3.9222_FLOAT,3.6370_FLOAT,3.8647_FLOAT,3.6256_FLOAT,
     *3.4595_FLOAT,3.9874_FLOAT,4.1938_FLOAT,4.0679_FLOAT,4.0430_FLOAT,
     *3.9781_FLOAT,3.3886_FLOAT,3.9008_FLOAT,3.8463_FLOAT,3.8288_FLOAT,
     *3.7950_FLOAT,3.7790_FLOAT/
      DATA T0AB(2241:2310)/3.7472_FLOAT,3.7117_FLOAT,3.9371_FLOAT,
     *3.6873_FLOAT,3.9846_FLOAT,3.7709_FLOAT,3.6210_FLOAT,4.1812_FLOAT,
     *4.3750_FLOAT,4.2044_FLOAT,4.1340_FLOAT,3.6459_FLOAT,3.5929_FLOAT,
     *3.5036_FLOAT,3.4577_FLOAT,3.4528_FLOAT,3.4146_FLOAT,3.3904_FLOAT,
     *3.9014_FLOAT,3.9031_FLOAT,4.1443_FLOAT,3.8961_FLOAT,4.2295_FLOAT,
     *4.0227_FLOAT,3.8763_FLOAT,4.4086_FLOAT,4.7097_FLOAT,4.4064_FLOAT,
     *4.4488_FLOAT,4.4243_FLOAT,4.4029_FLOAT,4.3842_FLOAT,4.3655_FLOAT,
     *4.3514_FLOAT,4.1162_FLOAT,4.2205_FLOAT,4.1953_FLOAT,4.2794_FLOAT,
     *2.8032_FLOAT,3.0805_FLOAT,3.4519_FLOAT,3.6700_FLOAT,3.4827_FLOAT,
     *3.5050_FLOAT,3.2799_FLOAT,3.3482_FLOAT,3.1233_FLOAT,2.9747_FLOAT,
     *3.5971_FLOAT,3.8586_FLOAT,3.7461_FLOAT,3.9100_FLOAT,3.6228_FLOAT,
     *3.8535_FLOAT,3.6147_FLOAT,3.4490_FLOAT,3.9764_FLOAT,4.1773_FLOAT,
     *4.0511_FLOAT,4.0270_FLOAT,3.9614_FLOAT,3.3754_FLOAT,3.8836_FLOAT,
     *3.8291_FLOAT,3.8121_FLOAT,3.7780_FLOAT,3.7619_FLOAT,3.7300_FLOAT,
     *3.6965_FLOAT,3.9253_FLOAT/
      DATA T0AB(2311:2380)/3.6734_FLOAT,3.9733_FLOAT,3.7597_FLOAT,
     *3.6099_FLOAT,4.1683_FLOAT,4.3572_FLOAT,4.1862_FLOAT,4.1153_FLOAT,
     *3.6312_FLOAT,3.5772_FLOAT,3.4881_FLOAT,3.4429_FLOAT,3.4395_FLOAT,
     *3.4009_FLOAT,3.3766_FLOAT,3.8827_FLOAT,3.8868_FLOAT,4.1316_FLOAT,
     *3.8807_FLOAT,4.2164_FLOAT,4.0092_FLOAT,3.8627_FLOAT,4.3936_FLOAT,
     *4.6871_FLOAT,4.3882_FLOAT,4.4316_FLOAT,4.4073_FLOAT,4.3858_FLOAT,
     *4.3672_FLOAT,4.3485_FLOAT,4.3344_FLOAT,4.0984_FLOAT,4.2036_FLOAT,
     *4.1791_FLOAT,4.2622_FLOAT,4.2450_FLOAT,2.7967_FLOAT,3.0689_FLOAT,
     *3.4445_FLOAT,3.6581_FLOAT,3.4717_FLOAT,3.4951_FLOAT,3.2694_FLOAT,
     *3.3397_FLOAT,3.1147_FLOAT,2.9661_FLOAT,3.5898_FLOAT,3.8468_FLOAT,
     *3.7358_FLOAT,3.9014_FLOAT,3.6129_FLOAT,3.8443_FLOAT,3.6054_FLOAT,
     *3.4396_FLOAT,3.9683_FLOAT,4.1656_FLOAT,4.0394_FLOAT,4.0158_FLOAT,
     *3.9498_FLOAT,3.3677_FLOAT,3.8718_FLOAT,3.8164_FLOAT,3.8005_FLOAT,
     *3.7662_FLOAT,3.7500_FLOAT,3.7181_FLOAT,3.6863_FLOAT,3.9170_FLOAT,
     *3.6637_FLOAT,3.9641_FLOAT/
      DATA T0AB(2381:2450)/3.7503_FLOAT,3.6004_FLOAT,4.1590_FLOAT,
     *4.3448_FLOAT,4.1739_FLOAT,4.1029_FLOAT,3.6224_FLOAT,3.5677_FLOAT,
     *3.4785_FLOAT,3.4314_FLOAT,3.4313_FLOAT,3.3923_FLOAT,3.3680_FLOAT,
     *3.8698_FLOAT,3.8758_FLOAT,4.1229_FLOAT,3.8704_FLOAT,4.2063_FLOAT,
     *3.9987_FLOAT,3.8519_FLOAT,4.3832_FLOAT,4.6728_FLOAT,4.3759_FLOAT,
     *4.4195_FLOAT,4.3952_FLOAT,4.3737_FLOAT,4.3551_FLOAT,4.3364_FLOAT,
     *4.3223_FLOAT,4.0861_FLOAT,4.1911_FLOAT,4.1676_FLOAT,4.2501_FLOAT,
     *4.2329_FLOAT,4.2208_FLOAT,2.7897_FLOAT,3.0636_FLOAT,3.4344_FLOAT,
     *3.6480_FLOAT,3.4626_FLOAT,3.4892_FLOAT,3.2626_FLOAT,3.3344_FLOAT,
     *3.1088_FLOAT,2.9597_FLOAT,3.5804_FLOAT,3.8359_FLOAT,3.7251_FLOAT,
     *3.8940_FLOAT,3.6047_FLOAT,3.8375_FLOAT,3.5990_FLOAT,3.4329_FLOAT,
     *3.9597_FLOAT,4.1542_FLOAT,4.0278_FLOAT,4.0048_FLOAT,3.9390_FLOAT,
     *3.3571_FLOAT,3.8608_FLOAT,3.8056_FLOAT,3.7899_FLOAT,3.7560_FLOAT,
     *3.7400_FLOAT,3.7081_FLOAT,3.6758_FLOAT,3.9095_FLOAT,3.6552_FLOAT,
     *3.9572_FLOAT,3.7436_FLOAT/
      DATA T0AB(2451:2520)/3.5933_FLOAT,4.1508_FLOAT,4.3337_FLOAT,
     *4.1624_FLOAT,4.0916_FLOAT,3.6126_FLOAT,3.5582_FLOAT,3.4684_FLOAT,
     *3.4212_FLOAT,3.4207_FLOAT,3.3829_FLOAT,3.3586_FLOAT,3.8604_FLOAT,
     *3.8658_FLOAT,4.1156_FLOAT,3.8620_FLOAT,4.1994_FLOAT,3.9917_FLOAT,
     *3.8446_FLOAT,4.3750_FLOAT,4.6617_FLOAT,4.3644_FLOAT,4.4083_FLOAT,
     *4.3840_FLOAT,4.3625_FLOAT,4.3439_FLOAT,4.3253_FLOAT,4.3112_FLOAT,
     *4.0745_FLOAT,4.1807_FLOAT,4.1578_FLOAT,4.2390_FLOAT,4.2218_FLOAT,
     *4.2097_FLOAT,4.1986_FLOAT,2.8395_FLOAT,3.0081_FLOAT,3.3171_FLOAT,
     *3.4878_FLOAT,3.5360_FLOAT,3.5145_FLOAT,3.2809_FLOAT,3.3307_FLOAT,
     *3.1260_FLOAT,2.9940_FLOAT,3.4741_FLOAT,3.6675_FLOAT,3.7832_FLOAT,
     *3.6787_FLOAT,3.6156_FLOAT,3.8041_FLOAT,3.5813_FLOAT,3.4301_FLOAT,
     *3.8480_FLOAT,3.9849_FLOAT,3.9314_FLOAT,3.8405_FLOAT,3.8029_FLOAT,
     *3.2962_FLOAT,3.7104_FLOAT,3.6515_FLOAT,3.6378_FLOAT,3.6020_FLOAT,
     *3.5849_FLOAT,3.5550_FLOAT,3.7494_FLOAT,3.6893_FLOAT,3.6666_FLOAT,
     *3.9170_FLOAT,3.7150_FLOAT/
      DATA T0AB(2521:2590)/3.5760_FLOAT,4.0268_FLOAT,4.1596_FLOAT,
     *4.1107_FLOAT,3.9995_FLOAT,3.5574_FLOAT,3.5103_FLOAT,3.4163_FLOAT,
     *3.3655_FLOAT,3.3677_FLOAT,3.3243_FLOAT,3.2975_FLOAT,3.7071_FLOAT,
     *3.9047_FLOAT,3.8514_FLOAT,3.8422_FLOAT,3.8022_FLOAT,3.9323_FLOAT,
     *3.7932_FLOAT,4.2343_FLOAT,4.4583_FLOAT,4.3115_FLOAT,4.2457_FLOAT,
     *4.2213_FLOAT,4.1945_FLOAT,4.1756_FLOAT,4.1569_FLOAT,4.1424_FLOAT,
     *4.0620_FLOAT,4.0494_FLOAT,3.9953_FLOAT,4.0694_FLOAT,4.0516_FLOAT,
     *4.0396_FLOAT,4.0280_FLOAT,4.0130_FLOAT,2.9007_FLOAT,2.9674_FLOAT,
     *3.8174_FLOAT,3.5856_FLOAT,3.6486_FLOAT,3.5339_FLOAT,3.2832_FLOAT,
     *3.3154_FLOAT,3.1144_FLOAT,2.9866_FLOAT,3.9618_FLOAT,3.8430_FLOAT,
     *3.9980_FLOAT,3.8134_FLOAT,3.6652_FLOAT,3.7985_FLOAT,3.5756_FLOAT,
     *3.4207_FLOAT,4.4061_FLOAT,4.2817_FLOAT,4.1477_FLOAT,4.0616_FLOAT,
     *3.9979_FLOAT,3.6492_FLOAT,3.8833_FLOAT,3.8027_FLOAT,3.7660_FLOAT,
     *3.7183_FLOAT,3.6954_FLOAT,3.6525_FLOAT,3.9669_FLOAT,3.8371_FLOAT,
     *3.7325_FLOAT,3.9160_FLOAT/
      DATA T0AB(2591:2660)/3.7156_FLOAT,3.5714_FLOAT,4.6036_FLOAT,
     *4.4620_FLOAT,4.3092_FLOAT,4.2122_FLOAT,3.8478_FLOAT,3.7572_FLOAT,
     *3.6597_FLOAT,3.5969_FLOAT,3.5575_FLOAT,3.5386_FLOAT,3.5153_FLOAT,
     *3.7818_FLOAT,4.1335_FLOAT,4.0153_FLOAT,3.9177_FLOAT,3.8603_FLOAT,
     *3.9365_FLOAT,3.7906_FLOAT,4.7936_FLOAT,4.7410_FLOAT,4.5461_FLOAT,
     *4.5662_FLOAT,4.5340_FLOAT,4.5059_FLOAT,4.4832_FLOAT,4.4604_FLOAT,
     *4.4429_FLOAT,4.2346_FLOAT,4.4204_FLOAT,4.3119_FLOAT,4.3450_FLOAT,
     *4.3193_FLOAT,4.3035_FLOAT,4.2933_FLOAT,4.1582_FLOAT,4.2450_FLOAT,
     *2.8559_FLOAT,2.9050_FLOAT,3.8325_FLOAT,3.5442_FLOAT,3.5077_FLOAT,
     *3.4905_FLOAT,3.2396_FLOAT,3.2720_FLOAT,3.0726_FLOAT,2.9467_FLOAT,
     *3.9644_FLOAT,3.8050_FLOAT,3.8981_FLOAT,3.7762_FLOAT,3.6216_FLOAT,
     *3.7531_FLOAT,3.5297_FLOAT,3.3742_FLOAT,4.3814_FLOAT,4.2818_FLOAT,
     *4.1026_FLOAT,4.0294_FLOAT,3.9640_FLOAT,3.6208_FLOAT,3.8464_FLOAT,
     *3.7648_FLOAT,3.7281_FLOAT,3.6790_FLOAT,3.6542_FLOAT,3.6117_FLOAT,
     *3.8650_FLOAT,3.8010_FLOAT/
      DATA T0AB(2661:2730)/3.6894_FLOAT,3.8713_FLOAT,3.6699_FLOAT,
     *3.5244_FLOAT,4.5151_FLOAT,4.4517_FLOAT,4.2538_FLOAT,4.1483_FLOAT,
     *3.8641_FLOAT,3.7244_FLOAT,3.6243_FLOAT,3.5589_FLOAT,3.5172_FLOAT,
     *3.4973_FLOAT,3.4715_FLOAT,3.7340_FLOAT,4.0316_FLOAT,3.9958_FLOAT,
     *3.8687_FLOAT,3.8115_FLOAT,3.8862_FLOAT,3.7379_FLOAT,4.7091_FLOAT,
     *4.7156_FLOAT,4.5199_FLOAT,4.5542_FLOAT,4.5230_FLOAT,4.4959_FLOAT,
     *4.4750_FLOAT,4.4529_FLOAT,4.4361_FLOAT,4.1774_FLOAT,4.3774_FLOAT,
     *4.2963_FLOAT,4.3406_FLOAT,4.3159_FLOAT,4.3006_FLOAT,4.2910_FLOAT,
     *4.1008_FLOAT,4.1568_FLOAT,4.0980_FLOAT,2.8110_FLOAT,2.8520_FLOAT,
     *3.7480_FLOAT,3.5105_FLOAT,3.4346_FLOAT,3.3461_FLOAT,3.1971_FLOAT,
     *3.2326_FLOAT,3.0329_FLOAT,2.9070_FLOAT,3.8823_FLOAT,3.7928_FLOAT,
     *3.8264_FLOAT,3.7006_FLOAT,3.5797_FLOAT,3.7141_FLOAT,3.4894_FLOAT,
     *3.3326_FLOAT,4.3048_FLOAT,4.2217_FLOAT,4.0786_FLOAT,3.9900_FLOAT,
     *3.9357_FLOAT,3.6331_FLOAT,3.8333_FLOAT,3.7317_FLOAT,3.6957_FLOAT,
     *3.6460_FLOAT,3.6197_FLOAT/
      DATA T0AB(2731:2800)/3.5779_FLOAT,3.7909_FLOAT,3.7257_FLOAT,
     *3.6476_FLOAT,3.5729_FLOAT,3.6304_FLOAT,3.4834_FLOAT,4.4368_FLOAT,
     *4.3921_FLOAT,4.2207_FLOAT,4.1133_FLOAT,3.8067_FLOAT,3.7421_FLOAT,
     *3.6140_FLOAT,3.5491_FLOAT,3.5077_FLOAT,3.4887_FLOAT,3.4623_FLOAT,
     *3.6956_FLOAT,3.9568_FLOAT,3.8976_FLOAT,3.8240_FLOAT,3.7684_FLOAT,
     *3.8451_FLOAT,3.6949_FLOAT,4.6318_FLOAT,4.6559_FLOAT,4.4533_FLOAT,
     *4.4956_FLOAT,4.4641_FLOAT,4.4366_FLOAT,4.4155_FLOAT,4.3936_FLOAT,
     *4.3764_FLOAT,4.1302_FLOAT,4.3398_FLOAT,4.2283_FLOAT,4.2796_FLOAT,
     *4.2547_FLOAT,4.2391_FLOAT,4.2296_FLOAT,4.0699_FLOAT,4.1083_FLOAT,
     *4.0319_FLOAT,3.9855_FLOAT,2.7676_FLOAT,2.8078_FLOAT,3.6725_FLOAT,
     *3.4804_FLOAT,3.3775_FLOAT,3.2411_FLOAT,3.1581_FLOAT,3.1983_FLOAT,
     *2.9973_FLOAT,2.8705_FLOAT,3.8070_FLOAT,3.7392_FLOAT,3.7668_FLOAT,
     *3.6263_FLOAT,3.5402_FLOAT,3.6807_FLOAT,3.4545_FLOAT,3.2962_FLOAT,
     *4.2283_FLOAT,4.1698_FLOAT,4.0240_FLOAT,3.9341_FLOAT,3.8711_FLOAT,
     *3.5489_FLOAT,3.7798_FLOAT/
      DATA T0AB(2801:2870)/3.7000_FLOAT,3.6654_FLOAT,3.6154_FLOAT,
     *3.5882_FLOAT,3.5472_FLOAT,3.7289_FLOAT,3.6510_FLOAT,3.6078_FLOAT,
     *3.5355_FLOAT,3.5963_FLOAT,3.4480_FLOAT,4.3587_FLOAT,4.3390_FLOAT,
     *4.1635_FLOAT,4.0536_FLOAT,3.7193_FLOAT,3.6529_FLOAT,3.5512_FLOAT,
     *3.4837_FLOAT,3.4400_FLOAT,3.4191_FLOAT,3.3891_FLOAT,3.6622_FLOAT,
     *3.8934_FLOAT,3.8235_FLOAT,3.7823_FLOAT,3.7292_FLOAT,3.8106_FLOAT,
     *3.6589_FLOAT,4.5535_FLOAT,4.6013_FLOAT,4.3961_FLOAT,4.4423_FLOAT,
     *4.4109_FLOAT,4.3835_FLOAT,4.3625_FLOAT,4.3407_FLOAT,4.3237_FLOAT,
     *4.0863_FLOAT,4.2835_FLOAT,4.1675_FLOAT,4.2272_FLOAT,4.2025_FLOAT,
     *4.1869_FLOAT,4.1774_FLOAT,4.0126_FLOAT,4.0460_FLOAT,3.9815_FLOAT,
     *3.9340_FLOAT,3.8955_FLOAT,2.6912_FLOAT,2.7604_FLOAT,3.6037_FLOAT,
     *3.4194_FLOAT,3.3094_FLOAT,3.1710_FLOAT,3.0862_FLOAT,3.1789_FLOAT,
     *2.9738_FLOAT,2.8427_FLOAT,3.7378_FLOAT,3.6742_FLOAT,3.6928_FLOAT,
     *3.5512_FLOAT,3.4614_FLOAT,3.4087_FLOAT,3.4201_FLOAT,3.2607_FLOAT,
     *4.1527_FLOAT,4.0977_FLOAT/
      DATA T0AB(2871:2940)/3.9523_FLOAT,3.8628_FLOAT,3.8002_FLOAT,
     *3.4759_FLOAT,3.7102_FLOAT,3.6466_FLOAT,3.6106_FLOAT,3.5580_FLOAT,
     *3.5282_FLOAT,3.4878_FLOAT,3.6547_FLOAT,3.5763_FLOAT,3.5289_FLOAT,
     *3.5086_FLOAT,3.5593_FLOAT,3.4099_FLOAT,4.2788_FLOAT,4.2624_FLOAT,
     *4.0873_FLOAT,3.9770_FLOAT,3.6407_FLOAT,3.5743_FLOAT,3.5178_FLOAT,
     *3.4753_FLOAT,3.3931_FLOAT,3.3694_FLOAT,3.3339_FLOAT,3.6002_FLOAT,
     *3.8164_FLOAT,3.7478_FLOAT,3.7028_FLOAT,3.6952_FLOAT,3.7669_FLOAT,
     *3.6137_FLOAT,4.4698_FLOAT,4.5488_FLOAT,4.3168_FLOAT,4.3646_FLOAT,
     *4.3338_FLOAT,4.3067_FLOAT,4.2860_FLOAT,4.2645_FLOAT,4.2478_FLOAT,
     *4.0067_FLOAT,4.2349_FLOAT,4.0958_FLOAT,4.1543_FLOAT,4.1302_FLOAT,
     *4.1141_FLOAT,4.1048_FLOAT,3.9410_FLOAT,3.9595_FLOAT,3.8941_FLOAT,
     *3.8465_FLOAT,3.8089_FLOAT,3.7490_FLOAT,2.7895_FLOAT,2.5849_FLOAT,
     *3.6484_FLOAT,3.0162_FLOAT,3.1267_FLOAT,3.2125_FLOAT,3.0043_FLOAT,
     *2.9572_FLOAT,2.8197_FLOAT,2.7261_FLOAT,3.7701_FLOAT,3.2446_FLOAT,
     *3.5239_FLOAT,3.4696_FLOAT/
      DATA T0AB(2941:3010)/3.4261_FLOAT,3.3508_FLOAT,3.1968_FLOAT,
     *3.0848_FLOAT,4.1496_FLOAT,3.6598_FLOAT,3.5111_FLOAT,3.4199_FLOAT,
     *3.3809_FLOAT,3.5382_FLOAT,3.2572_FLOAT,3.2100_FLOAT,3.1917_FLOAT,
     *3.1519_FLOAT,3.1198_FLOAT,3.1005_FLOAT,3.5071_FLOAT,3.5086_FLOAT,
     *3.5073_FLOAT,3.4509_FLOAT,3.3120_FLOAT,3.2082_FLOAT,4.2611_FLOAT,
     *3.8117_FLOAT,3.6988_FLOAT,3.5646_FLOAT,3.6925_FLOAT,3.6295_FLOAT,
     *3.5383_FLOAT,3.4910_FLOAT,3.4625_FLOAT,3.4233_FLOAT,3.4007_FLOAT,
     *3.2329_FLOAT,3.6723_FLOAT,3.6845_FLOAT,3.6876_FLOAT,3.6197_FLOAT,
     *3.4799_FLOAT,3.3737_FLOAT,4.4341_FLOAT,4.0525_FLOAT,3.9011_FLOAT,
     *3.8945_FLOAT,3.8635_FLOAT,3.8368_FLOAT,3.8153_FLOAT,3.7936_FLOAT,
     *3.7758_FLOAT,3.4944_FLOAT,3.4873_FLOAT,3.9040_FLOAT,3.7110_FLOAT,
     *3.6922_FLOAT,3.6799_FLOAT,3.6724_FLOAT,3.5622_FLOAT,3.6081_FLOAT,
     *3.5426_FLOAT,3.4922_FLOAT,3.4498_FLOAT,3.3984_FLOAT,3.4456_FLOAT,
     *2.7522_FLOAT,2.5524_FLOAT,3.5742_FLOAT,2.9508_FLOAT,3.0751_FLOAT,
     *3.0158_FLOAT,2.9644_FLOAT/
      DATA T0AB(3011:3080)/2.8338_FLOAT,2.7891_FLOAT,2.6933_FLOAT,
     *3.6926_FLOAT,3.1814_FLOAT,3.4528_FLOAT,3.4186_FLOAT,3.3836_FLOAT,
     *3.2213_FLOAT,3.1626_FLOAT,3.0507_FLOAT,4.0548_FLOAT,3.5312_FLOAT,
     *3.4244_FLOAT,3.3409_FLOAT,3.2810_FLOAT,3.4782_FLOAT,3.1905_FLOAT,
     *3.1494_FLOAT,3.1221_FLOAT,3.1128_FLOAT,3.0853_FLOAT,3.0384_FLOAT,
     *3.4366_FLOAT,3.4562_FLOAT,3.4638_FLOAT,3.3211_FLOAT,3.2762_FLOAT,
     *3.1730_FLOAT,4.1632_FLOAT,3.6825_FLOAT,3.5822_FLOAT,3.4870_FLOAT,
     *3.6325_FLOAT,3.5740_FLOAT,3.4733_FLOAT,3.4247_FLOAT,3.3969_FLOAT,
     *3.3764_FLOAT,3.3525_FLOAT,3.1984_FLOAT,3.5989_FLOAT,3.6299_FLOAT,
     *3.6433_FLOAT,3.4937_FLOAT,3.4417_FLOAT,3.3365_FLOAT,4.3304_FLOAT,
     *3.9242_FLOAT,3.7793_FLOAT,3.7623_FLOAT,3.7327_FLOAT,3.7071_FLOAT,
     *3.6860_FLOAT,3.6650_FLOAT,3.6476_FLOAT,3.3849_FLOAT,3.3534_FLOAT,
     *3.8216_FLOAT,3.5870_FLOAT,3.5695_FLOAT,3.5584_FLOAT,3.5508_FLOAT,
     *3.4856_FLOAT,3.5523_FLOAT,3.4934_FLOAT,3.4464_FLOAT,3.4055_FLOAT,
     *3.3551_FLOAT,3.3888_FLOAT/
      DATA T0AB(3081:3150)/3.3525_FLOAT,2.7202_FLOAT,2.5183_FLOAT,
     *3.4947_FLOAT,2.8731_FLOAT,3.0198_FLOAT,3.1457_FLOAT,2.9276_FLOAT,
     *2.7826_FLOAT,2.7574_FLOAT,2.6606_FLOAT,3.6090_FLOAT,3.0581_FLOAT,
     *3.3747_FLOAT,3.3677_FLOAT,3.3450_FLOAT,3.1651_FLOAT,3.1259_FLOAT,
     *3.0147_FLOAT,3.9498_FLOAT,3.3857_FLOAT,3.2917_FLOAT,3.2154_FLOAT,
     *3.1604_FLOAT,3.4174_FLOAT,3.0735_FLOAT,3.0342_FLOAT,3.0096_FLOAT,
     *3.0136_FLOAT,2.9855_FLOAT,2.9680_FLOAT,3.3604_FLOAT,3.4037_FLOAT,
     *3.4243_FLOAT,3.2633_FLOAT,3.1810_FLOAT,3.1351_FLOAT,4.0557_FLOAT,
     *3.5368_FLOAT,3.4526_FLOAT,3.3699_FLOAT,3.5707_FLOAT,3.5184_FLOAT,
     *3.4085_FLOAT,3.3595_FLOAT,3.3333_FLOAT,3.3143_FLOAT,3.3041_FLOAT,
     *3.1094_FLOAT,3.5193_FLOAT,3.5745_FLOAT,3.6025_FLOAT,3.4338_FLOAT,
     *3.3448_FLOAT,3.2952_FLOAT,4.2158_FLOAT,3.7802_FLOAT,3.6431_FLOAT,
     *3.6129_FLOAT,3.5853_FLOAT,3.5610_FLOAT,3.5406_FLOAT,3.5204_FLOAT,
     *3.5036_FLOAT,3.2679_FLOAT,3.2162_FLOAT,3.7068_FLOAT,3.4483_FLOAT,
     *3.4323_FLOAT,3.4221_FLOAT/
      DATA T0AB(3151:3220)/3.4138_FLOAT,3.3652_FLOAT,3.4576_FLOAT,
     *3.4053_FLOAT,3.3618_FLOAT,3.3224_FLOAT,3.2711_FLOAT,3.3326_FLOAT,
     *3.2950_FLOAT,3.2564_FLOAT,2.5315_FLOAT,2.6104_FLOAT,3.2734_FLOAT,
     *3.2299_FLOAT,3.1090_FLOAT,2.9942_FLOAT,2.9159_FLOAT,2.8324_FLOAT,
     *2.8350_FLOAT,2.7216_FLOAT,3.3994_FLOAT,3.4475_FLOAT,3.4354_FLOAT,
     *3.3438_FLOAT,3.2807_FLOAT,3.2169_FLOAT,3.2677_FLOAT,3.1296_FLOAT,
     *3.7493_FLOAT,3.8075_FLOAT,3.6846_FLOAT,3.6104_FLOAT,3.5577_FLOAT,
     *3.2052_FLOAT,3.4803_FLOAT,3.4236_FLOAT,3.3845_FLOAT,3.3640_FLOAT,
     *3.3365_FLOAT,3.3010_FLOAT,3.3938_FLOAT,3.3624_FLOAT,3.3440_FLOAT,
     *3.3132_FLOAT,3.4035_FLOAT,3.2754_FLOAT,3.8701_FLOAT,3.9523_FLOAT,
     *3.8018_FLOAT,3.7149_FLOAT,3.3673_FLOAT,3.3199_FLOAT,3.2483_FLOAT,
     *3.2069_FLOAT,3.1793_FLOAT,3.1558_FLOAT,3.1395_FLOAT,3.4097_FLOAT,
     *3.5410_FLOAT,3.5228_FLOAT,3.5116_FLOAT,3.4921_FLOAT,3.4781_FLOAT,
     *3.4690_FLOAT,4.0420_FLOAT,4.1759_FLOAT,4.0078_FLOAT,4.0450_FLOAT,
     *4.0189_FLOAT,3.9952_FLOAT/
      DATA T0AB(3221:3290)/3.9770_FLOAT,3.9583_FLOAT,3.9434_FLOAT,
     *3.7217_FLOAT,3.8228_FLOAT,3.7826_FLOAT,3.8640_FLOAT,3.8446_FLOAT,
     *3.8314_FLOAT,3.8225_FLOAT,3.6817_FLOAT,3.7068_FLOAT,3.6555_FLOAT,
     *3.6159_FLOAT,3.5831_FLOAT,3.5257_FLOAT,3.2133_FLOAT,3.1689_FLOAT,
     *3.1196_FLOAT,3.3599_FLOAT,2.9852_FLOAT,2.7881_FLOAT,3.5284_FLOAT,
     *3.3493_FLOAT,3.6958_FLOAT,3.3642_FLOAT,3.1568_FLOAT,3.0055_FLOAT,
     *2.9558_FLOAT,2.8393_FLOAT,3.6287_FLOAT,3.5283_FLOAT,4.1511_FLOAT,
     *3.8259_FLOAT,3.6066_FLOAT,3.4527_FLOAT,3.3480_FLOAT,3.2713_FLOAT,
     *3.9037_FLOAT,3.8361_FLOAT,3.8579_FLOAT,3.7311_FLOAT,3.6575_FLOAT,
     *3.5176_FLOAT,3.5693_FLOAT,3.5157_FLOAT,3.4814_FLOAT,3.4559_FLOAT,
     *3.4445_FLOAT,3.4160_FLOAT,4.1231_FLOAT,3.8543_FLOAT,3.6816_FLOAT,
     *3.5602_FLOAT,3.4798_FLOAT,3.4208_FLOAT,4.0542_FLOAT,4.0139_FLOAT,
     *4.0165_FLOAT,3.9412_FLOAT,3.7698_FLOAT,3.6915_FLOAT,3.6043_FLOAT,
     *3.5639_FLOAT,3.5416_FLOAT,3.5247_FLOAT,3.5153_FLOAT,3.5654_FLOAT,
     *4.2862_FLOAT,4.0437_FLOAT/
      DATA T0AB(3291:3360)/3.8871_FLOAT,3.7741_FLOAT,3.6985_FLOAT,
     *3.6413_FLOAT,4.2345_FLOAT,4.3663_FLOAT,4.3257_FLOAT,4.0869_FLOAT,
     *4.0612_FLOAT,4.0364_FLOAT,4.0170_FLOAT,3.9978_FLOAT,3.9834_FLOAT,
     *3.9137_FLOAT,3.8825_FLOAT,3.8758_FLOAT,3.9143_FLOAT,3.8976_FLOAT,
     *3.8864_FLOAT,3.8768_FLOAT,3.9190_FLOAT,4.1613_FLOAT,4.0566_FLOAT,
     *3.9784_FLOAT,3.9116_FLOAT,3.8326_FLOAT,3.7122_FLOAT,3.6378_FLOAT,
     *3.5576_FLOAT,3.5457_FLOAT,4.3127_FLOAT,3.1160_FLOAT,2.8482_FLOAT,
     *4.0739_FLOAT,3.3599_FLOAT,3.5698_FLOAT,3.5366_FLOAT,3.2854_FLOAT,
     *3.1039_FLOAT,2.9953_FLOAT,2.9192_FLOAT,4.1432_FLOAT,3.5320_FLOAT,
     *3.9478_FLOAT,4.0231_FLOAT,3.7509_FLOAT,3.5604_FLOAT,3.4340_FLOAT,
     *3.3426_FLOAT,4.3328_FLOAT,3.8288_FLOAT,3.7822_FLOAT,3.7909_FLOAT,
     *3.6907_FLOAT,3.6864_FLOAT,3.5793_FLOAT,3.5221_FLOAT,3.4883_FLOAT,
     *3.4649_FLOAT,3.4514_FLOAT,3.4301_FLOAT,3.9256_FLOAT,4.0596_FLOAT,
     *3.8307_FLOAT,3.6702_FLOAT,3.5651_FLOAT,3.4884_FLOAT,4.4182_FLOAT,
     *4.2516_FLOAT,3.9687_FLOAT/
      DATA T0AB(3361:3430)/3.9186_FLOAT,3.9485_FLOAT,3.8370_FLOAT,
     *3.7255_FLOAT,3.6744_FLOAT,3.6476_FLOAT,3.6295_FLOAT,3.6193_FLOAT,
     *3.5659_FLOAT,4.0663_FLOAT,4.2309_FLOAT,4.0183_FLOAT,3.8680_FLOAT,
     *3.7672_FLOAT,3.6923_FLOAT,4.5240_FLOAT,4.4834_FLOAT,4.1570_FLOAT,
     *4.3204_FLOAT,4.2993_FLOAT,4.2804_FLOAT,4.2647_FLOAT,4.2481_FLOAT,
     *4.2354_FLOAT,3.8626_FLOAT,3.8448_FLOAT,4.2267_FLOAT,4.1799_FLOAT,
     *4.1670_FLOAT,3.8738_FLOAT,3.8643_FLOAT,3.8796_FLOAT,4.0575_FLOAT,
     *4.0354_FLOAT,3.9365_FLOAT,3.8611_FLOAT,3.7847_FLOAT,3.7388_FLOAT,
     *3.6826_FLOAT,3.6251_FLOAT,3.5492_FLOAT,4.0889_FLOAT,4.2764_FLOAT,
     *3.1416_FLOAT,2.8325_FLOAT,3.7735_FLOAT,3.3787_FLOAT,3.4632_FLOAT,
     *3.5923_FLOAT,3.3214_FLOAT,3.1285_FLOAT,3.0147_FLOAT,2.9366_FLOAT,
     *3.8527_FLOAT,3.5602_FLOAT,3.8131_FLOAT,3.8349_FLOAT,3.7995_FLOAT,
     *3.5919_FLOAT,3.4539_FLOAT,3.3540_FLOAT,4.0654_FLOAT,3.8603_FLOAT,
     *3.7972_FLOAT,3.7358_FLOAT,3.7392_FLOAT,3.8157_FLOAT,3.6055_FLOAT,
     *3.5438_FLOAT,3.5089_FLOAT/
      DATA T0AB(3431:3500)/3.4853_FLOAT,3.4698_FLOAT,3.4508_FLOAT,
     *3.7882_FLOAT,3.8682_FLOAT,3.8837_FLOAT,3.7055_FLOAT,3.5870_FLOAT,
     *3.5000_FLOAT,4.1573_FLOAT,4.0005_FLOAT,3.9568_FLOAT,3.8936_FLOAT,
     *3.9990_FLOAT,3.9433_FLOAT,3.8172_FLOAT,3.7566_FLOAT,3.7246_FLOAT,
     *3.7033_FLOAT,3.6900_FLOAT,3.5697_FLOAT,3.9183_FLOAT,4.0262_FLOAT,
     *4.0659_FLOAT,3.8969_FLOAT,3.7809_FLOAT,3.6949_FLOAT,4.2765_FLOAT,
     *4.2312_FLOAT,4.1401_FLOAT,4.0815_FLOAT,4.0580_FLOAT,4.0369_FLOAT,
     *4.0194_FLOAT,4.0017_FLOAT,3.9874_FLOAT,3.8312_FLOAT,3.8120_FLOAT,
     *3.9454_FLOAT,3.9210_FLOAT,3.9055_FLOAT,3.8951_FLOAT,3.8866_FLOAT,
     *3.8689_FLOAT,3.9603_FLOAT,3.9109_FLOAT,3.9122_FLOAT,3.8233_FLOAT,
     *3.7438_FLOAT,3.7436_FLOAT,3.6981_FLOAT,3.6555_FLOAT,3.5452_FLOAT,
     *3.9327_FLOAT,4.0658_FLOAT,4.1175_FLOAT,2.9664_FLOAT,2.8209_FLOAT,
     *3.5547_FLOAT,3.3796_FLOAT,3.3985_FLOAT,3.3164_FLOAT,3.2364_FLOAT,
     *3.1956_FLOAT,3.0370_FLOAT,2.9313_FLOAT,3.6425_FLOAT,3.5565_FLOAT,
     *3.7209_FLOAT,3.7108_FLOAT/
      DATA T0AB(3501:3570)/3.6639_FLOAT,3.6484_FLOAT,3.4745_FLOAT,
     *3.3492_FLOAT,3.8755_FLOAT,4.2457_FLOAT,3.7758_FLOAT,3.7161_FLOAT,
     *3.6693_FLOAT,3.6155_FLOAT,3.5941_FLOAT,3.5643_FLOAT,3.5292_FLOAT,
     *3.4950_FLOAT,3.4720_FLOAT,3.4503_FLOAT,3.6936_FLOAT,3.7392_FLOAT,
     *3.7388_FLOAT,3.7602_FLOAT,3.6078_FLOAT,3.4960_FLOAT,3.9800_FLOAT,
     *4.3518_FLOAT,4.2802_FLOAT,3.8580_FLOAT,3.8056_FLOAT,3.7527_FLOAT,
     *3.7019_FLOAT,3.6615_FLOAT,3.5768_FLOAT,3.5330_FLOAT,3.5038_FLOAT,
     *3.5639_FLOAT,3.8192_FLOAT,3.8883_FLOAT,3.9092_FLOAT,3.9478_FLOAT,
     *3.7995_FLOAT,3.6896_FLOAT,4.1165_FLOAT,4.5232_FLOAT,4.4357_FLOAT,
     *4.4226_FLOAT,4.4031_FLOAT,4.3860_FLOAT,4.3721_FLOAT,4.3580_FLOAT,
     *4.3466_FLOAT,4.2036_FLOAT,4.2037_FLOAT,3.8867_FLOAT,4.2895_FLOAT,
     *4.2766_FLOAT,4.2662_FLOAT,4.2598_FLOAT,3.8408_FLOAT,3.9169_FLOAT,
     *3.8681_FLOAT,3.8250_FLOAT,3.7855_FLOAT,3.7501_FLOAT,3.6753_FLOAT,
     *3.5499_FLOAT,3.4872_FLOAT,3.5401_FLOAT,3.8288_FLOAT,3.9217_FLOAT,
     *3.9538_FLOAT,4.0054_FLOAT/
      DATA T0AB(3571:3640)/2.8388_FLOAT,2.7890_FLOAT,3.4329_FLOAT,
     *3.5593_FLOAT,3.3488_FLOAT,3.2486_FLOAT,3.1615_FLOAT,3.1000_FLOAT,
     *3.0394_FLOAT,2.9165_FLOAT,3.5267_FLOAT,3.7479_FLOAT,3.6650_FLOAT,
     *3.6263_FLOAT,3.5658_FLOAT,3.5224_FLOAT,3.4762_FLOAT,3.3342_FLOAT,
     *3.7738_FLOAT,4.0333_FLOAT,3.9568_FLOAT,3.8975_FLOAT,3.8521_FLOAT,
     *3.4929_FLOAT,3.7830_FLOAT,3.7409_FLOAT,3.7062_FLOAT,3.6786_FLOAT,
     *3.6471_FLOAT,3.6208_FLOAT,3.6337_FLOAT,3.6519_FLOAT,3.6363_FLOAT,
     *3.6278_FLOAT,3.6110_FLOAT,3.4825_FLOAT,3.8795_FLOAT,4.1448_FLOAT,
     *4.0736_FLOAT,4.0045_FLOAT,3.6843_FLOAT,3.6291_FLOAT,3.5741_FLOAT,
     *3.5312_FLOAT,3.4974_FLOAT,3.4472_FLOAT,3.4034_FLOAT,3.7131_FLOAT,
     *3.7557_FLOAT,3.7966_FLOAT,3.8005_FLOAT,3.8068_FLOAT,3.8015_FLOAT,
     *3.6747_FLOAT,4.0222_FLOAT,4.3207_FLOAT,4.2347_FLOAT,4.2191_FLOAT,
     *4.1990_FLOAT,4.1811_FLOAT,4.1666_FLOAT,4.1521_FLOAT,4.1401_FLOAT,
     *3.9970_FLOAT,3.9943_FLOAT,3.9592_FLOAT,4.0800_FLOAT,4.0664_FLOAT,
     *4.0559_FLOAT,4.0488_FLOAT/
      DATA T0AB(3641:3710)/3.9882_FLOAT,4.0035_FLOAT,3.9539_FLOAT,
     *3.9138_FLOAT,3.8798_FLOAT,3.8355_FLOAT,3.5359_FLOAT,3.4954_FLOAT,
     *3.3962_FLOAT,3.5339_FLOAT,3.7595_FLOAT,3.8250_FLOAT,3.8408_FLOAT,
     *3.8600_FLOAT,3.8644_FLOAT,2.7412_FLOAT,2.7489_FLOAT,3.3374_FLOAT,
     *3.3950_FLOAT,3.3076_FLOAT,3.1910_FLOAT,3.0961_FLOAT,3.0175_FLOAT,
     *3.0280_FLOAT,2.8929_FLOAT,3.4328_FLOAT,3.5883_FLOAT,3.6227_FLOAT,
     *3.5616_FLOAT,3.4894_FLOAT,3.4241_FLOAT,3.3641_FLOAT,3.3120_FLOAT,
     *3.6815_FLOAT,3.8789_FLOAT,3.8031_FLOAT,3.7413_FLOAT,3.6939_FLOAT,
     *3.4010_FLOAT,3.6225_FLOAT,3.5797_FLOAT,3.5443_FLOAT,3.5139_FLOAT,
     *3.4923_FLOAT,3.4642_FLOAT,3.5860_FLOAT,3.5849_FLOAT,3.5570_FLOAT,
     *3.5257_FLOAT,3.4936_FLOAT,3.4628_FLOAT,3.7874_FLOAT,3.9916_FLOAT,
     *3.9249_FLOAT,3.8530_FLOAT,3.5932_FLOAT,3.5355_FLOAT,3.4757_FLOAT,
     *3.4306_FLOAT,3.3953_FLOAT,3.3646_FLOAT,3.3390_FLOAT,3.5637_FLOAT,
     *3.7053_FLOAT,3.7266_FLOAT,3.7177_FLOAT,3.6996_FLOAT,3.6775_FLOAT,
     *3.6558_FLOAT,3.9331_FLOAT/
      DATA T0AB(3711:3780)/4.1655_FLOAT,4.0879_FLOAT,4.0681_FLOAT,
     *4.0479_FLOAT,4.0299_FLOAT,4.0152_FLOAT,4.0006_FLOAT,3.9883_FLOAT,
     *3.8500_FLOAT,3.8359_FLOAT,3.8249_FLOAT,3.9269_FLOAT,3.9133_FLOAT,
     *3.9025_FLOAT,3.8948_FLOAT,3.8422_FLOAT,3.8509_FLOAT,3.7990_FLOAT,
     *3.7570_FLOAT,3.7219_FLOAT,3.6762_FLOAT,3.4260_FLOAT,3.3866_FLOAT,
     *3.3425_FLOAT,3.5294_FLOAT,3.7022_FLOAT,3.7497_FLOAT,3.7542_FLOAT,
     *3.7494_FLOAT,3.7370_FLOAT,3.7216_FLOAT,3.4155_FLOAT,3.0522_FLOAT,
     *4.2541_FLOAT,3.8218_FLOAT,4.0438_FLOAT,3.5875_FLOAT,3.3286_FLOAT,
     *3.1682_FLOAT,3.0566_FLOAT,2.9746_FLOAT,4.3627_FLOAT,4.0249_FLOAT,
     *4.6947_FLOAT,4.1718_FLOAT,3.8639_FLOAT,3.6735_FLOAT,3.5435_FLOAT,
     *3.4479_FLOAT,4.6806_FLOAT,4.3485_FLOAT,4.2668_FLOAT,4.1690_FLOAT,
     *4.1061_FLOAT,4.1245_FLOAT,4.0206_FLOAT,3.9765_FLOAT,3.9458_FLOAT,
     *3.9217_FLOAT,3.9075_FLOAT,3.8813_FLOAT,3.9947_FLOAT,4.1989_FLOAT,
     .3.9507_FLOAT,3.7960_FLOAT,3.6925_FLOAT,3.6150_FLOAT,4.8535_FLOAT,
     *4.5642_FLOAT,4.4134_FLOAT/
      DATA T0AB(3781:3850)/4.3688_FLOAT,4.3396_FLOAT,4.2879_FLOAT,
     *4.2166_FLOAT,4.1888_FLOAT,4.1768_FLOAT,4.1660_FLOAT,4.1608_FLOAT,
     *4.0745_FLOAT,4.2289_FLOAT,4.4863_FLOAT,4.2513_FLOAT,4.0897_FLOAT,
     *3.9876_FLOAT,3.9061_FLOAT,5.0690_FLOAT,5.0446_FLOAT,4.6186_FLOAT,
     *4.6078_FLOAT,4.5780_FLOAT,4.5538_FLOAT,4.5319_FLOAT,4.5101_FLOAT,
     *4.4945_FLOAT,4.1912_FLOAT,4.2315_FLOAT,4.5534_FLOAT,4.4373_FLOAT,
     *4.4224_FLOAT,4.4120_FLOAT,4.4040_FLOAT,4.2634_FLOAT,4.7770_FLOAT,
     *4.6890_FLOAT,4.6107_FLOAT,4.5331_FLOAT,4.4496_FLOAT,4.4082_FLOAT,
     *4.3095_FLOAT,4.2023_FLOAT,4.0501_FLOAT,4.2595_FLOAT,4.5497_FLOAT,
     *4.3056_FLOAT,4.1506_FLOAT,4.0574_FLOAT,3.9725_FLOAT,5.0796_FLOAT,
     *3.0548_FLOAT,3.3206_FLOAT,3.8132_FLOAT,3.9720_FLOAT,3.7675_FLOAT,
     *3.7351_FLOAT,3.5167_FLOAT,3.5274_FLOAT,3.3085_FLOAT,3.1653_FLOAT,
     *3.9500_FLOAT,4.1730_FLOAT,4.0613_FLOAT,4.1493_FLOAT,3.8823_FLOAT,
     *4.0537_FLOAT,3.8200_FLOAT,3.6582_FLOAT,4.3422_FLOAT,4.5111_FLOAT,
     *4.3795_FLOAT,4.3362_FLOAT/
      DATA T0AB(3851:3920)/4.2751_FLOAT,3.7103_FLOAT,4.1973_FLOAT,
     *4.1385_FLOAT,4.1129_FLOAT,4.0800_FLOAT,4.0647_FLOAT,4.0308_FLOAT,
     *4.0096_FLOAT,4.1619_FLOAT,3.9360_FLOAT,4.1766_FLOAT,3.9705_FLOAT,
     *3.8262_FLOAT,4.5348_FLOAT,4.7025_FLOAT,4.5268_FLOAT,4.5076_FLOAT,
     *3.9562_FLOAT,3.9065_FLOAT,3.8119_FLOAT,3.7605_FLOAT,3.7447_FLOAT,
     *3.7119_FLOAT,3.6916_FLOAT,4.1950_FLOAT,4.2110_FLOAT,4.3843_FLOAT,
     *4.1631_FLOAT,4.4427_FLOAT,4.2463_FLOAT,4.1054_FLOAT,4.7693_FLOAT,
     *5.0649_FLOAT,4.7365_FLOAT,4.7761_FLOAT,4.7498_FLOAT,4.7272_FLOAT,
     *4.7076_FLOAT,4.6877_FLOAT,4.6730_FLOAT,4.4274_FLOAT,4.5473_FLOAT,
     *4.5169_FLOAT,4.5975_FLOAT,4.5793_FLOAT,4.5667_FLOAT,4.5559_FLOAT,
     *4.3804_FLOAT,4.6920_FLOAT,4.6731_FLOAT,4.6142_FLOAT,4.5600_FLOAT,
     *4.4801_FLOAT,4.0149_FLOAT,3.8856_FLOAT,3.7407_FLOAT,4.1545_FLOAT,
     *4.2253_FLOAT,4.4229_FLOAT,4.1923_FLOAT,4.5022_FLOAT,4.3059_FLOAT,
     *4.1591_FLOAT,4.7883_FLOAT,4.9294_FLOAT,3.3850_FLOAT,3.4208_FLOAT,
     *3.7004_FLOAT,3.8800_FLOAT/
      DATA T0AB(3921:3990)/3.9886_FLOAT,3.9040_FLOAT,3.6719_FLOAT,
     *3.6547_FLOAT,3.4625_FLOAT,3.3370_FLOAT,3.8394_FLOAT,4.0335_FLOAT,
     *4.2373_FLOAT,4.3023_FLOAT,4.0306_FLOAT,4.1408_FLOAT,3.9297_FLOAT,
     *3.7857_FLOAT,4.1907_FLOAT,4.3230_FLOAT,4.2664_FLOAT,4.2173_FLOAT,
     *4.1482_FLOAT,3.6823_FLOAT,4.0711_FLOAT,4.0180_FLOAT,4.0017_FLOAT,
     *3.9747_FLOAT,3.9634_FLOAT,3.9383_FLOAT,4.1993_FLOAT,4.3205_FLOAT,
     *4.0821_FLOAT,4.2547_FLOAT,4.0659_FLOAT,3.9359_FLOAT,4.3952_FLOAT,
     *4.5176_FLOAT,4.3888_FLOAT,4.3607_FLOAT,3.9583_FLOAT,3.9280_FLOAT,
     *3.8390_FLOAT,3.7971_FLOAT,3.7955_FLOAT,3.7674_FLOAT,3.7521_FLOAT,
     *4.1062_FLOAT,4.3633_FLOAT,4.2991_FLOAT,4.2767_FLOAT,4.4857_FLOAT,
     *4.3039_FLOAT,4.1762_FLOAT,4.6197_FLOAT,4.8654_FLOAT,4.6633_FLOAT,
     *4.5878_FLOAT,4.5640_FLOAT,4.5422_FLOAT,4.5231_FLOAT,4.5042_FLOAT,
     *4.4901_FLOAT,4.3282_FLOAT,4.3978_FLOAT,4.3483_FLOAT,4.4202_FLOAT,
     *4.4039_FLOAT,4.3926_FLOAT,4.3807_FLOAT,4.2649_FLOAT,4.6135_FLOAT,
     *4.5605_FLOAT,4.5232_FLOAT/
      DATA T0AB(3991:4060)/4.4676_FLOAT,4.3948_FLOAT,4.0989_FLOAT,
     *3.9864_FLOAT,3.8596_FLOAT,4.0942_FLOAT,4.2720_FLOAT,4.3270_FLOAT,
     *4.3022_FLOAT,4.5410_FLOAT,4.3576_FLOAT,4.2235_FLOAT,4.6545_FLOAT,
     *4.7447_FLOAT,4.7043_FLOAT,3.0942_FLOAT,3.2075_FLOAT,3.5152_FLOAT,
     *3.6659_FLOAT,3.8289_FLOAT,3.7459_FLOAT,3.5156_FLOAT,3.5197_FLOAT,
     *3.3290_FLOAT,3.2069_FLOAT,3.6702_FLOAT,3.8448_FLOAT,4.0340_FLOAT,
     *3.9509_FLOAT,3.8585_FLOAT,3.9894_FLOAT,3.7787_FLOAT,3.6365_FLOAT,
     *4.1425_FLOAT,4.1618_FLOAT,4.0940_FLOAT,4.0466_FLOAT,3.9941_FLOAT,
     *3.5426_FLOAT,3.8952_FLOAT,3.8327_FLOAT,3.8126_FLOAT,3.7796_FLOAT,
     *3.7635_FLOAT,3.7356_FLOAT,4.0047_FLOAT,3.9655_FLOAT,3.9116_FLOAT,
     *4.1010_FLOAT,3.9102_FLOAT,3.7800_FLOAT,4.2964_FLOAT,4.3330_FLOAT,
     *4.2622_FLOAT,4.2254_FLOAT,3.8195_FLOAT,3.7560_FLOAT,3.6513_FLOAT,
     *3.5941_FLOAT,3.5810_FLOAT,3.5420_FLOAT,3.5178_FLOAT,3.8861_FLOAT,
     *4.1459_FLOAT,4.1147_FLOAT,4.0772_FLOAT,4.3120_FLOAT,4.1207_FLOAT,
     *3.9900_FLOAT,4.4733_FLOAT/
      DATA T0AB(4061:4130)/4.6157_FLOAT,4.4580_FLOAT,4.4194_FLOAT,
     *4.3954_FLOAT,4.3739_FLOAT,4.3531_FLOAT,4.3343_FLOAT,4.3196_FLOAT,
     *4.2140_FLOAT,4.2339_FLOAT,4.1738_FLOAT,4.2458_FLOAT,4.2278_FLOAT,
     *4.2158_FLOAT,4.2039_FLOAT,4.1658_FLOAT,4.3595_FLOAT,4.2857_FLOAT,
     *4.2444_FLOAT,4.1855_FLOAT,4.1122_FLOAT,3.7839_FLOAT,3.6879_FLOAT,
     *3.5816_FLOAT,3.8633_FLOAT,4.1585_FLOAT,4.1402_FLOAT,4.1036_FLOAT,
     *4.3694_FLOAT,4.1735_FLOAT,4.0368_FLOAT,4.5095_FLOAT,4.5538_FLOAT,
     *4.5240_FLOAT,4.4252_FLOAT,3.0187_FLOAT,3.1918_FLOAT,3.5127_FLOAT,
     *3.6875_FLOAT,3.7404_FLOAT,3.6943_FLOAT,3.4702_FLOAT,3.4888_FLOAT,
     *3.2914_FLOAT,3.1643_FLOAT,3.6669_FLOAT,3.8724_FLOAT,3.9940_FLOAT,
     *4.0816_FLOAT,3.8054_FLOAT,3.9661_FLOAT,3.7492_FLOAT,3.6024_FLOAT,
     *4.0428_FLOAT,4.1951_FLOAT,4.1466_FLOAT,4.0515_FLOAT,4.0075_FLOAT,
     *3.5020_FLOAT,3.9158_FLOAT,3.8546_FLOAT,3.8342_FLOAT,3.8008_FLOAT,
     *3.7845_FLOAT,3.7549_FLOAT,3.9602_FLOAT,3.8872_FLOAT,3.8564_FLOAT,
     *4.0793_FLOAT,3.8835_FLOAT/
      DATA T0AB(4131:4200)/3.7495_FLOAT,4.2213_FLOAT,4.3704_FLOAT,
     *4.3300_FLOAT,4.2121_FLOAT,3.7643_FLOAT,3.7130_FLOAT,3.6144_FLOAT,
     *3.5599_FLOAT,3.5474_FLOAT,3.5093_FLOAT,3.4853_FLOAT,3.9075_FLOAT,
     *4.1115_FLOAT,4.0473_FLOAT,4.0318_FLOAT,4.2999_FLOAT,4.1050_FLOAT,
     *3.9710_FLOAT,4.4320_FLOAT,4.6706_FLOAT,4.5273_FLOAT,4.4581_FLOAT,
     *4.4332_FLOAT,4.4064_FLOAT,4.3873_FLOAT,4.3684_FLOAT,4.3537_FLOAT,
     *4.2728_FLOAT,4.2549_FLOAT,4.2032_FLOAT,4.2794_FLOAT,4.2613_FLOAT,
     *4.2491_FLOAT,4.2375_FLOAT,4.2322_FLOAT,4.3665_FLOAT,4.3061_FLOAT,
     *4.2714_FLOAT,4.2155_FLOAT,4.1416_FLOAT,3.7660_FLOAT,3.6628_FLOAT,
     *3.5476_FLOAT,3.8790_FLOAT,4.1233_FLOAT,4.0738_FLOAT,4.0575_FLOAT,
     *4.3575_FLOAT,4.1586_FLOAT,4.0183_FLOAT,4.4593_FLOAT,4.5927_FLOAT,
     *4.4865_FLOAT,4.3813_FLOAT,4.4594_FLOAT,2.9875_FLOAT,3.1674_FLOAT,
     *3.4971_FLOAT,3.6715_FLOAT,3.7114_FLOAT,3.6692_FLOAT,3.4446_FLOAT,
     *3.4676_FLOAT,3.2685_FLOAT,3.1405_FLOAT,3.6546_FLOAT,3.8579_FLOAT,
     *3.9637_FLOAT,4.0581_FLOAT/
      DATA T0AB(4201:4270)/3.7796_FLOAT,3.9463_FLOAT,3.7275_FLOAT,
     *3.5792_FLOAT,4.0295_FLOAT,4.1824_FLOAT,4.1247_FLOAT,4.0357_FLOAT,
     *3.9926_FLOAT,3.4827_FLOAT,3.9007_FLOAT,3.8392_FLOAT,3.8191_FLOAT,
     *3.7851_FLOAT,3.7687_FLOAT,3.7387_FLOAT,3.9290_FLOAT,3.8606_FLOAT,
     *3.8306_FLOAT,4.0601_FLOAT,3.8625_FLOAT,3.7269_FLOAT,4.2062_FLOAT,
     *4.3566_FLOAT,4.3022_FLOAT,4.1929_FLOAT,3.7401_FLOAT,3.6888_FLOAT,
     *3.5900_FLOAT,3.5350_FLOAT,3.5226_FLOAT,3.4838_FLOAT,3.4594_FLOAT,
     *3.8888_FLOAT,4.0813_FLOAT,4.0209_FLOAT,4.0059_FLOAT,4.2810_FLOAT,
     *4.0843_FLOAT,3.9486_FLOAT,4.4162_FLOAT,4.6542_FLOAT,4.5005_FLOAT,
     *4.4444_FLOAT,4.4196_FLOAT,4.3933_FLOAT,4.3741_FLOAT,4.3552_FLOAT,
     *4.3406_FLOAT,4.2484_FLOAT,4.2413_FLOAT,4.1907_FLOAT,4.2656_FLOAT,
     *4.2474_FLOAT,4.2352_FLOAT,4.2236_FLOAT,4.2068_FLOAT,4.3410_FLOAT,
     *4.2817_FLOAT,4.2479_FLOAT,4.1921_FLOAT,4.1182_FLOAT,3.7346_FLOAT,
     *3.6314_FLOAT,3.5168_FLOAT,3.8582_FLOAT,4.0927_FLOAT,4.0469_FLOAT,
     *4.0313_FLOAT,4.3391_FLOAT/
      DATA T0AB(4271:4340)/4.1381_FLOAT,3.9962_FLOAT,4.4429_FLOAT,
     *4.5787_FLOAT,4.4731_FLOAT,4.3588_FLOAT,4.4270_FLOAT,4.3957_FLOAT,
     *2.9659_FLOAT,3.1442_FLOAT,3.4795_FLOAT,3.6503_FLOAT,3.6814_FLOAT,
     *3.6476_FLOAT,3.4222_FLOAT,3.4491_FLOAT,3.2494_FLOAT,3.1209_FLOAT,
     *3.6324_FLOAT,3.8375_FLOAT,3.9397_FLOAT,3.8311_FLOAT,3.7581_FLOAT,
     *3.9274_FLOAT,3.7085_FLOAT,3.5598_FLOAT,4.0080_FLOAT,4.1641_FLOAT,
     *4.1057_FLOAT,4.0158_FLOAT,3.9726_FLOAT,3.4667_FLOAT,3.8802_FLOAT,
     *3.8188_FLOAT,3.7989_FLOAT,3.7644_FLOAT,3.7474_FLOAT,3.7173_FLOAT,
     *3.9049_FLOAT,3.8424_FLOAT,3.8095_FLOAT,4.0412_FLOAT,3.8436_FLOAT,
     *3.7077_FLOAT,4.1837_FLOAT,4.3366_FLOAT,4.2816_FLOAT,4.1686_FLOAT,
     *3.7293_FLOAT,3.6709_FLOAT,3.5700_FLOAT,3.5153_FLOAT,3.5039_FLOAT,
     *3.4684_FLOAT,3.4437_FLOAT,3.8663_FLOAT,4.0575_FLOAT,4.0020_FLOAT,
     *3.9842_FLOAT,4.2612_FLOAT,4.0643_FLOAT,3.9285_FLOAT,4.3928_FLOAT,
     *4.6308_FLOAT,4.4799_FLOAT,4.4244_FLOAT,4.3996_FLOAT,4.3737_FLOAT,
     *4.3547_FLOAT,4.3358_FLOAT/
      DATA T0AB(4341:4410)/4.3212_FLOAT,4.2275_FLOAT,4.2216_FLOAT,
     *4.1676_FLOAT,4.2465_FLOAT,4.2283_FLOAT,4.2161_FLOAT,4.2045_FLOAT,
     *4.1841_FLOAT,4.3135_FLOAT,4.2562_FLOAT,4.2226_FLOAT,4.1667_FLOAT,
     *4.0932_FLOAT,3.7134_FLOAT,3.6109_FLOAT,3.4962_FLOAT,3.8352_FLOAT,
     *4.0688_FLOAT,4.0281_FLOAT,4.0099_FLOAT,4.3199_FLOAT,4.1188_FLOAT,
     *3.9768_FLOAT,4.4192_FLOAT,4.5577_FLOAT,4.4516_FLOAT,4.3365_FLOAT,
     *4.4058_FLOAT,4.3745_FLOAT,4.3539_FLOAT,2.8763_FLOAT,3.1294_FLOAT,
     *3.5598_FLOAT,3.7465_FLOAT,3.5659_FLOAT,3.5816_FLOAT,3.3599_FLOAT,
     *3.4024_FLOAT,3.1877_FLOAT,3.0484_FLOAT,3.7009_FLOAT,3.9451_FLOAT,
     *3.8465_FLOAT,3.9873_FLOAT,3.7079_FLOAT,3.9083_FLOAT,3.6756_FLOAT,
     *3.5150_FLOAT,4.0829_FLOAT,4.2780_FLOAT,4.1511_FLOAT,4.1260_FLOAT,
     *4.0571_FLOAT,3.4865_FLOAT,3.9744_FLOAT,3.9150_FLOAT,3.8930_FLOAT,
     *3.8578_FLOAT,3.8402_FLOAT,3.8073_FLOAT,3.7977_FLOAT,4.0036_FLOAT,
     *3.7604_FLOAT,4.0288_FLOAT,3.8210_FLOAT,3.6757_FLOAT,4.2646_FLOAT,
     *4.4558_FLOAT,4.2862_FLOAT/
      DATA T0AB(4411:4465)/4.2122_FLOAT,3.7088_FLOAT,3.6729_FLOAT,
     *3.5800_FLOAT,3.5276_FLOAT,3.5165_FLOAT,3.4783_FLOAT,3.4539_FLOAT,
     *3.9553_FLOAT,3.9818_FLOAT,4.2040_FLOAT,3.9604_FLOAT,4.2718_FLOAT,
     *4.0689_FLOAT,3.9253_FLOAT,4.4869_FLOAT,4.7792_FLOAT,4.4918_FLOAT,
     *4.5342_FLOAT,4.5090_FLOAT,4.4868_FLOAT,4.4680_FLOAT,4.4486_FLOAT,
     *4.4341_FLOAT,4.2023_FLOAT,4.3122_FLOAT,4.2710_FLOAT,4.3587_FLOAT,
     *4.3407_FLOAT,4.3281_FLOAT,4.3174_FLOAT,4.1499_FLOAT,4.3940_FLOAT,
     *4.3895_FLOAT,4.3260_FLOAT,4.2725_FLOAT,4.1961_FLOAT,3.7361_FLOAT,
     *3.6193_FLOAT,3.4916_FLOAT,3.9115_FLOAT,3.9914_FLOAT,3.9809_FLOAT,
     *3.9866_FLOAT,4.3329_FLOAT,4.1276_FLOAT,3.9782_FLOAT,4.5097_FLOAT,
     *4.6769_FLOAT,4.5158_FLOAT,4.3291_FLOAT,4.3609_FLOAT,4.3462_FLOAT,
     *4.3265_FLOAT,4.4341_FLOAT/
      ANGTOAU=1._FLOAT/PAR(32)
      k=0
      do i=1,max_elem
         do j=1,i
            k=k+1
            R0=T0ab(k)*ANGTOAU
            r(i,j)=R0
            r(j,i)=R0
         enddo
      enddo

      end subroutine setr0ab

      subroutine loadoldpar(max_elem,maxc,c6ab,r0ab)
      USE NUMBERS
      USE PARINF_MODULE
      IMPLICIT REAL(FLOAT)(A-H,O-Z)
      REAL(FLOAT) :: r0ab(max_elem,max_elem),c6ab(max_elem,max_elem,
     *maxc,maxc,3),c6(86),r0(86)
      DATA R0/0.91_FLOAT,0.92_FLOAT,0.75_FLOAT,1.28_FLOAT,1.35_FLOAT,
     *1.32_FLOAT,1.27_FLOAT,1.22_FLOAT,1.17_FLOAT,1.13_FLOAT,
     *1.04_FLOAT,1.24_FLOAT,1.49_FLOAT,1.56_FLOAT,1.55_FLOAT,1.53_FLOAT,
     *1.49_FLOAT,1.45_FLOAT,1.35_FLOAT,1.34_FLOAT,
     *1.42_FLOAT,1.42_FLOAT,1.42_FLOAT,1.42_FLOAT,1.42_FLOAT,1.42_FLOAT,
     *1.42_FLOAT,1.42_FLOAT,1.42_FLOAT,1.42_FLOAT,1.50_FLOAT,1.57_FLOAT,
     *1.60_FLOAT,1.61_FLOAT,1.59_FLOAT,1.57_FLOAT,1.48_FLOAT,1.46_FLOAT,
     *1.49_FLOAT,1.49_FLOAT,1.49_FLOAT,1.49_FLOAT,1.49_FLOAT,1.49_FLOAT,
     *1.49_FLOAT,1.49_FLOAT,1.49_FLOAT,1.49_FLOAT,1.52_FLOAT,1.64_FLOAT,
     *1.71_FLOAT,1.72_FLOAT,1.72_FLOAT,1.71_FLOAT,
     *1.638_FLOAT,1.602_FLOAT,1.564_FLOAT,1.594_FLOAT,1.594_FLOAT,
     *1.594_FLOAT,1.594_FLOAT,1.594_FLOAT,1.594_FLOAT,1.594_FLOAT,
     *1.594_FLOAT,1.594_FLOAT,1.594_FLOAT,1.594_FLOAT,1.594_FLOAT,
     *1.594_FLOAT,1.594_FLOAT,
     *1.625_FLOAT,1.611_FLOAT,1.611_FLOAT,1.611_FLOAT,1.611_FLOAT,
     *1.611_FLOAT,1.611_FLOAT,1.611_FLOAT,1.598_FLOAT,1.805_FLOAT,
     *1.767_FLOAT,1.725_FLOAT,1.823_FLOAT,1.810_FLOAT,1.749_FLOAT/
      DATA c6 /0.14_FLOAT,0.08_FLOAT,1.61_FLOAT,1.61_FLOAT,3.13_FLOAT,
     *1.75_FLOAT,1.23_FLOAT,0.70_FLOAT,0.75_FLOAT,0.63_FLOAT,5.71_FLOAT,
     *5.71_FLOAT,10.79_FLOAT,9.23_FLOAT,7.84_FLOAT,5.57_FLOAT,5.07_FLOAT
     *,4.61_FLOAT,10.8_FLOAT,10.8_FLOAT,10.8_FLOAT,10.8_FLOAT,10.8_FLOAT
     *,10.8_FLOAT,10.8_FLOAT,10.8_FLOAT,10.8_FLOAT,10.8_FLOAT,10.8_FLOAT
     *,10.8_FLOAT,16.99_FLOAT,17.10_FLOAT,16.37_FLOAT,12.64_FLOAT,
     *12.47_FLOAT,12.01_FLOAT,24.67_FLOAT,24.67_FLOAT,24.67_FLOAT,
     *24.67_FLOAT,24.67_FLOAT,24.67_FLOAT,24.67_FLOAT,24.67_FLOAT,
     *24.67_FLOAT,24.67_FLOAT,24.67_FLOAT,24.67_FLOAT,37.32_FLOAT,
     *38.71_FLOAT,38.44_FLOAT,31.74_FLOAT,31.50_FLOAT,29.99_FLOAT,
     *315.275_FLOAT,226.994_FLOAT,176.252_FLOAT,140.68_FLOAT,
     *140.68_FLOAT,140.68_FLOAT,140.68_FLOAT,140.68_FLOAT,140.68_FLOAT,
     *140.68_FLOAT,140.68_FLOAT,140.68_FLOAT,140.68_FLOAT,140.68_FLOAT,
     *140.68_FLOAT,140.68_FLOAT,140.68_FLOAT,105.112_FLOAT,81.24_FLOAT,
     *81.24_FLOAT,81.24_FLOAT,81.24_FLOAT,81.24_FLOAT,81.24_FLOAT,
     *81.24_FLOAT,57.364_FLOAT,57.254_FLOAT,63.162_FLOAT,63.540_FLOAT,
     *55.283_FLOAT,57.171_FLOAT,56.64_FLOAT/
      c6ab = -1._FLOAT
      ANGTOAU=1._FLOAT/PAR(32)
      do i=1,86
         do j=1,i
            R1=(r0(i)+r0(j))*ANGTOAU
            r0ab(i,j)=R1
            r0ab(j,i)=R1
            R1=sqrt(c6(i)*c6(j))
            c6ab(i,j,1,1,1)=R1
            c6ab(j,i,1,1,1)=R1
         enddo
      enddo
      end subroutine loadoldpar

      subroutine JGB_limit(iat,jat,iadr,jadr)
      USE NUMBERS
      IMPLICIT REAL(FLOAT)(A-H,O-Z)
      iadr=1
      jadr=1
      i=100
 10   if(iat.gt.100) then
         iat=iat-100
         iadr=iadr+1
         goto 10
      endif

      i=100
 20   if(jat.gt.100) then
         jat=jat-100
         jadr=jadr+1
         goto 20
      endif
      end subroutine JGB_limit
