Subroutine gcp_input
  use strings
  Use Numbers
  Use Parinf_module
  Use Basato_module
  Use vdwd3_module
  use gcp_module
  Use Symbat_module
  USE MEMORY_USE
  Implicit Real(float)(a-h,o-z)
  PARAMETER(NDIMARG=20,NMETHOD=81)
  character(len=9),PARAMETER :: znamz="GCP_INPUT"
  character(len=100) :: vdwmode,arglist(NDIMARG),arg,args(NDIMARG),funclist(NDIMARG)
  logical :: knownmethod,parm,ex
  character(len=25) :: gcpmethoddum, gcpfunc,tmp,gcp_methodlist(NMETHOD)
  integer,parameter :: MAX_ELEM=94
 REAL(FLOAT),DIMENSION(:),ALLOCATABLE :: EMISS
 INTEGER,DIMENSION(:),ALLOCATABLE :: NBAS
  DIMENSION p(4),xx(10)
  character(LEN=80) :: atmp,ftmp,dtmp,basname
  character(LEN=12) :: scalestring
  DATA gcp_methodlist/ &
  "hf/sv","hf/svp","hf/631gd","hf/minis","hf/minix","dft/lanl",        &
  "b3lyp/lanl","dft/sv","b3lyp/sv","dft/sv(p)","b3lyp/sv(p)","dft/svx",&
  "b3lyp/svx","dft/svp","b3lyp/svp","dft/631gd","b3lyp/631gd",         &
  "dft/minis","b3lyp/minis","dft/tz","b3lyp/tz","dft/pobtz",           &
  "b3lyp/pobtz","blyp/minis","gga/minis","gga/sv","blyp/sv","gga/svp", &
  "blyp/svp","gga/tz","blyp/tz","tpss/minis","tpss/svp","pw6b95/minis",&
  "pw6b95/svp","hf/2g","hf/sv(p)",                                     &
  "hf/def2tzvp","hf/deftzvp","hf/tzvp","hf/ccdz", "hf/ccpvdz",         &
  "hf/accdz","hf/augccpvdz","hf/sv(p)","hf/def2sv(p)","hf/sv_p",       &
  "hf/def2svp","b3lyp/def2sv(p)","pbe/def2sv(p)","dft/def2sv(p)",      &
  "dft/sv_p","dft/def2tzvp","b3lyp/def2tzvp","dft/deftzvp",            &
  "b3lyp/deftzvp", "dft/tzvp", "b3lyp/tzvp","dft/ccdz","dft/ccpvdz",   &
  "b3lyp/ccdz","b3lyp/ccpvdz","dft/accdz","dft/augccpvdz",             &
  "b3lyp/accdz","b3lyp/augccpvdz","hf/dzp","dft/dzp","rsh2c/moddz",    &
  "pbeh3c","b3pbe3c","hse3c","b973c","pbesol03c","hsesol3c",           &
  "hfsol/minix","hf/631gdp","r2scansol3c","r2scan0sol3c","r2scanpob3c",&
  "r2scan0pob3c"/
  NAF=inf(24)
  if(.not.allocated(gcp_grad)) CALL CRYALLOC(GCP_GRAD,3,NAF,'GCP_INPUT','GCP_GRAD')
  if(.not.allocated(gcp_grad2))CALL CRYALLOC(GCP_GRAD2,NAF*3,'GCP_INPUT','GCP_GRAD2')
  knownmethod=.false.
  presentsigma=.false.
  presenteta=.false.
  presentalpha=.false.
  presentbeta=.false.
  gcpmethoddum=gcpmethod
  gcpradius=60._FLOAT
 CALL CRYALLOC(EMISS,MAX_ELEM,ZNAMZ,'EMISS')
 CALL CRYALLOC(NBAS,MAX_ELEM,ZNAMZ,'NBAS')
if(inf(308).eq.2) then      
  ! Get gcp arguments
  II=1
  do
    read(iin,'(a)')arg
    if(arg.eq.'END')exit
    if(II.gt.40)exit
    call JGB_parse(arg,' ',args,nargs)
    do k=1,nargs
       arglist(II)=trim(args(k))
      II=II+1 
    enddo
  enddo
  do i=1,II
    call convert_lu(arglist(i), arglist(i))
    if(arglist(i).ne.'') then
    if(arglist(i).eq.'METHOD') then
      gcpmethoddum=trim(adjustl(arglist(i+1)))
      presentmethod=.true.
    elseif(arglist(i).eq.'SIGMA') then
      call JGB_value(arglist(i+1),gcpsigma,ios)
      presentsigma=.true.
    elseif(arglist(i).eq.'ETA') then
      call JGB_value(arglist(i+1),gcpeta,ios)
      presenteta=.true.
    elseif(arglist(i).eq.'ALPHA') then
      call JGB_value(arglist(i+1),gcpalpha,ios)
      presentalpha=.true.
    elseif(arglist(i).eq.'BETA') then
      call JGB_value(arglist(i+1),gcpbeta,ios)
      presentbeta=.true.
    elseif(arglist(i).eq.'RADIUS') then
      call JGB_value(arglist(i+1),gcpradius,ios)
      presentradius=.true.
    elseif(arglist(i).eq.'SCALE') then
      call JGB_value(arglist(i+1),gcpscale,ios)
      presentscale=.true.
   elseif(arglist(i).eq.'PRINTEMISS') then
      lprint(128)=1
    endif
    endif  
  enddo

  call errvrs(2,znamz,'PERFORM GEOMETRIC BSSE CORRECTION GCP')
  
elseif(inf(308).eq.5) then
  inf(308)=2
  call errvrs(2,znamz,'WITH AUTOMATIC PARAMETER SETUP')
endif

if(inf(308).eq.2)then
  WRITE(IOUT,10)
  if(inf(208).ne.2) then
  call errvrs(1,znamz,'GCP SHOULD BE ALWAYS COMBINED WITH AN ACCURATE DISPERSION CORRECTION ')  
  call errvrs(1,znamz,'CONSIDER USING THE DFTD3 METHOD')
  endif
endif


if(inf(170).eq.0)vdwfunc='hf'

gcpfunc=lowercase(vdwfunc)
gcpfunc=trim(gcpfunc)
gcpbasis=lowercase(gcpbasis)
gcpbasis=trim(gcpbasis)


do while (index(gcpfunc,'-').ne.0)
  i=scan(gcpfunc,'-')
  tmp=trim(gcpfunc(:(i-1)))//trim(gcpfunc((i+1):))
  gcpfunc=trim(tmp)
enddo
do while (index(gcpbasis,'-').ne.0)
  i=scan(gcpbasis,'-')
  tmp=trim(gcpbasis(:(i-1)))//trim(gcpbasis((i+1):))
  gcpbasis=trim(tmp)
enddo
if((gcpbasis.eq.'pobtzvp').or.(gcpbasis.eq.'tzvp')) then
  gcpbasis='pobtz'
elseif((gcpbasis.eq.'pobdzvp').or.(gcpbasis.eq.'dzvp')) then
  gcpbasis='pobdz'
elseif(gcpbasis.eq.'def2svp') then
  gcpbasis='svp'
elseif(gcpbasis.eq.'def2sv(p)') then
  gcpbasis='sv(p)'
elseif(gcpbasis.eq.'def2tzvp') then
  gcpbasis='tz'
elseif(gcpbasis.eq.'def2tzvpp') then
  gcpbasis='tz'
elseif(gcpbasis.eq.'631g**') then
  gcpbasis='631gd'
else
  if(.not.presentmethod) then
    call errvrs(1,znamz,'gCP BASIS NOT KNOWN')
    call errvrs(1,znamz,'PROVIDE GCPMETHOD IN INPUT BLOCK')
    call errvrs(1,znamz,'PARAMETRIZED METHODS')
    do i=1,NMETHOD
      write(IOUT,'(a20,5x)')  gcp_methodlist(i)
      if(mod(i,5).eq.0) write(IOUT,'('' '')')
    enddo
    write(IOUT,'('' '')')  
    call errvrs(0,znamz,'ERROR IN GCP INPUT')
  endif
endif


gcpmethod=trim(adjustl(gcpfunc)) // '/'
gcpmethod=trim(adjustl(gcpmethod)) // trim(adjustl(gcpbasis))
gcpmethod=trim(adjustl(gcpmethod))

call lower_case(gcpmethod)
do while (index(gcpmethod,'-').ne.0)
  i=scan(gcpmethod,'-')
  gcpmethod=trim(gcpmethod(:(i-1)))//trim(gcpmethod((i+1):))
  gcpmethod=trim(gcpmethod)
enddo

call lower_case(gcpmethoddum)
do while (index(gcpmethoddum,'-').ne.0)
  i=scan(gcpmethoddum,'-')
  gcpmethoddum=trim(gcpmethoddum(:(i-1)))//trim(gcpmethoddum((i+1):))
  gcpmethoddum=trim(gcpmethoddum)
enddo

if(presentmethod) then
  if(gcpmethod.ne.gcpmethoddum) then
    call errvrs(1,znamz,'FUNCTIONAL/BASIS COMBINATION AND GCPMETHOD NOT IDENTICAL')
    call errvrs(1,znamz,'TAKE GCPMETHOD FOR PARAMETER SETUP')
    gcpmethod=lowercase(gcpmethoddum)
  endif
endif
  
  do i=1,NMETHOD
    if(lowercase(gcpmethod).eq.lowercase(gcp_methodlist(i)))knownmethod=.true.
  enddo

if((.not.knownmethod).and.(.not.presentmethod)) then
  if((inf(170).ne.0).and.(PAR(48).gt.0._FLOAT).and.(PAR(48).lt.100._FLOAT)) then
    gcpmethod='dft/'
    gcpmethod=trim(adjustl(gcpmethod)) // trim(adjustl(gcpbasis))
    gcpmethod=trim(adjustl(gcpmethod))
    do i=1,NMETHOD
      if(gcpmethod.eq.gcp_methodlist(i))knownmethod=.true.
    enddo
  elseif(inf(170).ne.0) then
    gcpmethod='gga/'
    gcpmethod=trim(adjustl(gcpmethod)) // trim(adjustl(gcpbasis))
    gcpmethod=trim(adjustl(gcpmethod))
    do i=1,NMETHOD
      if(gcpmethod.eq.gcp_methodlist(i))knownmethod=.true.
    enddo
  endif
endif


if(gcpmethod.eq.'b973c') then
   call errvrs(2,znamz,'SWITCH TO SRB CORRECTION FOR B97-3c')
else   
if(knownmethod) then
    call JGB_setparam(emiss,nbas,p,gcpmethod)
else
  if(presentsigma.and.presenteta.and.presentalpha.and.presentbeta) then
    call JGB_setparam(emiss,nbas,p,gcpmethod)
  else
    call JGB_setparam(emiss,nbas,p,gcpmethod)
    call errvrs(1,znamz,'FUNCTIONAL NOT PARAMETRIZED')
    call errvrs(1,znamz,'PROVIDE GCP PARAMETERS IN INPUT BLOCK')
    call errvrs(1,znamz,'PARAMETRIZED METHODS')
    do i=1,NMETHOD
      write(IOUT,'(a20,5x)')  gcp_methodlist(i)
      if(mod(i,5).eq.0) write(IOUT,'('' '')')
    enddo
    write(IOUT,'('' '')')  
    call errvrs(0,znamz,'ERROR IN GCP INPUT')
  end if
end if
endif

if(presentsigma) then
  gcpsigma=gcpsigma 
  call errvrs(1,znamz,'USE MANUALLY DEFINED SIGMA')
else
  gcpsigma=p(1)
endif  
if(presenteta) then
  gcpeta=gcpeta
  call errvrs(1,znamz,'USE MANUALLY DEFINED ETA')
else
  gcpeta=p(2)
endif  
if(presentalpha) then
  gcpalpha=gcpalpha
  call errvrs(1,znamz,'USE MANUALLY DEFINED ALPHA')
else
  gcpalpha=p(3)
endif
if(presentbeta) then
  gcpbeta=gcpbeta
  call errvrs(1,znamz,'USE MANUALLY DEFINED BETA')
else
  gcpbeta=p(4)
endif

if(presentscale) then
  call writenum(gcpscale,scalestring,'f5.3')
  call errvrs(1,znamz,'SCALE GCP BY FACTOR '//trim(scalestring))
else
  gcpscale=1._FLOAT
endif
 CALL CRYDEALLOC(NBAS,ZNAMZ,'NBAS')
 CALL CRYDEALLOC(EMISS,ZNAMZ,'EMISS')
 RETURN
10 FORMAT(/' PLEASE CITE'/' (1) GCP Reference:'/'  H. Kruse and S. Grimme J. Chem. Phys, 2012,136,154101'/ &
   ' (2) GCP for Periodic Systems'/'  J. G. Brandenburg, M. Alessio, B. Civalleri, M. F. Peintinger,', &
   ' T. Bredow, and S.Grimme,'/'  J. Phys. Chem. A, 2013,117,9282'/)
end Subroutine gcp_input

subroutine gcp_free
Use gcp_module
Use Symbat_module
USE MEMORY_USE
  CALL CRYDEALLOC(GCP_GRAD2,'GCP_INPUT','GCP_GRAD2')
  CALL CRYDEALLOC(GCP_GRAD,'GCP_INPUT','GCP_GRAD')
endsubroutine gcp_free

subroutine gcp_calc_energy(echo)
  use Numbers
  use vdwd3_module
  use strings
  Use Parinf_module
  Use Basato_module
  Use Gvect_module
  USE XYVDIM_MODULE
  use gcp_module
  USE MEMORY_USE
  Implicit Real(float)(a-h,o-z)
  character(len=12),PARAMETER :: znamz="GCP  " 
  logical :: echo
  integer,DIMENSION(:),ALLOCATABLE :: iz
  NAF=inf(24)
  CALL CRYALLOC(IZ,NAF,'GCP_CALC_ENERGY','IZ')
  do i=1,3
    do j=1,3
      cell(i,j) = PARET(j,i)
    enddo
  enddo

  IZ(1:NAF)=MOD(NAT(1:NAF),100)
  if(echo)WRITE(IOUT,'(/'' GCP ENERGY CALCULATION'')')
  
  call pbcgcp(xa,iz,cell,NAF,gcp_energy,gcp_grad,gcp_cellgrad,  &
     .false.,gcpmethod,echo,gcpsigma,gcpeta,gcpalpha,gcpbeta)
  gcp_energy=gcp_energy*gcpscale
  gcp_grad=gcp_grad*gcpscale
  gcp_cellgrad=gcp_cellgrad*gcpscale

  if(allocated(iz)) CALL CRYDEALLOC(IZ,'GCP_CALC_ENERGY','IZ')
  if(echo)call errvrs(2,znamz,'NORMAL TERMINATION OF GCP MODULE')
  
end subroutine gcp_calc_energy

subroutine gcp_calc_energy_grad
  use Numbers
  use vdwd3_module
  use strings
  Use Parame_module
  Use Parinf_module
  Use Basato_module
  Use Gvect_module
  USE XYVDIM_MODULE
  use gcp_module
  USE MEMORY_USE
  Implicit Real(float)(a-h,o-z)
  character(len=12),PARAMETER :: znamz="GCP  " 
  !number of atoms
  logical echo
  integer,DIMENSION(:),ALLOCATABLE :: iz
  Real(float) :: gcp_cellgrad_dum(3,3)
  NAF=inf(24)
  ldim=inf(10)
  if(lprint(128).eq.1) then
     echo=.true.
  else
     echo=.false.
  end if

  CALL CRYALLOC(IZ,NAF,'GCP_CALC_ENERGY_GRAD','IZ')
  do i=1,3
    do j=1,3
      cell(i,j) = PARET(j,i)
    enddo
  enddo

  IZ(1:NAF)=MOD(NAT(1:NAF),100)
  if(echo) then
   call errvrs(2,znamz,' ')
   call errvrs(2,znamz,'CALCULATE GCP ENERGY AND GCP FORCE')
   WRITE(IOUT,'(/'' GCP ENERGY AND FORCE CALCULATION'')')
  endif

   call pbcgcp(xa,iz,cell,NAF,gcp_energy,gcp_grad,gcp_cellgrad,  &
        .true.,gcpmethod,.false.,gcpsigma,gcpeta,gcpalpha,gcpbeta)
 
   
   gcp_energy=gcp_energy*gcpscale
   gcp_grad=gcp_grad*gcpscale
   gcp_cellgrad=gcp_cellgrad*gcpscale


  k=1
  do i=1,NAF
    do j=1,3
      gcp_grad2(k)=-gcp_grad(j,i)
      k=k+1
    enddo
  enddo  
  
  gcp_cellgrad_dum=gcp_cellgrad
  do i=1,3
     do j=1,3
        gcp_cellgrad(i,j)=gcp_cellgrad_dum(j,i)
     enddo
  enddo  
  
   gcp_cellgrad2(1)=gcp_cellgrad(1,1)
  select case (ldim)
  case(2)
   gcp_cellgrad2(2)=gcp_cellgrad(1,2)
   gcp_cellgrad2(3)=gcp_cellgrad(2,1)
   gcp_cellgrad2(4)=gcp_cellgrad(2,2)
  case(3)
   gcp_cellgrad2(2)=gcp_cellgrad(1,2)
   gcp_cellgrad2(3)=gcp_cellgrad(1,3)
   gcp_cellgrad2(4)=gcp_cellgrad(2,1)
   gcp_cellgrad2(5)=gcp_cellgrad(2,2)
   gcp_cellgrad2(6)=gcp_cellgrad(2,3)
   gcp_cellgrad2(7)=gcp_cellgrad(3,1)
   gcp_cellgrad2(8)=gcp_cellgrad(3,2)
   gcp_cellgrad2(9)=gcp_cellgrad(3,3)
  end select  

  
  CALL CRYDEALLOC(IZ,'GCP_CALC_ENERGY','IZ')
  if(echo)  call errvrs(2,znamz,'NORMAL TERMINATION OF GCP MODULE')
end subroutine gcp_calc_energy_grad
subroutine pbcgcp(xyz,iz,Hlat,n,gcp,g,glat,grad,method,echo,sigma,eta,alpha,beta)
 USE NUMBERS
 USE PARINF_MODULE
use strings
  USE MEMORY_USE
 use gcp_module
 IMPLICIT DOUBLE PRECISION(A-H,O-Z)
parameter (autokcal=627.509541_FLOAT,max_elem=94,max_para=36)
 DIMENSION xyz(3,n),g(3,n),HLAT(3,3),GLAT(3,3),P(4),iz(n)
 REAL(FLOAT),DIMENSION(:),ALLOCATABLE :: EMISS,XVA,XVB
 INTEGER,DIMENSION(:),ALLOCATABLE :: IZ2,IFREZ,NEL,NBAS
logical echo,grad,warn,damp,srb
character(LEN=20) :: method,ctmp
character(LEN=2) :: esym
character(LEN=80) :: getlevel
 AUTOANG=PAR(32)
p(1)=sigma
p(2)=eta
p(3)=alpha
p(4)=beta
warn=.false.
damp=.false.
srb=.false.
 CALL CRYALLOC(EMISS,MAX_ELEM,'PBCGCP','EMISS')
 CALL CRYALLOC(NBAS,MAX_ELEM,'PBCGCP','NBAS')
 CALL CRYALLOC(NEL,MAX_ELEM,'PBCGCP','NEL')
 CALL CRYALLOC(IZ2,N,'PBCGCP','IZ2')
call lower_case(method)
do while (index(method,'-').ne.0) 
  i=scan(method,'-')
  ctmp=trim(method(:(i-1)))//trim(method((i+1):))
  method=trim(ctmp)
enddo
select case(method)
case('pbeh3c')
damp=.true.
gcp_dmp_scal=4.0_FLOAT
gcp_dmp_exp=6.0_FLOAT
case('hse3c')
damp=.true.
gcp_dmp_scal=4.0_FLOAT
gcp_dmp_exp=6.0_FLOAT
case('b973c')
srb=.true.
rscal=10.00_FLOAT
qscal=0.08_FLOAT
case('pbesol03c')
damp=.true.
gcp_dmp_scal=4.0_FLOAT
gcp_dmp_exp=6.0_FLOAT
case('hsesol3c')
damp=.true.
gcp_dmp_scal=4.0_FLOAT
gcp_dmp_exp=6.0_FLOAT
case('r2scansol3c')
damp=.true.
gcp_dmp_scal=4.0_FLOAT
gcp_dmp_exp=6.0_FLOAT
case('r2scan0sol3c')
damp=.true.
gcp_dmp_scal=4.0_FLOAT
gcp_dmp_exp=6.0_FLOAT
case('r2scanpob3c')
damp=.true.
gcp_dmp_scal=4.0_FLOAT
gcp_dmp_exp=6.0_FLOAT
case('r2scan0pob3c')
damp=.true.
gcp_dmp_scal=4.0_FLOAT
gcp_dmp_exp=6.0_FLOAT
case default
end select
if(.not.srb) then
   call JGB_setparam(emiss,nbas,p,method)
   p(1)=sigma
   p(2)=eta
   p(3)=alpha
   p(4)=beta   
iz2=iz
do i=1,n
izz=iz(i)
select case (izz)
 case(37:54)
  iz(i)=iz(i)-18
if(echo) WRITE(IOUT,10)esym(izz),esym(iz(i))
  warn=.true.
 case(55:57)
  iz(i)=iz(i)-18*2
if(echo)  WRITE(IOUT,10)esym(izz),esym(iz(i))
  warn=.true.
 case(58:71,90:94)
  iz(i)=21
if(echo)  WRITE(IOUT,10)esym(izz),esym(iz(i))
  warn=.true.
 case(72:89)
  iz(i)=iz(i)-2*18-14
if(echo)  WRITE(IOUT,10)esym(izz),esym(iz(i))
  warn=.true.
end select 
enddo
do i=1,maxval(iz)
  nel(i)=i
enddo
if(index(method,'vmb').ne.0) then
 CALL ERRVRS(2,'PBCGCP','VMB basis has ECPs. Adjusting electrons')
do i=5,10
nel(i)=nel(i)-21
enddo
do i=11,18
nel(i)=nel(i)-10
enddo
endif
nb=0
do nn=1,n
   nb=nb+nbas(iz(nn))
enddo
if(nb.lt.1)CALL ERRVRS(0,'PBCGCP','Nbf setup gone wrong')
 CALL CRYALLOC(XVA,N,'PBCGCP','XVA')
 CALL CRYALLOC(XVB,N,'PBCGCP','XVB')
 CALL CRYALLOC(IFREZ,N,'PBCGCP','IFREZ')
 xva=0._FLOAT
do i=1,n
  xva(i)=(nbas(iz(i))-nel(iz(i))*.5_FLOAT)  
  if(emiss(iz(i)).gt.1e-7_FLOAT.and.xva(i).lt.0._FLOAT)then
  WRITE(IOUT,11)esym(iz(i)),emiss(iz(i)),xva(i),nel(iz(i)),nbas(iz(i))
   CALL ERRVRS(0,'PBCGCP','negative number of virtual orbitals. Something is wrong in the parameters!')
  endif
  if(echo)then
   if(xva(i).lt.0.5_FLOAT)then
   WRITE(IOUT,12)esym(iz(i))
   CALL ERRVRS(1,'PBCGCP','ELEMENT has no virtual orbitals (no contribution)!')
   warn=.true.
  endif
  endif
enddo
xvb=xva
endif
if(echo)then
   if(method.eq."#".or.method.eq."hf/2g")then
      write(IOUT,'(2x,''level '',3x,a12)')adjustr(trim(method))
   else
      write(IOUT,'(2x,''level '',3x,a12,2x,a24)')adjustr(trim(method)),adjustr(trim(getlevel(method)))
   endif
   write(IOUT,'(2x,''Nbf   '',I12)')nb
   write(IOUT,'(2x,''Atoms '',I12)')n
   if(.not.srb) then
      write(IOUT,'(/2x,a12)')'Parameters: '
      write(IOUT,'(2x,a6,f10.4)')'sigma ',p(1)
      write(IOUT,'(2x,a6,f10.4)')'eta   ',p(2)
      write(IOUT,'(2x,a6,f10.4)')'alpha ',p(3)
      write(IOUT,'(2x,a6,f10.4/)')'beta  ',p(4)
      if(damp) then
         dmp_scal=gcp_dmp_scal
         dmp_exp=gcp_dmp_exp
         write(IOUT,'(2x,a,2x,f10.4)') 'dmp_scal',dmp_scal
         write(IOUT,'(2x,a,2x,f10.4)') 'dmp_exp ',dmp_exp
      endif
   else
      write(IOUT,'(2x,''Atoms '',I12)')n
      write(IOUT,'(/2x,a22)')'Parameters for SRB: '
      write(IOUT,'(2x,a,2x,f10.4)') 'rscal',rscal
      write(IOUT,'(2x,a,2x,f10.4)') 'qscal',qscal
   endif
endif

      g(1:3,1:n)=0._FLOAT
if(srb) then
   call srb_egrad2(xyz,iz,Hlat,n,gcp,g,glat,grad,rscal,qscal,echo)
else
      call JGB_e(n,max_elem,emiss,xyz,iz,p,gcp,g,glat,grad,echo,xva,xvb,Hlat,damp)
endif
iz=iz2


    if(echo)then
      if(srb) then
      write(IOUT,*)'** gCP correction **'
      write(IOUT,'(2x,a7,F18.10,'' / (a.u.) || '',1x,F11.4,'' / (kcal/mol)'')')'Egcp:  ',gcp,gcp*AUTOKCAL
       if(grad)write(IOUT,*)'|G|=',sum(abs(g(1:3,1:n)))
      else
         write(IOUT,*)'** gCP correction **'
         write(IOUT,'(2x,a7,F18.10,'' / (a.u.) || '',1x,F11.4,'' / (kcal/mol)'')')'Egcp:  ',gcp,gcp*AUTOKCAL
         if(grad)write(IOUT,*)'|G|=',sum(abs(g(1:3,1:n)))
      endif
   endif

    if(echo) then
    if(warn)CALL ERRVRS(1,'PBCGCP','Carefully read the notifications/warnings given at loadup')
    endif
if(allocated(IFREZ)) CALL CRYDEALLOC(IFREZ,'PBCGCP','IFREZ')
if(allocated(XVB)) CALL CRYDEALLOC(XVB,'PBCGCP','XVB')
if(allocated(XVA)) CALL CRYDEALLOC(XVA,'PBCGCP','XVA')
if(allocated(IZ2)) CALL CRYDEALLOC(IZ2,'PBCGCP','IZ2')
if(allocated(NEL)) CALL CRYDEALLOC(NEL,'PBCGCP','NEL')
if(allocated(NBAS)) CALL CRYDEALLOC(NBAS,'PBCGCP','NBAS')
if(allocated(EMISS)) CALL CRYDEALLOC(EMISS,'PBCGCP','EMISS')
 RETURN
10 FORMAT('NOTE: ELEMENT ',A2,' WILL BE TREATED AS ',A2)
11 FORMAT(A2,1P,2E16.6,2I8)
12 FORMAT('ELEMENT ',A2)
 end subroutine pbcgcp


subroutine JGB_E(n,max_elem,emiss,xyz,iz,p,gcp,g,cellgrad,grad,echo,xva,xvb,Hlat,damp)
 USE NUMBERS
 USE PARINF_MODULE
 USE PARAL1_MODULE
 use gcp_module 
 USE MEMORY_USE
 IMPLICIT DOUBLE PRECISION(A-H,O-Z)
 PARAMETER(autokcal=627.509541_FLOAT,max_para=36)
 REAL(FLOAT),DIMENSION(:,:),ALLOCATABLE :: XYZJAT
 REAL(FLOAT),DIMENSION(:),ALLOCATABLE :: ZA,ZB
 DIMENSION XYZ(3,N),G(3,N),EMISS(MAX_ELEM),XVA(N),XVB(N),P(4),GS(3)
 DIMENSION R0AB(MAX_ELEM,MAX_ELEM),Hlat(3,3),IZ(N),ITAU_MAX(3),Hlat_1(3,3)
 real(FLOAT) :: cellgrad(3,3), stress(3,3)
logical echo,grad,damp
thrR=gcpradius            ! X bohr
thrE=epsilon(1._FLOAT)
 cellgrad=0.0_FLOAT
 stress=0.0_FLOAT
if(echo)then
write(IOUT,'('' cutoffs: ''/''   distance [bohr]'',F5.1/''   noise    [a.u.]'',Es8.1)')thrR,thrE
write(IOUT,'('' SR-damping  '',L2)')damp
if(lprint(128).eq.1) then
write(IOUT,'(/2x,a5,2x,a7,2x,a5,4x,a15)') &
        '#','Nvirt','Emiss','BSSE [kcal/mol]'     
endif
endif
if(damp) then
  call setr0ab(max_elem,r0ab)
  dmp_scal=gcp_dmp_scal
  dmp_exp=gcp_dmp_exp
endif
g=0._FLOAT
ecp=0._FLOAT
tg=0._FLOAT
p1=abs(p(1))
p2=abs(p(2))
p3=abs(p(3))
p4=abs(p(4))
 CALL CRYALLOC(ZA,MAX_PARA,'JGB_E','ZA')
 CALL CRYALLOC(ZB,MAX_PARA,'JGB_E','ZB')
call setzet(p2,za,zb)
 CALL CRYALLOC(XYZJAT,3,N,'JGB_E','XYZJAT')
call JGB_criteria(thrR, Hlat,Itau_max)
do iat=1,n
   va=xva(iat)
   ea=0._FLOAT
   np=0
   ITAU1=Itau_max(1)
   FATI=-ITAU1
   ITAU2=Itau_max(2)
   FATJ0=-ITAU2
   ITAU3=Itau_max(3)
   FATK0=-ITAU3
   do i=-Itau1,Itau1
   FATJ=FATJ0
      do j=-Itau2,Itau2
   FATK=FATK0
         do k=-Itau3,Itau3
            do jat=1,n
               iselftest=0
               ! Test for equal atoms, remove selfinteraction
               if(iat.eq.jat) then
                  iselftest=iselftest+1
                  if(i.eq.0)iselftest=iselftest+1
                  if(j.eq.0)iselftest=iselftest+1
                  if(k.eq.0)iselftest=iselftest+1
               end if
               if(iselftest.eq.4)cycle
               dx=xyz(1,iat)-xyz(1,jat)+i*Hlat(1,1)+j*Hlat(1,2)+k*Hlat(1,3)
               dy=xyz(2,iat)-xyz(2,jat)+i*Hlat(2,1)+j*Hlat(2,2)+k*Hlat(2,3)
               dz=xyz(3,iat)-xyz(3,jat)+i*Hlat(3,1)+j*Hlat(3,2)+k*Hlat(3,3)

               r=sqrt(dx*dx+dy*dy+dz*dz)
               vb=xvb(jat)
               if(vb.lt.0.5_FLOAT) cycle
               if(r.gt.thrR) cycle
               call ssovl(r,iat,jat,iz,za(iz(iat)),zb(iz(jat)),sab)
               if(abs(sab).lt.thrE) cycle
               ene_old_num=exp(-p3*r**p4)
               ene_old_den=sqrt(vb*Sab)
               ene_old=ene_old_num/ene_old_den
               if(abs(ene_old).lt.thrE) cycle
               if(damp) then
                r0abij=r0ab(iz(iat),iz(jat))
                rscal=r/r0abij
                rscalexp=rscal**dmp_exp
                dampval=(1.0_FLOAT-1.0_FLOAT/(1.0_FLOAT+dmp_scal*rscalexp))
                ene_dmp=ene_old*dampval
                ea=ea+emiss(iz(iat))*ene_dmp
              else
                ea=ea+emiss(iz(iat))*ene_old
              endif
               np=np+1
               if(grad)then
                  call cpu_time(t0)
                  gab=0.0_FLOAT
                  call gsovl(r,iat,jat,iz,za(iz(iat)),zb(iz(jat)),gab)
                  gs(1)=gab*dx
                  gs(2)=gab*dy
                  gs(3)=gab*dz

                  dum=-exp(-p3*r**p4)*.5_FLOAT
                  dum2=p3*p4*r**p4*sab/r*2._FLOAT
                  dum22=r*sab**1.5_FLOAT
                  tmpb=dum22*sqrt(vb)
                  if(damp) then
                    rscalexpm1=rscal**(dmp_exp-1.0_FLOAT)
                    grd_dmp=dmp_scal*dmp_exp*rscalexpm1/r0abij
                    grd_dmp=grd_dmp/((dmp_scal*rscalexp+1.0_FLOAT)**2)
                  endif
                  tmpa=dum2*dx+gs(1)
                  tmp=dum*tmpa/tmpb
                  if(damp)tmp=tmp*dampval+ene_old*grd_dmp*(dx/r)
                  g(1,iat)=g(1,iat)+tmp*emiss(iz(iat))
                  stress(1,1)=stress(1,1)+dx* tmp*emiss(iz(iat))
                  stress(2,1)=stress(2,1)+dy* tmp*emiss(iz(iat))
                  stress(3,1)=stress(3,1)+dz* tmp*emiss(iz(iat))
            
                  tmpa=dum2*dy+gs(2)
                  tmp=dum*tmpa/tmpb
                  if(damp)tmp=tmp*dampval+ene_old*grd_dmp*(dy/r)
                  g(2,iat)=g(2,iat)+tmp*emiss(iz(iat))
                  stress(1,2)=stress(1,2)+dx* tmp*emiss(iz(iat))
                  stress(2,2)=stress(2,2)+dy* tmp*emiss(iz(iat))
                  stress(3,2)=stress(3,2)+dz* tmp*emiss(iz(iat))

                  tmpa=dum2*dz+gs(3)
                  tmp=dum*tmpa/tmpb
                  if(damp)tmp=tmp*dampval+ene_old*grd_dmp*(dz/r)
                  g(3,iat)=g(3,iat)+tmp*emiss(iz(iat))
                  stress(1,3)=stress(1,3)+dx* tmp*emiss(iz(iat))
                  stress(2,3)=stress(2,3)+dy* tmp*emiss(iz(iat))
                  stress(3,3)=stress(3,3)+dz* tmp*emiss(iz(iat))

                  if(va.lt.0.5_FLOAT) cycle
                  if(damp) then
                    ene_old_den=sqrt(va*Sab)
                    ene_old=ene_old_num/ene_old_den
                  endif
                  tmpb=dum22*sqrt(va)

                  tmpa=-(dum2*dx+gs(1))
                  tmp=dum*tmpa/tmpb
                  if(damp)tmp=tmp*dampval+ene_old*grd_dmp*(-dx/r)
                  g(1,iat)=g(1,iat)-tmp*emiss(iz(jat))
                  stress(1,1)=stress(1,1)-dx* tmp*emiss(iz(jat))
                  stress(2,1)=stress(2,1)-dy* tmp*emiss(iz(jat))
                  stress(3,1)=stress(3,1)-dz* tmp*emiss(iz(jat))

                  tmpa=-(dum2*dy+gs(2))
                  tmp=dum*tmpa/tmpb
                  if(damp)tmp=tmp*dampval+ene_old*grd_dmp*(-dy/r)
                  g(2,iat)=g(2,iat)-tmp*emiss(iz(jat))
                  stress(1,2)=stress(1,2)-dx* tmp*emiss(iz(jat))
                  stress(2,2)=stress(2,2)-dy* tmp*emiss(iz(jat))
                  stress(3,2)=stress(3,2)-dz* tmp*emiss(iz(jat))


                  tmpa=-(dum2*dz+gs(3))
                  tmp=dum*tmpa/tmpb
                  if(damp)tmp=tmp*dampval+ene_old*grd_dmp*(-dz/r)
                  g(3,iat)=g(3,iat)-tmp*emiss(iz(jat))
                  stress(1,3)=stress(1,3)-dx* tmp*emiss(iz(jat))
                  stress(2,3)=stress(2,3)-dy* tmp*emiss(iz(jat))
                  stress(3,3)=stress(3,3)-dz* tmp*emiss(iz(jat))

                  
                  tg=tg-t0
               endif
            enddo
         FATK=FATK+1._FLOAT
         enddo
         FATJ=FATJ+1._FLOAT
      end do
         FATI=FATI+1._FLOAT
   end do
   if(lprint(128).eq.1) then
   if(echo)write(IOUT,'(2x,I5,2x,F5.1,2x,F10.4,2x,F10.4,2x,a)')  &
        iz(iat),va,emiss(iz(iat)),ea*p(1)*AUTOKCAL
   endif
   ecp=ecp+ea
enddo

call MINV3(Hlat,Hlat_1,Hdet)
 cellgrad=0.0_FLOAT
 CALL MXMBN(stress,1,3,Hlat_1,3,1,cellgrad,1,3,3,3,3)
 cellgrad=0.5_FLOAT*cellgrad

 CALL CRYDEALLOC(XYZJAT,'JGB_E','XYZJAT')
 CALL CRYDEALLOC(ZB,'JGB_E','ZB')
 CALL CRYDEALLOC(ZA,'JGB_E','ZA')
 stress=stress*p(1)
 cellgrad=cellgrad*p(1)
 gcp=ecp*p(1)
 g=g*p(1)
 RETURN
 end

subroutine gcp_numcellgrad2
use parame_module
use parinf_module
use basato_module
USE PARAL1_MODULE
use gradient_memory
USE NUMDFT_MODULE
USE ORDER_DFT_BATCHES
use MEMORY_G
use memory_opt
use xyvdim_module
use retic_module
use gvect_module
use symmetd_module
use saed_module
use mp2ene
use gcp_module
USE MEMORY_USE
IMPLICIT REAL(FLOAT) (A-H,O-Z)
character(len=12)::znamz='NUMEGRAD_GEN'
REAL(FLOAT),DIMENSION(:,:),ALLOCATABLE :: xg_r,trasv_r
REAL(FLOAT),DIMENSION(:),ALLOCATABLE :: gmodus_r,CODISP
INTEGER,DIMENSION(:),ALLOCATABLE :: nm_r
DIMENSION ECB(2),paret_r(3,3)
inf5=inf(5)
inf79=inf(79)
mvf=inf(2)
NAF=INF(24)
inf197=inf(197)
inf(197)=0
ncart=nsadir+nsaed
icz1=nsadir+1
icz2=ncart
CALL CRYALLOC(CODISP,ncart,ZNAMZ,'CODISP')
CALL CRYALLOC(xg_r,3,inf79,ZNAMZ,'xg_r')
CALL CRYALLOC(gmodus_r,inf5,ZNAMZ,'gmodus_r')
CALL CRYALLOC(nm_r,inf5,ZNAMZ,'nm_r')
xg_r(1:3,1:inf79)=xg(1:3,1:inf79)
gmodus_r(1:inf5)=gmodus(1:inf5)
nm_r(1:inf5)=nm(1:inf5)
ECB0=gcp_energy
disp=par(51)
DISP_INV=.5_float/disp
CODISP(1:ncart)=0._FLOAT
CALL CRYALLOC(TRASV_R,3,MVF,ZNAMZ,'TRASV_R')
trasv_r(1:3,1:mvf)=trasv(1:3,1:mvf)
paret_r(1:3,1:3)=paret(1:3,1:3)
icount=0
DO i=icz1,icz2
   icount=icount+1
   CODISP(I)=disp
   DO IST=1,2
      call gcoor_move(codisp)
!energy call
      call gcp_calc_energy(.false.)
      ECB(IST)=gcp_energy
      xa(1:3,1:naf)=xyza(1:3,1:naf)
      trasv(1:3,1:mvf)=trasv_r(1:3,1:mvf)
      paret(1:3,1:3)=paret_r(1:3,1:3)
      CODISP(I)=-disp
   ENDDO
   FORZE(I)=FORZE(I)+(ECB(2)-ECB(1))*DISP_INV
   ECB(1:2)=0._FLOAT
   CODISP(I)=0._FLOAT
ENDDO
CALL CRYDEALLOC(TRASV_R,ZNAMZ,'TRASV_R')
DO LA=1,INF(20)
   XL(1:3,LA)=XA(1:3,LATOAT(LA))
ENDDO
xg(1:3,1:inf79)=xg_r(1:3,1:inf79)
gmodus(1:inf5)=gmodus_r(1:inf5)
nm(1:inf5)=nm_r(1:inf5)
CALL CRYdealloc(nm_r,ZNAMZ,'NM_R')
CALL CRYdealloc(gmodus_r,ZNAMZ,'GMODUS_R')
CALL CRYdealloc(xg_r,ZNAMZ,'XG_R')
gcp_energy=ECB0
inf(197)=inf197
!LD
if(inf(408).ne.2.or.inf(408).ne.4) then
!LD
IF(printforces)WRITE(IOUT_current,103)(J,FORZE(J),J=icz1,icz2)
endif
CALL CRYDEALLOC(CODISP,ZNAMZ,'CODISP')
return
103   FORMAT(/' SYMMETRY ALLOWED FORCES (NUMERICAL)',&
     ' (DIRECTION, FORCE)'//(4(I5,1P,E15.7)))
end subroutine gcp_numcellgrad2


subroutine sgcp(n,max_elem,emiss,xyz,iz,p,grad,echo,xva,xvb,Hlat,sigma)
 USE NUMBERS
 USE PARINF_MODULE
  USE MEMORY_USE
 IMPLICIT DOUBLE PRECISION(A-H,O-Z)
 PARAMETER(DELTA=1E-5_FLOAT,DELTAM=-DELTA,DELTA2=5.E4_FLOAT)
 REAL(FLOAT),DIMENSION(:,:),ALLOCATABLE :: ABC,XYZDUM,GDUM
 DIMENSION XYZ(3,N),G(3,N),EMISS(MAX_ELEM),XVA(N),XVB(N),P(4),Hlat(3,3),HlatDUM(3,3),SIGMA(3,3),IZ(N)
 DIMENSION cellgrad(3,3)
logical echo,grad,damp !added clau 28/07/2016
integer a, b, c, Itau_max(3), abcsum, selftest, i, j,k
!cell volume
Hlatdum=0._FLOAT
detlat =Hlat(1,1)*Hlat(2,2)*Hlat(3,3)+Hlat(1,2)*Hlat(2,3)*Hlat(3,1)+Hlat(1,3)*Hlat(2,1)*Hlat(3,2)
detlat=detlat-Hlat(3,1)*Hlat(2,2)*Hlat(1,3)-Hlat(3,2)*Hlat(2,3)*Hlat(1,1)-Hlat(3,3)*Hlat(2,1)*Hlat(1,1)
Hlatdum=Hlat
 CALL CRYALLOC(ABC,3,N,'SGCP','ABC')
 CALL CRYALLOC(XYZDUM,3,N,'SGCP','XYZDUM')
 CALL CRYALLOC(GDUM,3,N,'SGCP','GDUM')
call xyz2abc(n,xyz,abc,Hlat)
do i=1,3
   do j=1,3
      Hlatdum(i,j)=Hlat(i,j)+deltaM
      call abc2xyz(n,xyzdum,abc,Hlatdum)
      call JGB_e(n,max_elem,emiss,xyzdum,iz,p,gcpl,gdum,cellgrad,.false.,.false.,xva,xvb,Hlatdum,damp)
      Hlatdum(i,j)=Hlat(i,j)+delta
      call abc2xyz(n,xyzdum,abc,Hlatdum)
      call JGB_e(n,max_elem,emiss,xyzdum,iz,p,gcpr,gdum,cellgrad,.false.,.false.,xva,xvb,Hlatdum,damp)
      ! symmetric numerical derivative
      sigma(i,j)=(gcpr-gcpl)*DELTA2
      !sigma(i,j)=sigma(i,j)/vol
   end do
end do
 CALL CRYDEALLOC(GDUM,'SGCP','GDUM')
 CALL CRYDEALLOC(XYZDUM,'SGCP','XYZDUM')
 CALL CRYDEALLOC(ABC,'SGCP','ABC')
 RETURN
end subroutine sgcp

subroutine JGB_setparam(emiss,nbas,p,method)
 USE NUMBERS
 USE PARINF_MODULE
 IMPLICIT DOUBLE PRECISION(A-H,O-Z)
 PARAMETER(MAX_ELEM=94,MAX_PARA=36)
character(len=12),PARAMETER :: znamz="JGB_SETPARAM"
 DIMENSION EMISS(MAX_ELEM),P(4),NBAS(MAX_ELEM)
 character(LEN=20) :: method
 REAL(FLOAT),DIMENSION(MAX_PARA) :: HFsv,HFminis,HF631gd,HFsvp,HFtz,HFvmb, &
 HFminisd,oldHFsvp,HFpobtz,HFpobdzvp,HF2g,HFdef1tzvp,HFccdz,HFaccdz,HFdzp, &
 HFmoddz,HFmsvp,HFdef2mtzvp,HFmsvp_ld,HFrevminix,HFpobtzrev2,HF631GDP
 REAL(FLOAT),DIMENSION(10) :: HFlanl2
 INTEGER,DIMENSION(MAX_PARA) :: BASsv,BASminis,BAS631gd,BAStz,BASsvp,BASvmb, &
 BASminisd,oldBASsvp,BASpobtz,BASpobdzvp,BAS2g,BASdef1tzvp,BASccdz,BASaccdz, &
 BASdzp,BASmoddz,BASmsvp,BASdef2mtzvp,BASrevminix,BASpobtzrev2,BAS631gdp
 INTEGER,DIMENSION(10) :: BASlanl2
data HFsv/ & !H-Kr (no 3d)
0.009037_FLOAT,0.008843_FLOAT,&  ! He_FLOAT,He
0.204189_FLOAT,0.107747_FLOAT,0.049530_FLOAT,0.055482_FLOAT,0.072823_FLOAT,0.100847_FLOAT,0.134029_FLOAT,0.174222_FLOAT,&  ! Li-Ne
0.315616_FLOAT,0.261123_FLOAT,0.168568_FLOAT,0.152287_FLOAT,0.146909_FLOAT,0.168248_FLOAT,0.187882_FLOAT,0.211160_FLOAT,&  !Na -Ar
0.374252_FLOAT,0.460972_FLOAT,&  ! K-Ca
0.444886_FLOAT,0.404993_FLOAT,0.378406_FLOAT,0.373439_FLOAT,0.361245_FLOAT,0.360014_FLOAT,0.362928_FLOAT,&
0.243801_FLOAT,0.405299_FLOAT,0.396510_FLOAT,&   ! 3d-TM
0.362671_FLOAT,0.360457_FLOAT,0.363355_FLOAT,0.384170_FLOAT,0.399698_FLOAT,0.417307_FLOAT/ !Ga-Kr
data HFminis/ &
0.042400_FLOAT,0.028324_FLOAT,&
0.252661_FLOAT,0.197201_FLOAT,0.224237_FLOAT,0.279950_FLOAT,0.357906_FLOAT,0.479012_FLOAT,0.638518_FLOAT,0.832349_FLOAT, &
1.232920_FLOAT,1.343390_FLOAT,1.448280_FLOAT,1.613360_FLOAT,1.768140_FLOAT,1.992010_FLOAT,2.233110_FLOAT,2.493230_FLOAT, &
3.029640_FLOAT,3.389980_FLOAT,&  ! H-Ca
10*0._FLOAT,6*0._FLOAT/
data HF631GD/ &! H-Ca + Br (no 3d)
0.010083_FLOAT,0.008147_FLOAT,&
0.069260_FLOAT,0.030540_FLOAT,0.032736_FLOAT,0.021407_FLOAT,0.024248_FLOAT,0.036746_FLOAT,0.052733_FLOAT,0.075120_FLOAT,& 
0.104255_FLOAT,0.070691_FLOAT,0.100260_FLOAT,0.072534_FLOAT,0.054099_FLOAT,0.056408_FLOAT,0.056025_FLOAT,0.057578_FLOAT,&
0.079198_FLOAT,0.161462_FLOAT,&
10*0.0_FLOAT, &
0.000000_FLOAT,0.000000_FLOAT,0.000000_FLOAT,0.000000_FLOAT,0.381049_FLOAT,0.000000_FLOAT/
data oldHFsvp / & ! Li_FLOAT,Na_FLOAT,Mg_FLOAT,K had wrong parameters
0.008107_FLOAT,0.008045_FLOAT,&
0.142957_FLOAT,0.028371_FLOAT,0.049369_FLOAT,0.055376_FLOAT,0.072785_FLOAT,0.100310_FLOAT,0.133273_FLOAT,0.173600_FLOAT,&
0.191109_FLOAT,0.222839_FLOAT,0.167188_FLOAT,0.149843_FLOAT,0.145396_FLOAT,0.164308_FLOAT,0.182990_FLOAT,0.205668_FLOAT,&
0.221189_FLOAT,0.299661_FLOAT,&
0.325995_FLOAT,0.305488_FLOAT,0.291723_FLOAT,0.293801_FLOAT,0.29179_FLOAT,0.296729_FLOAT,0.304603_FLOAT,0.242041_FLOAT,&
0.354186_FLOAT,0.350715_FLOAT,&
0.350021_FLOAT,0.345779_FLOAT,0.349532_FLOAT,0.367305_FLOAT,0.382008_FLOAT,0.399709_FLOAT/
data HFsvp /  & ! H-Kr
0.008107_FLOAT,0.008045_FLOAT,&
0.113583_FLOAT,0.028371_FLOAT,0.049369_FLOAT,0.055376_FLOAT,0.072785_FLOAT,0.100310_FLOAT,0.133273_FLOAT,0.173600_FLOAT,&
0.181140_FLOAT,0.125558_FLOAT,0.167188_FLOAT,0.149843_FLOAT,0.145396_FLOAT,0.164308_FLOAT,0.182990_FLOAT,0.205668_FLOAT,&
0.200956_FLOAT,0.299661_FLOAT, &
0.325995_FLOAT,0.305488_FLOAT,0.291723_FLOAT,0.293801_FLOAT,0.29179_FLOAT,0.296729_FLOAT,0.304603_FLOAT,&
0.242041_FLOAT,0.354186_FLOAT,0.350715_FLOAT,&
0.350021_FLOAT,0.345779_FLOAT,0.349532_FLOAT,0.367305_FLOAT,0.382008_FLOAT,0.399709_FLOAT/
data HFtz /&  ! H-Kr
0.007577_FLOAT,0.003312_FLOAT,&
0.086763_FLOAT,0.009962_FLOAT,0.013964_FLOAT,0.005997_FLOAT,0.004731_FLOAT,0.005687_FLOAT,0.006367_FLOAT,0.007511_FLOAT,&
0.077721_FLOAT,0.050003_FLOAT,0.068317_FLOAT,0.041830_FLOAT,0.025796_FLOAT,0.025512_FLOAT,0.023345_FLOAT,0.022734_FLOAT,&
0.097241_FLOAT,0.099167_FLOAT,&
0.219194_FLOAT,0.189098_FLOAT,0.164378_FLOAT,0.147238_FLOAT,0.137298_FLOAT,0.12751_FLOAT,0.118589_FLOAT,0.0318653_FLOAT,&
0.120985_FLOAT,0.0568313_FLOAT, &
0.090996_FLOAT,0.071820_FLOAT,0.063562_FLOAT,0.064241_FLOAT,0.061848_FLOAT,0.061021_FLOAT/
data HFdef2mtzvp /&  ! m def2-TZVP, no f for B-Ne
0.007930_FLOAT,0.003310_FLOAT,&
0.086760_FLOAT,0.009960_FLOAT,0.013960_FLOAT,0.006000_FLOAT,0.003760_FLOAT,0.004430_FLOAT,0.005380_FLOAT,0.006750_FLOAT,&
0.077720_FLOAT,0.050000_FLOAT,0.068320_FLOAT,0.041830_FLOAT,0.025800_FLOAT,0.025510_FLOAT,0.023340_FLOAT,0.022730_FLOAT,&
0.097240_FLOAT,0.099170_FLOAT,&
0.219190_FLOAT,0.189100_FLOAT,0.164380_FLOAT,0.147240_FLOAT,0.137300_FLOAT,0.127510_FLOAT,0.118590_FLOAT,0.031870_FLOAT,&
0.120990_FLOAT,0.056830_FLOAT,&
0.091000_FLOAT,0.071820_FLOAT,0.063560_FLOAT,0.064240_FLOAT,0.061850_FLOAT,0.061020_FLOAT/
data HFvmb/&
0.042400_FLOAT,0.028324_FLOAT,&
0.252661_FLOAT,0.197201_FLOAT,0.156009_FLOAT,0.164586_FLOAT,0.169273_FLOAT,0.214704_FLOAT,0.729138_FLOAT,0.336072_FLOAT,&
0.262329_FLOAT,0.231722_FLOAT,0.251169_FLOAT,0.287795_FLOAT,0.213590_FLOAT,0.250524_FLOAT,0.728579_FLOAT,0.260658_FLOAT, &
2*0._FLOAT,&
16*0._FLOAT/
data HFminisd/& !Al-Ar MINIS + Ahlrichs "P" funktions
0.0_FLOAT,0.0_FLOAT,&
8*0.0_FLOAT,&
2*0.0_FLOAT,1.446950_FLOAT,1.610980_FLOAT,1.766610_FLOAT,1.988230_FLOAT,2.228450_FLOAT,2.487960_FLOAT,&
2*0.0_FLOAT,&
16*0._FLOAT/
data HFlanl2/ & !  LANL2TZ+ vs LANL2DZ (ORCA)_FLOAT, only Sc-Zn
0.102545_FLOAT,0.0719529_FLOAT,0.0491798_FLOAT,0.0362976_FLOAT,0.0266369_FLOAT,0.0235484_FLOAT,0.0171578_FLOAT,0.0438906_FLOAT,&
0.0100259_FLOAT,0.016208_FLOAT/
data HFpobtz /  & ! H-Kr_FLOAT, no RG
0.010077_FLOAT,0.000000_FLOAT,&
0.173239_FLOAT,0.101973_FLOAT,0.131181_FLOAT,0.032234_FLOAT,0.011630_FLOAT,0.008475_FLOAT,0.011673_FLOAT,0.000000_FLOAT,&
0.240653_FLOAT,0.661819_FLOAT,0.522306_FLOAT,0.14163_FLOAT,0.052456_FLOAT,0.184547_FLOAT,0.040837_FLOAT,0.000000_FLOAT,&
0.318136_FLOAT,0.564721_FLOAT,&
0.523194_FLOAT,0.767449_FLOAT,0.620122_FLOAT,0.390227_FLOAT,0.237571_FLOAT,0.247940_FLOAT,0.249589_FLOAT,0.117864_FLOAT,&
&0.325725_FLOAT,0.197183_FLOAT,&
0.264489_FLOAT,0.180375_FLOAT,0.111262_FLOAT,0.147239_FLOAT,0.081747_FLOAT,0.000000_FLOAT/
data HFpobdzvp /  & ! H-Kr_FLOAT, no RG
0.008043_FLOAT,0.000000_FLOAT,&
0.172918_FLOAT,0.114769_FLOAT,0.485944_FLOAT,0.216638_FLOAT,0.088129_FLOAT,0.099072_FLOAT,0.135393_FLOAT,0.000000_FLOAT,&
0.284580_FLOAT,0.210684_FLOAT,0.293510_FLOAT,0.239443_FLOAT,0.212742_FLOAT,0.157717_FLOAT,0.205039_FLOAT,0.000000_FLOAT,&
0.448751_FLOAT,0.495005_FLOAT,&
0.52108_FLOAT,0.374506_FLOAT,0.674866_FLOAT,0.605531_FLOAT,0.37337_FLOAT,0.68191_FLOAT,0.603594_FLOAT,0.193774_FLOAT,&
0.565269_FLOAT,0.566908_FLOAT,0.484427_FLOAT,0.514235_FLOAT,0.390789_FLOAT,0.497878_FLOAT,0.440349_FLOAT,0.000000_FLOAT/
data HF2g / &
0.065950_FLOAT,0.116136_FLOAT,&
0.138319_FLOAT,0.151109_FLOAT,0.219802_FLOAT,0.391898_FLOAT,0.595040_FLOAT,0.993557_FLOAT,1.529560_FLOAT,2.210120_FLOAT,&
0.132269_FLOAT,0.157467_FLOAT,0.152204_FLOAT,0.103184_FLOAT,0.129863_FLOAT,0.184039_FLOAT,0.247297_FLOAT,0.449657_FLOAT,&
0.000000_FLOAT,0.000000_FLOAT,&
0.000000_FLOAT,0.000000_FLOAT,0.000000_FLOAT,0.000000_FLOAT,0.000000_FLOAT,0.000000_FLOAT,0.000000_FLOAT,0.000000_FLOAT,&
0.000000_FLOAT,0.000000_FLOAT,&
0.000000_FLOAT,0.000000_FLOAT,0.000000_FLOAT,0.000000_FLOAT,0.000000_FLOAT,0.000000_FLOAT/
data HFdef1tzvp /&
0.007577_FLOAT,0.003312_FLOAT,&
0.136371_FLOAT,0.011163_FLOAT,0.017129_FLOAT,0.008140_FLOAT,0.005826_FLOAT,0.006770_FLOAT,0.007108_FLOAT,0.008132_FLOAT,&
0.134992_FLOAT,0.147417_FLOAT,0.085253_FLOAT,0.054238_FLOAT,0.033790_FLOAT,0.032862_FLOAT,0.029038_FLOAT,0.026555_FLOAT,&
0.141595_FLOAT,0.207980_FLOAT,& 
0.223252_FLOAT,0.193038_FLOAT,0.167892_FLOAT,0.148726_FLOAT,0.140473_FLOAT,0.130220_FLOAT,0.121166_FLOAT,0.113839_FLOAT,&
0.121855_FLOAT,0.107138_FLOAT,&
0.105637_FLOAT,0.086639_FLOAT,0.075084_FLOAT,0.075089_FLOAT,0.070868_FLOAT,0.068706_FLOAT/
data HFccdz / &
0.007907_FLOAT,0.008287_FLOAT,&
0.047380_FLOAT,0.014240_FLOAT,0.022133_FLOAT,0.014999_FLOAT,0.018148_FLOAT,0.028240_FLOAT,0.042261_FLOAT,0.061485_FLOAT,&
0.073185_FLOAT,0.056218_FLOAT,0.082660_FLOAT,0.052975_FLOAT,0.033874_FLOAT,0.034056_FLOAT,0.031433_FLOAT,0.030034_FLOAT,&
0.000000_FLOAT,0.078016_FLOAT,& !no k cc-pVDZ Basis
0.036885_FLOAT,0.038540_FLOAT,0.036474_FLOAT,0.036061_FLOAT,0.030289_FLOAT,0.027959_FLOAT,0.025177_FLOAT,0.022709_FLOAT,&
0.027386_FLOAT,0.015816_FLOAT,&
0.135176_FLOAT,0.115515_FLOAT,0.102761_FLOAT,0.102967_FLOAT,0.097891_FLOAT,0.097331_FLOAT/
data HFaccdz / & !for li_FLOAT,be_FLOAT,na_FLOAT,mg_FLOAT,k-zn energy below def2-QZVPD reference
0.001183_FLOAT,0.005948_FLOAT,&
0.000000_FLOAT,0.000000_FLOAT,0.005269_FLOAT,0.006380_FLOAT,0.011700_FLOAT,0.021199_FLOAT,0.034160_FLOAT,0.051481_FLOAT,&
0.000000_FLOAT,0.000000_FLOAT,0.016018_FLOAT,0.009268_FLOAT,0.010076_FLOAT,0.015153_FLOAT,0.016889_FLOAT,0.018563_FLOAT,&
0.000000_FLOAT,0.000000_FLOAT,&
0.000000_FLOAT,0.000000_FLOAT,0.000000_FLOAT,0.000000_FLOAT,0.000000_FLOAT,0.000000_FLOAT,0.000000_FLOAT,0.000000_FLOAT,&
0.000000_FLOAT,0.000000_FLOAT,&
0.069963_FLOAT,0.065687_FLOAT,0.072944_FLOAT,0.077585_FLOAT,0.078777_FLOAT,0.080746_FLOAT/
data HFdzp / &
0.008107_FLOAT,0.008045_FLOAT,&
0.136751_FLOAT,0.016929_FLOAT,0.026729_FLOAT,0.021682_FLOAT,0.027391_FLOAT,0.040841_FLOAT,0.058747_FLOAT,0.082680_FLOAT,&
0.153286_FLOAT,0.162296_FLOAT,0.102704_FLOAT,0.073144_FLOAT,0.056217_FLOAT,0.061333_FLOAT,0.065045_FLOAT,0.071398_FLOAT,&
0.145642_FLOAT,0.212865_FLOAT,&
0.232821_FLOAT,0.204796_FLOAT,0.182933_FLOAT,0.169554_FLOAT,0.164701_FLOAT,0.160112_FLOAT,0.157723_FLOAT,0.158037_FLOAT,&
0.179104_FLOAT,0.169782_FLOAT,0.159396_FLOAT,0.140611_FLOAT,0.129645_FLOAT,0.132664_FLOAT,0.132121_FLOAT,0.134081_FLOAT/
data HFmoddz / &
0.008107_FLOAT,0.000000_FLOAT,&
0.113583_FLOAT,0.028371_FLOAT,0.026921_FLOAT,0.021817_FLOAT,0.027249_FLOAT,0.0405982_FLOAT,0.058547_FLOAT,0.0000000_FLOAT,&
0.181140_FLOAT,0.125558_FLOAT,0.167188_FLOAT,0.149843_FLOAT,0.145396_FLOAT,0.1643080_FLOAT,0.182990_FLOAT,0.0000000_FLOAT,&
0.200956_FLOAT,0.299661_FLOAT, &
0.325995_FLOAT,0.305488_FLOAT,0.291723_FLOAT,0.293801_FLOAT,0.291790_FLOAT,0.2967290_FLOAT,0.304603_FLOAT,0.242041_FLOAT,&
0.354186_FLOAT,0.350715_FLOAT, 0.350021_FLOAT,0.345779_FLOAT,0.349532_FLOAT,0.367305_FLOAT,0.382008_FLOAT,0.0000000_FLOAT/
data HFmsvp / &     !H-Kr modified Ahlrichs DZ, supplemented by def2-SV(P)
0.000000_FLOAT,0.000000_FLOAT,& !RG,H set to zero,  F adjusted empirically, Be corrected due to ROHF problems
0.107750_FLOAT,0.020000_FLOAT,0.026850_FLOAT,0.021740_FLOAT,0.027250_FLOAT,0.039930_FLOAT,0.030000_FLOAT,0.000000_FLOAT,&
0.153290_FLOAT,0.162300_FLOAT,0.102700_FLOAT,0.073140_FLOAT,0.056220_FLOAT,0.061330_FLOAT,0.065040_FLOAT,0.000000_FLOAT,&
0.200960_FLOAT,0.299660_FLOAT,&
0.325990_FLOAT,0.305490_FLOAT,0.291720_FLOAT,0.293800_FLOAT,0.291790_FLOAT,0.296730_FLOAT,0.304600_FLOAT,0.242040_FLOAT,&
0.354190_FLOAT,0.350720_FLOAT,&
0.350020_FLOAT,0.345780_FLOAT,0.349530_FLOAT,0.367310_FLOAT,0.382010_FLOAT,0.000000_FLOAT/
!LD
data HFmsvp_ld / &     !H-Kr modified Ahlrichs DZ, supplemented by def2-SV(P)
0.000000_FLOAT,0.000000_FLOAT,& !RG,H set to zero,  F adjusted empirically, Be corrected due to ROHF problems
0.194704_FLOAT,0.098185_FLOAT,0.042941_FLOAT,0.021740_FLOAT,0.027250_FLOAT,0.039930_FLOAT,0.030000_FLOAT,0.000000_FLOAT,&
0.347280_FLOAT,0.311334_FLOAT,0.192999_FLOAT,0.166911_FLOAT,0.159052_FLOAT,0.177211_FLOAT,0.195154_FLOAT,0.000000_FLOAT,&
0.286553_FLOAT,0.391766_FLOAT,&
0.328375_FLOAT,0.280359_FLOAT,0.244807_FLOAT,0.213519_FLOAT,0.210952_FLOAT,0.200381_FLOAT,0.193756_FLOAT,0.191063_FLOAT,&
0.199032_FLOAT,0.199518_FLOAT,&
0.176893_FLOAT,0.152992_FLOAT,0.139637_FLOAT,0.143083_FLOAT,0.142691_FLOAT,0.000000_FLOAT/
data HFrevminix/ &
0.042400_FLOAT,0.028324_FLOAT,&
0.261883_FLOAT,0.193520_FLOAT,0.224237_FLOAT,0.279950_FLOAT,0.357906_FLOAT,0.479012_FLOAT,0.638518_FLOAT,0.832349_FLOAT, &
1.143067_FLOAT,1.252507_FLOAT,1.353291_FLOAT,1.610980_FLOAT,1.766610_FLOAT,1.988230_FLOAT,2.228450_FLOAT,2.487960_FLOAT, &
0.486819_FLOAT,0.621753_FLOAT, &
0.526123_FLOAT,0.476202_FLOAT,0.442710_FLOAT,0.293943_FLOAT,0.555710_FLOAT,0.412730_FLOAT,0.374195_FLOAT,&
0.418722_FLOAT,0.449744_FLOAT,0.444591_FLOAT,&
0.388443_FLOAT,0.376226_FLOAT,0.376473_FLOAT,0.396153_FLOAT,0.412569_FLOAT,0.399709_FLOAT/
data HFpobtzrev2 / &
0.008378_FLOAT,0.000000_FLOAT,&
0.185316_FLOAT,0.112305_FLOAT,0.071828_FLOAT,0.016747_FLOAT,0.011631_FLOAT,0.009703_FLOAT,0.014874_FLOAT,0.000000_FLOAT,&
0.255998_FLOAT,0.262530_FLOAT,0.207733_FLOAT,0.090531_FLOAT,0.053056_FLOAT,0.112989_FLOAT,0.042341_FLOAT,0.000000_FLOAT,&
0.341508_FLOAT,0.481942_FLOAT,&
0.529221_FLOAT,0.427008_FLOAT,0.370708_FLOAT,0.298788_FLOAT,0.271848_FLOAT,0.264929_FLOAT,0.291093_FLOAT,0.295882_FLOAT,&
0.291338_FLOAT,0.220736_FLOAT,&
0.212699_FLOAT,0.131929_FLOAT,0.113791_FLOAT,0.093396_FLOAT,0.095653_FLOAT,0.000000_FLOAT/
data HF631GDP/ &! H-Ca + Br (no 3d)
0.009619_FLOAT,0.008147_FLOAT,&
0.069260_FLOAT,0.030540_FLOAT,0.032736_FLOAT,0.021407_FLOAT,0.024248_FLOAT,0.036746_FLOAT,0.052733_FLOAT,0.075120_FLOAT,&
0.104255_FLOAT,0.070691_FLOAT,0.100260_FLOAT,0.072534_FLOAT,0.054099_FLOAT,0.056408_FLOAT,0.056025_FLOAT,0.057578_FLOAT,&
0.079198_FLOAT,0.161462_FLOAT,&
10*0.0_FLOAT, &
0.000000_FLOAT,0.000000_FLOAT,0.000000_FLOAT,0.000000_FLOAT,0.381049_FLOAT,0.000000_FLOAT/
data BASsv/2*2,2*3,6*9,2*7,6*13,2*11,10*21,6*27/           
data BASminis/2*1,2*2,6*5,2*6,6*9,2*10,16*0/
data BAS631gd/2,5,8*14,8*18,2*22,10*0,4*0,32,0/
data BASsvp/2*5,9,9,6*14,15,18,6*18,24,24,10*31,6*32/  
data oldBASsvp/2*5,6,9,6*14,2*10,6*18,14,24,10*31,6*32/
data BAStz/2*6,14,19,6*31,2*32,6*37,33,36,9*45,48,6*48/  
data BASdef2mtzvp/2*6,14,19,6*24,2*32,6*37,33,36,9*45,48,6*48/   !def2-TZVP, no f for B-Ne
data BASvmb/2*1,2*2,6*4,2*1,6*4,2*0,16*0/ ! minimal basis set with ECPs 
data BASminisd/2*0,2*0,6*0,2*0,6*14,2*0,16*0/
data BASlanl2/22, 22, 22, 22, 22, 22, 22, 22, 22, 18/ ! Sc-Zn LANL2DZ
data BASpobtz / & ! H-Kr no RG
6,0,&
7,7,18,18,18,18,18,0,&
19,19,22,22,22,22,22,0,&
23,23,&
40,40,40,40,40,40,40,40,40,40,&
41,41,41,41,41,0 /
data BASpobdzvp / & ! H-KR, no RG
5,0,&
6,6,14,14,14,14,14,0,&
15,18,18,18,18,18,18,0,&
19,19,&
31,31,31,31,31,31,31,31,31,31,&
32,32,32,32,32,0/
data BAS2g / &
1,1,&
5,5,5,5,5,5,5,5,&
9,9,9,9,9,9,9,9,&
0,0,&
0,0,0,0,0,0,0,0,0,0,&
0,0,0,0,0,0/
data BASdef1tzvp / &
6,6,&
8,11,19,19,19,19,19,19,&
14,14,22,22,22,22,22,22,&
33,33,33,33,33,33,33,33,33,33,&
18,28,&
36,36,36,36,36,36/
data BASccdz / &
5,5,&
14,14,14,14,14,14,14,14,&
18,18,18,18,18,18,18,18,&
0,27,&
43,43,43,43,43,43,43,43,43,43,&
27,27,27,27,27,27/
data BASaccdz / &
9,9,&
23,23,23,23,23,23,23,23,&
27,27,27,27,27,27,27,27,&
0,0,&
59,59,59,59,59,59,59,59,59,59,&
36,36,36,36,36,36/
data BASdzp / &
5,5,&
7,10,15,15,15,15,15,15,&
15,15,23,23,23,23,23,23,&
26,36,&
41,41,41,41,41,41,41,41,41,41,&
41,41,41,41,41,41/
data BASmoddz / &
2,2,&
9,9,10,10,15,15,15,15,&
15,18,6*18,24,24,10*31,6*32/
data BASmsvp / &  ! modified Ahlrichs DZ, supplemented by def2-SV(P)
2,2,&
10,10,15,15,15,15,15,15,&
15,18,18,18,18,18,18,18,&
24,24,&
31,31,31,31,31,31,31,31,31,31,&
32,32,32,32,32,32/
!LD
data BASrevminix / &
1,1,&
7,7,5,5,5,5,5,5,&
11,11,19,14,14,14,14,14,&
11,11,&
21,21,21,21,21,21,21,21,21,21,&
32,32,32,32,32,32/
data BASpobtzrev2 / &
6,0,&
7,7,18,18,18,18,18,0,&
19,19,22,22,22,22,22,0,&
23,23,&
40,40,40,40,40,40,40,40,40,40,&
41,41,41,41,41,0/
data BAS631gdp/5,5,8*14,8*18,2*22,10*0,4*0,32,0/
emiss(1:MAX_ELEM)=0._FLOAT
nbas(1:MAX_ELEM)=0
call lower_case(method)
do while (index(method,'-').ne.0)
  i=scan(method,'-')
  method=trim(method(:(i-1)))//trim(method((i+1):))
  method=trim(method)
enddo
select case (method)
  case ('hf/sv') ! RMS=0.3218975
     emiss(1:MAX_PARA)=HFsv(1:MAX_PARA)
     nbas(1:MAX_PARA)=BASsv(1:MAX_PARA)
     p(1)=0.1724_FLOAT
     p(2)=1.2804_FLOAT
     p(3)=0.8568_FLOAT
     p(4)=1.2342_FLOAT
  case ('hf/svp','hf/def2svp') ! RMS=0.4065  ! fit does not include Li,Na,Mg,K, so not change here
     emiss(1:MAX_PARA)=HFsvp(1:MAX_PARA)
     nbas(1:MAX_PARA)=BASsvp(1:MAX_PARA)
     p(1)=0.2054_FLOAT
     p(2)=1.3157_FLOAT
     p(3)=0.8136_FLOAT
     p(4)=1.2572_FLOAT
  case('hf/sv(p)','hf/def2sv(p)') !RMS=0.3502
     emiss(1:MAX_PARA)=HFsvp(1:MAX_PARA)
     nbas(1:MAX_PARA)=BASsvp(1:MAX_PARA)
     emiss(1)=HFsv(1)
     nbas(1)=BASsv(1)
     p(1)=0.1373_FLOAT
     p(2)=1.4271_FLOAT
     p(3)=0.8141_FLOAT
     p(4)=1.2760_FLOAT
  case ('hf/svp_old') ! RMS=0.4065
     emiss(1:MAX_PARA)=oldHFsvp(1:MAX_PARA)
     nbas(1:MAX_PARA)=oldBASsvp(1:MAX_PARA)
     p(1)=0.2054_FLOAT
     p(2)=1.3157_FLOAT
     p(3)=0.8136_FLOAT
     p(4)=1.2572_FLOAT
  case ('hf/631gd','hf/631gs') ! RMS= 0.40476
     emiss(1:MAX_PARA)=HF631gd(1:MAX_PARA)
     nbas(1:MAX_PARA)=BAS631gd(1:MAX_PARA)
     p(1)=0.2048_FLOAT
     p(2)=1.5652_FLOAT
     p(3)=0.9447_FLOAT
     p(4)=1.2100_FLOAT
  case ('hf/631gdp') ! RMS= 0.40476
     emiss(1:MAX_PARA)=HF631GDP(1:MAX_PARA)
     nbas(1:MAX_PARA)=BAS631gdp(1:MAX_PARA)
     p(1)=0.2048_FLOAT
     p(2)=1.5652_FLOAT
     p(3)=0.9447_FLOAT
     p(4)=1.2100_FLOAT
  case ('hf/minis') ! RMS= 0.3040
    emiss(1:MAX_PARA)=HFminis(1:MAX_PARA)
    nbas(1:MAX_PARA)=BASminis(1:MAX_PARA)
     p(1)=0.1290_FLOAT
     p(2)=1.1526_FLOAT
     p(3)=1.1549_FLOAT
     p(4)=1.1763_FLOAT
  case ('hf/minix') 
    emiss(1:12)=HFminis(1:12)
    nbas(1:12)=BASminis(1:12)
    emiss(13:18)=HFminisd(13:18)
    nbas(13:18)=BASminisd(13:18)
    emiss(19:30)=HFsv(19:30)
    nbas(19:30)=BASsv(19:30)
    emiss(31:36)=HFsvp(31:36)
    nbas(31:36)=BASsvp(31:36)
    emiss(3)=0.177871_FLOAT
    emiss(4)=0.171596_FLOAT
    nbas(3:4)=5
    emiss(11)=1.114110_FLOAT
    emiss(12)=1.271150_FLOAT
    nbas(11:12)=9
     p(1)=0.1290_FLOAT
     p(2)=1.1526_FLOAT
     p(3)=1.1549_FLOAT
     p(4)=1.1763_FLOAT
   case ('hfsol/minix')
    emiss(1:36)=HFrevminix(1:36)
    nbas(1:36)=BASrevminix(1:36)
     p(1)=0.1290_FLOAT
     p(2)=1.1526_FLOAT
     p(3)=1.1549_FLOAT
     p(4)=1.1763_FLOAT
  case ('hf/tz','hf/def2tzvp') !  RMS= 0.1150
     emiss(1:MAX_PARA)=HFtz(1:MAX_PARA)
     nbas(1:MAX_PARA)=BAStz(1:MAX_PARA)
     p(1)=0.3127_FLOAT
     p(2)=1.9914_FLOAT
     p(3)=1.0216_FLOAT
     p(4)=1.2833_FLOAT
  case ('hf/deftzvp', 'hf/tzvp') ! RMS=0.209
     emiss(1:MAX_PARA)=HFdef1tzvp(1:MAX_PARA)
     nbas(1:MAX_PARA)=BASdef1tzvp(1:MAX_PARA)
     p(1)=0.2600_FLOAT
     p(2)=2.2448_FLOAT
     p(3)=0.7998_FLOAT
     p(4)=1.4381_FLOAT
  case ('hf/ccdz', 'hf/ccpvdz') ! RMS=0.4968
     emiss(1:MAX_PARA)=HFccdz(1:MAX_PARA)
     nbas(1:MAX_PARA)=BASccdz(1:MAX_PARA)
     p(1)=0.4416_FLOAT
     p(2)=1.5185_FLOAT
     p(3)=0.6902_FLOAT
     p(4)=1.3713_FLOAT
  case ('hf/accdz','hf/augccpvdz') !RMS=0.2222
     emiss(1:MAX_PARA)=HFaccdz(1:MAX_PARA)
     nbas(1:MAX_PARA)=BASaccdz(1:MAX_PARA)
     p(1)=0.0748_FLOAT
     p(2)=0.0663_FLOAT
     p(3)=0.3811_FLOAT
     p(4)=1.0155_FLOAT
  case ('hf/2g')
     emiss(1:18)=HF2g(1:18)
     nbas(1:18)=BAS2g(1:18)
     emiss(19:30)=HFsv(19:30)
     nbas(19:30)=BASsv(19:30)
     emiss(31:36)=HFsvp(31:36)
     nbas(31:36)=BASsvp(31:36)
     p(1)=0.1061_FLOAT
     p(2)=1.4239_FLOAT
     p(3)=0.8699_FLOAT
     p(4)=1.5000_FLOAT
  case('hf/dzp') !RMS=0.4571
     emiss(1:MAX_PARA)=HFdzp(1:MAX_PARA)
     nbas(1:MAX_PARA)=BASdzp(1:MAX_PARA)
     p(1)=0.1443_FLOAT
     p(2)=1.4547_FLOAT
     p(3)=0.3711_FLOAT
     p(4)=1.6300_FLOAT
  case ('dft/lanl','b3lyp/lanl')
    emiss(1:MAX_PARA)=HF631gd(1:MAX_PARA)
    nbas(1:MAX_PARA)=BAS631gd(1:MAX_PARA)
     p(1)=0.3405_FLOAT
     p(2)=1.6127_FLOAT
     p(3)=0.8589_FLOAT
     p(4)=1.2830_FLOAT
     emiss(21:30)=HFlanl2(1:10)
     nbas(21:30)=BASlanl2(1:10)
  case ('dft/sv','b3lyp/sv') ! RMS= 0.557
     emiss(1:MAX_PARA)=HFsv(1:MAX_PARA)
     nbas(1:MAX_PARA)=BASsv(1:MAX_PARA)
     p(1)=0.4048_FLOAT
     p(2)=1.1626_FLOAT
     p(3)=0.8652_FLOAT
     p(4)=1.2375_FLOAT
  case ('dft/sv(p)','b3lyp/sv(p)','b3lyp/def2sv(p)','pbe/def2sv(p)','dft/def2sv(p)','dft/sv_p') ! RMS= 0.57 ! def2-SV(P)  
     emiss(1:MAX_PARA)=HFsvp(1:MAX_PARA)
     nbas(1:MAX_PARA)=BASsvp(1:MAX_PARA)
     emiss(1)=HFsv(1)
     nbas(1)=BASsv(1)
     p(1)=0.2424_FLOAT
     p(2)=1.2371_FLOAT
     p(3)=0.6076_FLOAT
     p(4)=1.4078_FLOAT
  case ('dft/svx','b3lyp/svx') ! RMS=  0.56 ! def2-SV(P/h,c)  = SV at h,c
     emiss(1:MAX_PARA)=HFsvp(1:MAX_PARA)
     nbas(1:MAX_PARA)=BASsvp(1:MAX_PARA)
     emiss(1)=HFsv(1)
     nbas(1)=BASsv(1)
     emiss(6)=HFsv(6)
     nbas(6)=BASsv(6)
     p(1)=0.1861_FLOAT
     p(2)=1.3200_FLOAT
     p(3)=0.6171_FLOAT
     p(4)=1.4019_FLOAT
  case ('dft/svp','b3lyp/svp') ! RMS=0.6498
     emiss(1:MAX_PARA)=HFsvp(1:MAX_PARA)
     nbas(1:MAX_PARA)=BASsvp(1:MAX_PARA)
     p(1)=0.2990_FLOAT
     p(2)=1.2605_FLOAT
     p(3)=0.6438_FLOAT
     p(4)=1.3694_FLOAT
  case ('dft/svp_old','b3lyp/svp_old') ! RMS=0.6498
     emiss(1:MAX_PARA)=oldHFsvp(1:MAX_PARA)
     nbas(1:MAX_PARA)=oldBASsvp(1:MAX_PARA)
     p(1)=0.2990_FLOAT
     p(2)=1.2605_FLOAT
     p(3)=0.6438_FLOAT
     p(4)=1.3694_FLOAT
  case ('dft/631gd','b3lyp/631gd','dft/631gs','b3lyp/631gs') ! RMS=  0.47856
    emiss(1:MAX_PARA)=HF631gd(1:MAX_PARA)
    nbas(1:MAX_PARA)=BAS631gd(1:MAX_PARA)
     p(1)=0.3405_FLOAT
     p(2)=1.6127_FLOAT
     p(3)=0.8589_FLOAT
     p(4)=1.2830_FLOAT
  case ('dft/minis','b3lyp/minis') ! RMS= 0.3400
    emiss(1:MAX_PARA)=HFminis(1:MAX_PARA)
    nbas(1:MAX_PARA)=BASminis(1:MAX_PARA)
     p(1)=0.2059_FLOAT
     p(2)=0.9722_FLOAT
     p(3)=1.1961_FLOAT
     p(4)=1.1456_FLOAT
  case('dft/minix','b3lyp/minix')
    emiss(1:12)=HFminis(1:12)
    nbas(1:12)=BASminis(1:12)
    emiss(13:18)=HFminisd(13:18)
    nbas(13:18)=BASminisd(13:18)
    emiss(19:30)=HFsv(19:30)
    nbas(19:30)=BASsv(19:30)
    emiss(31:36)=HFsvp(31:36)
    nbas(31:36)=BASsvp(31:36)
    emiss(3)=0.177871_FLOAT
    emiss(4)=0.171596_FLOAT
    nbas(3:4)=5
    emiss(11)=1.114110_FLOAT
    emiss(12)=1.271150_FLOAT
    nbas(11:12)=9
    p(1)=0.2059_FLOAT
    p(2)=0.9722_FLOAT
    p(3)=1.1961_FLOAT
    p(4)=1.1456_FLOAT
  case ('dft/tz','b3lyp/tz','dft/def2tzvp','b3lyp/def2tzvp') ! RMS=0.19648
     emiss(1:MAX_PARA)=HFtz(1:MAX_PARA)
     nbas(1:MAX_PARA)=BAStz(1:MAX_PARA) 
     p(1)=0.2905_FLOAT
     p(2)=2.2495_FLOAT
     p(3)=0.8120_FLOAT
     p(4)=1.4412_FLOAT
  case ('dft/deftzvp','b3lyp/deftzvp', 'dft/tzvp', 'b3lyp/tzvp') ! RMS=0.1817
     emiss(1:MAX_PARA)=HFdef1tzvp(1:MAX_PARA)
     nbas(1:MAX_PARA)=BASdef1tzvp(1:MAX_PARA)
     p(1)=0.2393_FLOAT
     p(2)=2.2247_FLOAT
     p(3)=0.8185_FLOAT
     p(4)=1.4298_FLOAT
  case ('dft/ccdz','dft/ccpvdz','b3lyp/ccdz','b3lyp/ccpvdz') !RMS=0.7610
     emiss(1:MAX_PARA)=HFccdz(1:MAX_PARA)
     nbas(1:MAX_PARA)=BASccdz(1:MAX_PARA)
     p(1)=0.5383_FLOAT
     p(2)=1.6482_FLOAT
     p(3)=0.6230_FLOAT
     p(4)=1.4523_FLOAT
  case ('dft/accdz','dft/augccpvdz','b3lyp/accdz','b3lyp/augccpvdz') !RMS=0.1840
     emiss(1:MAX_PARA)=HFaccdz(1:MAX_PARA)
     nbas(1:MAX_PARA)=BASaccdz(1:MAX_PARA)
     p(1)=0.1465_FLOAT
     p(2)=0.0500_FLOAT
     p(3)=0.6003_FLOAT
     p(4)=0.8761_FLOAT
  case ('dft/pobtz','b3lyp/pobtz') ! RMS=
     emiss(1:MAX_PARA)=HFpobtz(1:MAX_PARA)
     nbas(1:MAX_PARA)=BASpobtz(1:MAX_PARA)
     p(1)=0.1300_FLOAT
     p(2)=1.3743_FLOAT
     p(3)=0.4792_FLOAT
     p(4)=1.3962_FLOAT
  case ('dft/pobdzvp','b3lyp/pobdzvp') ! RMS=0.5851
     emiss(1:MAX_PARA)=HFpobdzvp(1:MAX_PARA)
     nbas(1:MAX_PARA)=BASpobdzvp(1:MAX_PARA)
     p(1)=0.1124_FLOAT
     p(2)=1.8146_FLOAT
     p(3)=0.7180_FLOAT
     p(4)=1.4134_FLOAT
  case('dft/dzp','b3lyp/dzp') !RMS=0.7184
     emiss(1:MAX_PARA)=HFdzp(1:MAX_PARA)
     nbas(1:MAX_PARA)=BASdzp(1:MAX_PARA)
     p(1)=0.2687_FLOAT
     p(2)=1.4634_FLOAT
     p(3)=0.3513_FLOAT
     p(4)=1.6880_FLOAT

  case ('blyp/minis','gga/minis') ! RMS= 0.3462
  emiss(1:MAX_PARA)=HFminis(1:MAX_PARA)
  nbas(1:MAX_PARA)=BASminis(1:MAX_PARA)
  p(1)=0.1566_FLOAT
  p(2)=1.0271_FLOAT
  p(3)=1.0732_FLOAT
  p(4)=1.1968_FLOAT
  case ('tpss/minis') ! RMS= 
  emiss(1:MAX_PARA)=HFminis(1:MAX_PARA)
  nbas(1:MAX_PARA)=BASminis(1:MAX_PARA)
  p(1)=0.22982_FLOAT
  p(2)=1.35401_FLOAT
  p(3)=1.47633_FLOAT
  p(4)=1.11300_FLOAT
  case ('tpss/svp','tpss/def2svp') ! RMS=  0.61
  emiss(1:MAX_PARA)=HFsvp(1:MAX_PARA)
  nbas(1:MAX_PARA)=BASsvp(1:MAX_PARA)
  p(1)=0.6647_FLOAT
  p(2)=1.3306_FLOAT
  p(3)=1.0792_FLOAT
  p(4)=1.1651_FLOAT
  case ('gga/svp','blyp/svp','gga/def2svp','blyp/def2svp') ! RMS=
  emiss(1:MAX_PARA)=HFsvp(1:MAX_PARA)
  nbas(1:MAX_PARA)=BASsvp(1:MAX_PARA)
  p(1)=0.6823_FLOAT
  p(2)=1.2491_FLOAT
  p(3)=0.8225_FLOAT
  p(4)=1.2811_FLOAT
  case ('blyp/sv','gga/sv') ! RMS = 0.6652
  emiss(1:MAX_PARA)=HFsv(1:MAX_PARA)
  nbas(1:MAX_PARA)=BASsv(1:MAX_PARA)
  p(1)=0.2727_FLOAT
  p(2)=1.4022_FLOAT
  p(3)=0.8055_FLOAT
  p(4)=1.3000_FLOAT
  case ('blyp/tz','gga/tz','blyp/def2tzvp','gga/def2tzvp') !RMS = 0.21408
  emiss(1:MAX_PARA)=HFtz(1:MAX_PARA)
  nbas(1:MAX_PARA)=BAStz(1:MAX_PARA)
  p(1)=0.1182_FLOAT
  p(2)=1.0631_FLOAT
  p(3)=1.0510_FLOAT
  p(4)=1.1287_FLOAT
  case ('pw6b95/minis')  ! RMS = 0.3279929
  emiss(1:MAX_PARA)=HFminis(1:MAX_PARA)
  nbas(1:MAX_PARA)=BASminis(1:MAX_PARA)
  p(1)=0.21054_FLOAT
  p(2)=1.25458_FLOAT
  p(3)=1.35003_FLOAT
  p(4)=1.14061_FLOAT
  case ('pw6b95/svp','pw6b95/def2svp')  ! RMS = 0.58312
  emiss(1:MAX_PARA)=HFsvp(1:MAX_PARA)
  nbas(1:MAX_PARA)=BASsvp(1:MAX_PARA)
  p(1)=0.3098_FLOAT
  p(2)=1.2373_FLOAT
  p(3)=0.6896_FLOAT
  p(4)=1.3347_FLOAT
  case ('pbeh3c', 'pbeh3c/msvp')
  emiss(1:MAX_PARA)=HFmsvp(1:MAX_PARA)
  emiss(19:MAX_PARA)=HFdzp(19:MAX_PARA)
  emiss(36)=0.0_FLOAT
  nbas(1:MAX_PARA)=BASmsvp(1:MAX_PARA)
  p(1)=1.00000_FLOAT
  p(2)=1.32492_FLOAT
  p(3)=0.27649_FLOAT
  p(4)=1.95600_FLOAT
  case ('hse3c', 'hse3c/msvp')
  emiss(1:MAX_PARA)=HFmsvp(1:MAX_PARA)
  emiss(19:MAX_PARA)=HFdzp(19:MAX_PARA)
  emiss(36)=0.0_FLOAT
  nbas(1:MAX_PARA)=BASmsvp(1:MAX_PARA)
  p(1)=1.00000_FLOAT
  p(2)=1.40858_FLOAT
  p(3)=0.29083_FLOAT
  p(4)=1.95260_FLOAT
  case('b3pbe3c')
  emiss(1:MAX_PARA)=HFdef2mtzvp(1:MAX_PARA)
  nbas(1:MAX_PARA)=BASdef2mtzvp(1:MAX_PARA)
  p(1)=1.0000_FLOAT
  p(2)=2.98561_FLOAT
  p(3)=0.3011_FLOAT
  p(4)=2.4405_FLOAT
case ('pbesol03c')
  emiss(1:MAX_PARA)=HFmsvp_ld(1:MAX_PARA)
  emiss(36)=0.0_FLOAT
  nbas(1:MAX_PARA)=BASmsvp(1:MAX_PARA)
  p(1)=1.00000_FLOAT
  p(2)=1.36866_FLOAT
  p(3)=0.27464_FLOAT
  p(4)=1.96477_FLOAT
 case ('hsesol3c')
  emiss(1:MAX_PARA)=HFmsvp_ld(1:MAX_PARA)
  emiss(36)=0.0_FLOAT
  nbas(1:MAX_PARA)=BASmsvp(1:MAX_PARA)
  p(1)=1.00000_FLOAT
  p(2)=1.42809_FLOAT
  p(3)=0.29364_FLOAT
  p(4)=1.94671_FLOAT
 case ('r2scansol3c')
  emiss(1:MAX_PARA)=HFmsvp_ld(1:MAX_PARA)
  emiss(36)=0.0_FLOAT
  nbas(1:MAX_PARA)=BASmsvp(1:MAX_PARA)
  p(1)=1.00000_FLOAT
  p(2)=1.33135_FLOAT
  p(3)=0.29710_FLOAT
  p(4)=1.95202_FLOAT
 case ('r2scan0sol3c')
  emiss(1:MAX_PARA)=HFmsvp_ld(1:MAX_PARA)
  emiss(36)=0.0_FLOAT
  nbas(1:MAX_PARA)=BASmsvp(1:MAX_PARA)
  p(1)=1.00000_FLOAT
  p(2)=1.34892_FLOAT
  p(3)=0.29759_FLOAT
  p(4)=1.95599_FLOAT
  case ('r2scanpob3c') ! RMS=
  emiss(1:MAX_PARA)=HFpobtzrev2(1:MAX_PARA)
  nbas(1:MAX_PARA)=BASpobtzrev2(1:MAX_PARA)
  p(1)=1.00000_FLOAT
  p(2)=3.29259_FLOAT
  p(3)=0.34377_FLOAT
  p(4)=2.36190_FLOAT
  case ('r2scan0pob3c') ! RMS=
  emiss(1:MAX_PARA)=HFpobtzrev2(1:MAX_PARA)
  nbas(1:MAX_PARA)=BASpobtzrev2(1:MAX_PARA)
  p(1)=1.00000_FLOAT
  p(2)=3.31761_FLOAT
  p(3)=0.36330_FLOAT
  p(4)=2.31331_FLOAT
case('sv')
     emiss(1:MAX_PARA)=HFsv(1:MAX_PARA)
     nbas(1:MAX_PARA)=BASsv(1:MAX_PARA) 
case('sv(p)','def2sv(p)')
     emiss(1:MAX_PARA)=HFsvp(1:MAX_PARA)
     nbas(1:MAX_PARA)=BASsvp(1:MAX_PARA)
     emiss(1)=0.009037_FLOAT
     nbas(1)=2
case ('svx') ! RMS=  ! def2-SV(P/h,c)  = SV at h,c
     emiss(1:MAX_PARA)=HFsvp(1:MAX_PARA)
     nbas(1:MAX_PARA)=BASsvp(1:MAX_PARA)
     emiss(1)=HFsv(1)
     nbas(1)=BASsv(1)
     emiss(6)=HFsv(6)
     nbas(6)=BASsv(6)
case('svp','def2svp')
     emiss(1:MAX_PARA)=HFsvp(1:MAX_PARA)
     nbas(1:MAX_PARA)=BASsvp(1:MAX_PARA) 
case('minis')
     emiss(1:MAX_PARA)=HFminis(1:MAX_PARA)
     nbas(1:MAX_PARA)=BASminis(1:MAX_PARA) 
case('631gd')
     emiss(1:MAX_PARA)=HF631gd(1:MAX_PARA)
     nbas(1:MAX_PARA)=BAS631gd(1:MAX_PARA) 
case('tz','def2tzvp')
     emiss(1:MAX_PARA)=HFtz(1:MAX_PARA)
     nbas(1:MAX_PARA)=BAStz(1:MAX_PARA) 
case('pobtz')
   emiss(1:MAX_PARA)=HFpobtz(1:MAX_PARA)
   nbas(1:MAX_PARA)=BASpobtz(1:MAX_PARA)
case('pobdzvp')
   emiss(1:MAX_PARA)=HFpobdzvp(1:MAX_PARA)
   nbas(1:MAX_PARA)=BASpobdzvp(1:MAX_PARA)
case ('minix') 
    emiss(1:12)=HFminis(1:12)
    nbas(1:12)=BASminis(1:12)
    emiss(13:18)=HFminisd(13:18)
    nbas(13:18)=BASminisd(13:18)
    emiss(19:30)=HFsv(19:30)
    nbas(19:30)=BASsv(19:30)
    emiss(31:36)=HFsvp(31:36)
    nbas(31:36)=BASsvp(31:36)
    emiss(3)=0.177871_FLOAT
    emiss(4)=0.171596_FLOAT
    nbas(3:4)=5
    emiss(11)=1.114110_FLOAT
    emiss(12)=1.271150_FLOAT
    nbas(11:12)=9
case ('2g')
    emiss(1:18)=HF2g(1:18)
    nbas(1:18)=BAS2g(1:18)
    emiss(19:30)=HFsv(19:30)
    nbas(19:30)=BASsv(19:30)
    emiss(31:36)=HFsvp(31:36)
    nbas(31:36)=BASsvp(31:36)
case('deftzvp')
     emiss(1:MAX_PARA)=HFdef1tzvp(1:MAX_PARA)
     nbas(1:MAX_PARA)=BASdef1tzvp(1:MAX_PARA)
case('ccdz')
     emiss(1:MAX_PARA)=HFccdz(1:MAX_PARA)
     nbas(1:MAX_PARA)=BASccdz(1:MAX_PARA)
case('accdz')
     emiss(1:MAX_PARA)=HFaccdz(1:MAX_PARA)
     nbas(1:MAX_PARA)=BASaccdz(1:MAX_PARA)
case('dzp')
     emiss(1:MAX_PARA)=HFdzp(1:MAX_PARA)
     nbas(1:MAX_PARA)=BASdzp(1:MAX_PARA)
case('moddz')
     emiss(1:MAX_PARA)=HFmoddz(1:MAX_PARA)
     nbas(1:MAX_PARA)=BASmoddz(1:MAX_PARA)
case default
    call errvrs(1,znamz,'BASIS SET NOT PARAMETRIZED')
    call errvrs(0,znamz,'ERROR IN GCP INPUT')
end select
 RETURN
end

subroutine ssovl(r,iat,jat,iz,xza,xzb,ovl)
 USE NUMBERS
 USE PARINF_MODULE
 IMPLICIT DOUBLE PRECISION(A-H,O-Z)
 PARAMETER(THIRD=1._FLOAT/3._FLOAT,RAD13=0.577350269189625731_FLOAT) ! SQRT(1/3)
 PARAMETER(QUA=1._FLOAT/480._FLOAT,TQUA=QUA*THIRD)
 DIMENSION ishell(72),IZ(*)
data ishell/                 &
          1,1               &
          ,2,2,2,2,2,2,2,2, &
          3,3,3,3,3,3,3,3,  & 
          54*3/
       za=xza
       zb=xzb
       na=iz(iat)
       nb=iz(jat)
ii=ishell(na)*ishell(nb)
R05=R*0.5_FLOAT
ax=(za+zb)*R05
bx=(zb-za)*R05

if(za.eq.zb.OR.abs(za-zb).lt.0.1_FLOAT) then
  select case (ii)
   case (1)
    ovl=0.25_FLOAT*sqrt((za*zb*R*R)**3)*(A2_JGB(ax)*Bint(bx,0)-Bint(bx,2)*A0_JGB(ax))
   case (2)  
    ovl =RAD13
    if(ishell(na).lt.ishell(nb)) then
      anorm=SQRT((ZA**3)*(ZB**5))*(R**4)*0.125_FLOAT
      ovl=ovl*anorm*(A3_JGB(ax)*Bint(bx,0)-Bint(bx,3)*A0_JGB(ax)+A2_JGB(ax)*Bint(bx,1)-Bint(bx,2)*A1_JGB(ax))
     else
      xx=za
      za=zb
      zb=xx
      ax=(za+zb)*R05
      bx=(zb-za)*R05
      anorm=SQRT((ZA**3)*(ZB**5))*(R**4)*0.125_FLOAT
      ovl=ovl*anorm*(A3_JGB(ax)*Bint(bx,0)-Bint(bx,3)*A0_JGB(ax)+A2_JGB(ax)*Bint(bx,1)-Bint(bx,2)*A1_JGB(ax))
    endif
   case (4)
    anorm=SQRT((ZA*ZB)**5)*(R**5)*0.0625_FLOAT
    ovl=anorm* (A4_JGB(ax)*Bint(bx,0)+Bint(bx,4)*A0_JGB(ax)-2._FLOAT*A2_JGB(ax)*Bint(bx,2))*THIRD
   case(3) 
    if(ishell(na).lt.ishell(nb)) then
      anorm=SQRT((ZA**3)*(ZB**7)/7.5_FLOAT)*(R**5)*0.0625_FLOAT
      ovl=anorm*(A4_JGB(ax)*Bint(bx,0)-Bint(bx,4)*A0_JGB(ax)+2._FLOAT*(A3_JGB(ax)*Bint(bx,1)-Bint(bx,3)*A1_JGB(ax)))*RAD13
    else
      xx=za
      za=zb
      zb=xx
      ax=(za+zb)*R05
      bx=(zb-za)*R05
      anorm=SQRT((ZA**3)*(ZB**7)/7.5_FLOAT)*(R**5)*0.0625_FLOAT
      ovl=anorm*(A4_JGB(ax)*Bint(bx,0)-Bint(bx,4)*A0_JGB(ax)+2._FLOAT*(A3_JGB(ax)*Bint(bx,1)-Bint(bx,3)*A1_JGB(ax)))*RAD13
    endif
   case(6) 
    if(ishell(na).lt.ishell(nb)) then
      anorm=SQRT((za**5)*(zb**7)/7.5_FLOAT)*(R**6)*0.03125_FLOAT
      ovl=anorm*(A5_JGB(ax)*Bint(bx,0)+A4_JGB(ax)*Bint(bx,1)-2._FLOAT*(A3_JGB(ax)*Bint(bx,2)+A2_JGB(ax)*Bint(bx,3))+&
           A1_JGB(ax)*Bint(bx,4)+A0_JGB(ax)*Bint(bx,5))*THIRD
    else
      xx=za
      za=zb
      zb=xx
      ax=(za+zb)*R05
      bx=(zb-za)*R05
      anorm=SQRT((za**5)*(zb**7)/7.5_FLOAT)*(R**6)*0.03125_FLOAT
      ovl=anorm*(A5_JGB(ax)*Bint(bx,0)+A4_JGB(ax)*Bint(bx,1)-2._FLOAT*(A3_JGB(ax)*Bint(bx,2)+A2_JGB(ax)*Bint(bx,3))+&
           A1_JGB(ax)*Bint(bx,4)+A0_JGB(ax)*Bint(bx,5))*THIRD
    endif
   case(9)
      anorm=sqrt((ZA*ZB*R*R)**7)*QUA
      ovl=anorm*(A6_JGB(ax)*Bint(bx,0)-3._FLOAT*(A4_JGB(ax)*Bint(bx,2)-A2_JGB(ax)*Bint(bx,4))-A0_JGB(ax)*Bint(bx,6))*THIRD
   end select
else 
   select case (ii)
   case (1)
      anorm=0.25_FLOAT*sqrt((za*zb*R*R)**3)
      ovl=(A2_JGB(ax)*B0_JGB(bx)-B2_JGB(bx)*A0_JGB(ax))*anorm
   case (2)
      ovl =RAD13
    if(ishell(na).lt.ishell(nb)) then
      anorm=SQRT((ZA**3)*(ZB**5))*(R**4)*0.125_FLOAT
      ovl=ovl*anorm*(A3_JGB(ax)*B0_JGB(bx)-B3_JGB(bx)*A0_JGB(ax)+A2_JGB(ax)*B1_JGB(bx)-B2_JGB(bx)*A1_JGB(ax))
     else
      xx=za
      za=zb
      zb=xx 
      ax=(za+zb)*R05
      bx=(zb-za)*R05
      anorm=SQRT((ZA**3)*(ZB**5))*(R**4)*0.125_FLOAT
      ovl=ovl*anorm*(A3_JGB(ax)*B0_JGB(bx)-B3_JGB(bx)*A0_JGB(ax)+A2_JGB(ax)*B1_JGB(bx)-B2_JGB(bx)*A1_JGB(ax))
    endif
   case (4) ! <2s|2s>
      anorm=SQRT((ZA*ZB)**5)*(R**5)*0.0625_FLOAT
      ovl=anorm*(A4_JGB(ax)*B0_JGB(bx)+B4_JGB(bx)*A0_JGB(ax)-2.0_FLOAT*A2_JGB(ax)*B2_JGB(bx))*THIRD
   case(3)  ! <1s|3s> + <3s|1s>
    if(ishell(na).lt.ishell(nb)) then
      anorm=SQRT((ZA**3)*(ZB**7)/7.5_FLOAT)*(R**5)*0.0625_FLOAT
      ovl=anorm*(A4_JGB(ax)*B0_JGB(bx)-B4_JGB(bx)*A0_JGB(ax)+2._FLOAT*(A3_JGB(ax)*B1_JGB(bx)-B3_JGB(bx)*A1_JGB(ax)))*RAD13
    else
      xx=za
      za=zb
      zb=xx
      ax=(za+zb)*R05
      bx=(zb-za)*R05
      anorm=SQRT((ZA**3)*(ZB**7)/7.5_FLOAT)*(R**5)*0.0625_FLOAT
      ovl=anorm*(A4_JGB(ax)*B0_JGB(bx)-B4_JGB(bx)*A0_JGB(ax)+2._FLOAT*(A3_JGB(ax)*B1_JGB(bx)-B3_JGB(bx)*A1_JGB(ax)))*RAD13
    endif
   case(6)  ! <2s|3s> + <3s|2s>
    if(ishell(na).lt.ishell(nb)) then
      anorm=SQRT((za**5)*(zb**7)/7.5_FLOAT)*(R**6)*0.03125_FLOAT
      ovl=anorm*(A5_JGB(ax)*B0_JGB(bx)+A4_JGB(ax)*B1_JGB(bx)-2_FLOAT*(A3_JGB(ax)*B2_JGB(bx)+A2_JGB(ax)*B3_JGB(bx))+&
           A1_JGB(ax)*B4_JGB(bx)+A0_JGB(ax)*B5_JGB(bx))*THIRD
    else
      xx=za
      za=zb
      zb=xx
      ax=(za+zb)*R05
      bx=(zb-za)*R05
      anorm=SQRT((za**5)*(zb**7)/7.5_FLOAT)*(R**6)*0.03125_FLOAT
      ovl=anorm*(A5_JGB(ax)*B0_JGB(bx)+A4_JGB(ax)*B1_JGB(bx)-2_FLOAT*(A3_JGB(ax)*B2_JGB(bx)+A2_JGB(ax)*B3_JGB(bx))+&
           A1_JGB(ax)*B4_JGB(bx)+A0_JGB(ax)*B5_JGB(bx))*THIRD
    endif
    case(9) ! <3s|3>
      anorm=sqrt((ZA*ZB*R*R)**7)*TQUA
!      ovl=anorm*(A6(ax)*B0(bx)-3.d0*(A4(ax)*B2(bx)-A2(ax)*B4(bx))-A0(ax)*Bint(bx,6))
      ovl=anorm*(A6_JGB(ax)*B0_JGB(bx)-3._FLOAT*(A4_JGB(ax)*B2_JGB(bx)-A2_JGB(ax)*B4_JGB(bx))-A0_JGB(ax)*B6_JGB(bx))
   end select
endif
return
end subroutine ssovl
 pure function A0_JGB(x)
 USE NUMBERS
 IMPLICIT DOUBLE PRECISION(A-H,O-Z)
real(FLOAT),intent(in) :: x
A0_JGB=exp(-x)/x
return 
end function
pure function A1_JGB(x)
 USE NUMBERS
 IMPLICIT DOUBLE PRECISION(A-H,O-Z)
real(FLOAT),intent(in) :: x
A1_JGB=((1._FLOAT+x)*exp(-x))/(x**2)
return
end function
pure function A2_JGB(x)
 USE NUMBERS
 IMPLICIT DOUBLE PRECISION(A-H,O-Z)
real(FLOAT),intent(in) :: x
A2_JGB=(2._FLOAT+x+x+x**2)*exp(-x)/x**3
return
end function
pure function A3_JGB(x)
 USE NUMBERS
 IMPLICIT DOUBLE PRECISION(A-H,O-Z)
real(FLOAT),intent(in) :: x
x2=x*x
x3=x2*x
x4=x3*x
xx=(6._FLOAT+6._FLOAT*x+3._FLOAT*x2+x3)
A3_JGB=(xx*exp(-x))/x4
return
end function
pure function A4_JGB(x)
 USE NUMBERS
 IMPLICIT DOUBLE PRECISION(A-H,O-Z)
real(FLOAT),intent(in) :: x
x2=x*x
x3=x2*x
x4=x3*x
x5=x4*x
xx=(24._FLOAT+24._FLOAT*x+12._FLOAT*x2+4._FLOAT*x3+x4)
A4_JGB=(xx*exp(-x))/x5
return
end function
pure function A5_JGB(x)
 USE NUMBERS
 IMPLICIT DOUBLE PRECISION(A-H,O-Z)
real(FLOAT),intent(in) :: x
x2=x*x
x3=x2*x
x4=x3*x
x5=x4*x
x6=x5*x
xx=(x*120._FLOAT+x2*60._FLOAT+x3*20._FLOAT+x4*5._FLOAT+x5+120._FLOAT)
A5_JGB=(xx*exp(-x))/x6
return
end function
real(8) pure function A6_JGB(x)
 USE NUMBERS
 IMPLICIT DOUBLE PRECISION(A-H,O-Z)
real(FLOAT),intent(in) :: x
x2=x*x
x3=x2*x
x4=x3*x
x5=x4*x
x6=x5*x
x7=x6*x
xx=(x*720._FLOAT+x2*360._FLOAT+x3*120._FLOAT+x4*30._FLOAT+x5*6._FLOAT+x6+720._FLOAT)
A6_JGB=(xx*exp(-x))/x7
return
end function
pure function B0_JGB(x)
 USE NUMBERS
 IMPLICIT DOUBLE PRECISION(A-H,O-Z)
real(FLOAT),intent(in) :: x
B0_JGB=(exp(x)-exp(-x))/x
return
end function
pure function B1_JGB(x)
 USE NUMBERS
 IMPLICIT DOUBLE PRECISION(A-H,O-Z)
real(FLOAT),intent(in) :: x
x2=x*x
x3=x2*x
B1_JGB=((1._FLOAT-x)*exp(x)-(1._FLOAT+x)*exp(-x))/x2
return
end function
pure function B2_JGB(x)
 USE NUMBERS
 IMPLICIT DOUBLE PRECISION(A-H,O-Z)
real(FLOAT),intent(in) :: x
x2=x*x
x3=x2*x
B2_JGB=(((2._FLOAT-x-x+x2)*exp(x)) - ((2._FLOAT+x+x+x2)*exp(-x)))/x3
return
end function
pure function B3_JGB(x)
 USE NUMBERS
 IMPLICIT DOUBLE PRECISION(A-H,O-Z)
real(FLOAT), intent(in) :: x
x2=x*x
x3=x2*x
x4=x3*x
xx=(6._FLOAT-6._FLOAT*x+3._FLOAT*x2-x3)*exp(x)/x4
yy=(6._FLOAT+6._FLOAT*x+3._FLOAT*x2+x3)*exp(-x)/x4
B3_JGB=xx-yy
return
end function


pure function B4_JGB(x)
 USE NUMBERS
 IMPLICIT DOUBLE PRECISION(A-H,O-Z)
real(FLOAT), intent(in) :: x
x2=x*x
x3=x2*x
x4=x3*x
x5=x4*x
xx=(24._FLOAT-24._FLOAT*x+12._FLOAT*x2-4._FLOAT*x3+x4)*exp(x)/x5
yy=(24._FLOAT+24._FLOAT*x+12._FLOAT*x2+4._FLOAT*x3+x4)*exp(-x)/x5
B4_JGB=xx-yy
return
end function
pure function B5_JGB(x)
 USE NUMBERS
 IMPLICIT DOUBLE PRECISION(A-H,O-Z)
real(FLOAT),intent(in) :: x
x2=x*x
x3=x2*x
x4=x3*x
x5=x4*x
x6=x5*x
xx=(120._FLOAT-120._FLOAT*x+60._FLOAT*x2-20._FLOAT*x3+5._FLOAT*x4-x5)*exp(x)/x6
yy=(120._FLOAT+120._FLOAT*x+60._FLOAT*x2+20._FLOAT*x3+5._FLOAT*x4+x5)*exp(-x)/x6
B5_JGB=xx-yy
return
end function
function B6_JGB(x)
 USE NUMBERS
 IMPLICIT DOUBLE PRECISION(A-H,O-Z)
real(FLOAT),intent(in) :: x
x2=x*x
x3=x2*x
x4=x3*x
x5=x4*x
x6=x5*x
x7=x6*x
xx=(720._FLOAT-720._FLOAT*x+360._FLOAT*x2-120._FLOAT*x3+30._FLOAT*x4-6._FLOAT*x5+x6)*exp(x)/x7
yy=(720._FLOAT+720._FLOAT*x+360._FLOAT*x2+120._FLOAT*x3+30._FLOAT*x4+6._FLOAT*x5+x6)*exp(-x)/x7
B6_JGB=xx-yy
return
end function
function bint(x,k)
 USE NUMBERS
 IMPLICIT DOUBLE PRECISION(A-H,O-Z)
real(FLOAT),intent(in) :: x
integer,intent(in) :: k
bint=0._FLOAT
if(abs(x).lt.1e-6_FLOAT)then
FAT=0._FLOAT
IFAT=0
do i=0,k
   bint=(1._FLOAT+(-1._FLOAT)**IFAT)/(FAT+1._FLOAT)
   FAT=FAT+1._FLOAT
   IFAT=IFAT+1
enddo
return
endif
FAT=K+1._FLOAT
do i=0,12
xx=1._FLOAT-(-1._FLOAT)**(k+i+1)
yy=fact_MPF(i)*FAT
FAT=FAT+1._FLOAT
bint=bint+xx/yy*(-x)**i
enddo
end function bint
function fact_MPF(N)
 USE NUMBERS
 IMPLICIT DOUBLE PRECISION(A-H,O-Z)
fact_MPF=1._FLOAT
FAC=2._FLOAT
do j=2,n
  fact_MPF=fact_MPF*FAC
  FAC=FAC+1._FLOAT
enddo
return 
end
subroutine gsovl(r,iat,jat,iz,xza,xzb,g)
 USE NUMBERS
 IMPLICIT DOUBLE PRECISION(A-H,O-Z)
 DIMENSION ishell(72),IZ(*)
data ishell/                 &                                                             
          1,1               &                                                             
          ,2,2,2,2,2,2,2,2, &                                                             
          3,3,3,3,3,3,3,3,  &                                                              
          54*3/
logical lsame
       za=xza
       zb=xzb
       na=iz(iat)
       nb=iz(jat)
ii=ishell(na)*ishell(nb)                                        
R05=R*0.5_FLOAT
ax=(za+zb)*R05
Fa=(za+zb)
bx=(zb-za)*R05
Fb=(zb-za)
lsame=.false.
if(za.eq.zb.OR.abs(za-zb).lt.0.1_FLOAT) then
lsame=.true.
  select case (ii)                                      
   case (1)                                             
     call JGB_g1s1s(za,zb,Fa,Fb,R,g,lsame)
   case (2)             
    if(ishell(na).lt.ishell(nb)) then
      call JGB_g2s1s(za,zb,Fa,Fb,R,g,.false.,lsame)
     else
      xx=za
      za=zb
      zb=xx
      call JGB_g2s1s(za,zb,Fa,Fb,R,g,.true.,lsame)
    endif
   case (4)                    
     call JGB_g2s2s(za,zb,Fa,Fb,R,g,lsame)
   case(3)                    
    if(ishell(na).lt.ishell(nb)) then
    call JGB_g1s3s(za,zb,Fa,Fb,R,g,.false.,lsame)
    else
      xx=za
      za=zb
      zb=xx
    call JGB_g1s3s(za,zb,Fa,Fb,R,g,.true.,lsame)
    endif
   case(6) 
    if(ishell(na).lt.ishell(nb)) then
    call JGB_g2s3s(za,zb,Fa,Fb,R,g,.false.,lsame)
    else
      xx=za
      za=zb
      zb=xx
    call JGB_g2s3s(za,zb,Fa,Fb,R,g,.true.,lsame)
    endif
   case(9)                                                         
    call JGB_g3s3s(za,zb,Fa,Fb,R,g,lsame)
   end select                                                                
else ! different elements
lsame=.false.                                                    
   select case (ii)                                                          
   case (1)                                                                  
     call JGB_g1s1s(za,zb,Fa,Fb,R,g,lsame)
   return
   case (2)  ! <1s|2s>                                                                  
    if(ishell(na).lt.ishell(nb)) then                                          
      call JGB_g2s1s(za,zb,Fa,Fb,R,g,.false.,lsame)
     else                                                                    
      xx=za                                                                  
      za=zb                                                                  
      zb=xx                                                                  
      call JGB_g2s1s(za,zb,Fa,Fb,R,g,.true.,lsame)
    endif                                                                    
   case (4) ! <2s|2s>                                                        
      call JGB_g2s2s(za,zb,Fa,Fb,R,g,lsame)
   case(3)  ! <1s|3s> + <3s|1s>                                              
    if(ishell(na).lt.ishell(nb)) then                                          
    call JGB_g1s3s(za,zb,Fa,Fb,R,g,.false.,lsame)
    else                                                                                  
      xx=za                                                                               
      za=zb                                                                               
      zb=xx                                                                               
    call JGB_g1s3s(za,zb,Fa,Fb,R,g,.true.,lsame)
    endif                                                                                 
   case(6)  ! <2s|3s> + <3s|2s>                                                           
    if(ishell(na).lt.ishell(nb)) then                                                       
    call JGB_g2s3s(za,zb,Fa,Fb,R,g,.false.,lsame)
    else                                                                                                       
      xx=za                                                                                                    
      za=zb                                                                                                    
      zb=xx                                                                                                    
    call JGB_g2s3s(za,zb,Fa,Fb,R,g,.true.,lsame)
    endif                                                                                                      
    case(9) ! <3s|3>                                                                                           
    call JGB_g3s3s(za,zb,Fa,Fb,R,g,lsame)
   end select                                                                                                  
endif                                                                                                          
return
end subroutine gsovl
subroutine JGB_g1s1s(za,zb,Fa,Fb,R,g,sameElement)
 USE NUMBERS
 IMPLICIT DOUBLE PRECISION(A-H,O-Z)
 PARAMETER(THIRD=-1._FLOAT/3._FLOAT)
logical sameElement
if(sameElement) then
  t1 = za ** 2
  t3 = zb ** 2
  t5 = t1 * za * t3 * zb
  t6 = R ** 2
  t7 = t6 ** 2
  t10 = Fa * R
  t14 = exp(-0.5_FLOAT*t10)
  t17 = sqrt(t5 * t7 * t6)
  g = THIRD * t5 * t7 / Fa * (2._FLOAT + t10) * t14 / t17
  return
else
  t1 = za ** 2
  t3 = zb ** 2
  t5 = t1 * za * t3 * zb
  t6 = Fb ** 2
  t7 = Fb * R
  t8 = 0.5_FLOAT * t7
  t9 = exp(t8)
  t12 = exp(-t8)
  t15 = t6 * Fa
  t22 = Fa ** 2
  t23 = t22 * t9
  t27 = t22 * t12
  t31 = t6 * Fb
  t32 = R * t31
  t37 = t22 * Fa
  t38 = R * t37
  t43 = R ** 2
  t44 = t43 * t31
  t51 = t43 * t37
  t56 = 4._FLOAT * t6 * t9 - 4._FLOAT * t6 * t12 + 2._FLOAT * t15 * R * t9 -          &
  2._FLOAT * t15 * R * t12 - 4._FLOAT * t23 + 2._FLOAT * t23 * t7 + 4._FLOAT * t27       &
  + 2._FLOAT * t27 * t7 - 2._FLOAT * t32 * t9 - 2._FLOAT * t32 * t12 - 2._FLOAT * t38 *  &
  t9 + 2._FLOAT * t38 * t12 - 1._FLOAT * t44 * Fa * t9 - 1._FLOAT                     &
  * t44 * Fa * t12 + t51 * t9 * Fb + t51 * t12 * Fb
  t61 = exp(-0.5_FLOAT * Fa * R)
  t62 = t43 ** 2
  t65 = sqrt(t5 * t62 * t43)
  g = -2._FLOAT * t5 * R * t56 * t61 / t65 / t31 / t37
  return
endif
end subroutine JGB_g1s1s

subroutine JGB_g2s1s(za,zb,Fa,Fb,R,g,switch,lsame)
 USE NUMBERS
 IMPLICIT DOUBLE PRECISION(A-H,O-Z)
 PARAMETER(VENTIQ=1._FLOAT/24._FLOAT)
logical switch
logical lsame
anorm=VENTIQ*sqrt(za**3*zb**5*3_FLOAT)
if(switch) then
Fb=-Fb
endif
if(lsame) then
      t1 = Fa * R
      t3 = exp(-0.5_FLOAT * t1)
      t6 = Fa ** 2
      t7 = R ** 2
      g = -0.1000000000e-8_FLOAT* R * t3 * (0.5333333380e10_FLOAT+ 0.2666666670e10_FLOAT&
      * t1 + 0.1333333333e10_FLOAT * t6 * t7) / t6
      g=g*anorm
else
      t3 = exp(-0.5_FLOAT * Fa * R)
      t4 = Fa ** 2
      t5 = t4 * Fa
      t6 = Fb * R
      t7 = 0.5_FLOAT * t6
      t8 = exp(t7)
      t9 = t5 * t8
      t11 = Fb ** 2
      t12 = t11 * Fa
      t15 = exp(-t7)
      t18 = t4 ** 2
      t19 = R * t18
      t22 = t11 ** 2
      t29 = Fb * t4
      t36 = R ** 2
      t37 = t36 * t18
      t44 = t36 * R
      t48 = -12._FLOAT * t9 + 4._FLOAT * t12 * t8 - 4._FLOAT * t12 * t15 &
      - 6._FLOAT * t19 * t8 - 6._FLOAT * t22 * t8 * R - 6._FLOAT * t22 * t15 &
      * R + 4._FLOAT * t29 * t15 - 4._FLOAT * t29 * t8 + 6._FLOAT * t19 * t15 &
      + 2._FLOAT * t37 * t8 * Fb + 4._FLOAT * t37 * t15 * Fb + t44 * t18 * t15 * t11
      t49 = t5 * t15
      t51 = t11 * Fb
      t58 = t51 * Fa
      t59 = R * t8
      t76 = t36 * t15
      t79 = t22 * Fa
      t87 = 12._FLOAT * t49 - 12._FLOAT * t51 * t15 - 1._FLOAT * t22 * t4 * t15 * t44 &
      + 4._FLOAT * t58 * t59 - 8._FLOAT * t58 * R * t15 + 4._FLOAT * t9 * t6 + 8._FLOAT * &
      t49 * t6 + 2._FLOAT * t49 * t11 * t36 + 4._FLOAT * t11 * t4 * t59 - 2._FLOAT * t51 &
      * t4 * t76 - 2._FLOAT * t79 * t36 * t8 - 4._FLOAT * t79 * t76 + 12._FLOAT * t51 * t8
      g = -16._FLOAT * t3 * (t48 + t87) / t36 / t22 / t18
      g=g*anorm
endif
return
end subroutine JGB_g2s1s

subroutine JGB_g2s2s(za,zb,Fa,Fb,R,g,SameElement)
 USE NUMBERS
 IMPLICIT DOUBLE PRECISION(A-H,O-Z)
 PARAMETER(QUARA8=1._FLOAT/48._FLOAT)
logical SameElement
 anorm=SQRT((ZA*ZB)**5)*QUARA8
if(SameElement) then
      t2 = R*R
      t5 = Fa*FA
      t9 = t5 * Fa
      t10 = t2*T2
      t16 = exp(-Fa * R*.5_FLOAT)
      g = (-42.66666666_FLOAT * R - 21.33333333_FLOAT * Fa * t2 - 2.133333333_FLOAT &
          * t5 * t2 * R - 1.066666666_FLOAT * t9 * t10) * t16 / t9
      g=g*anorm
return
else
      t1 = R*R
      t3 = 384._FLOAT * t1 * Fb
      t4 = t1 * R
      t5 = Fb*FB
      t7 = 64._FLOAT * t4 * t5
      t8 = 768._FLOAT * R
      t10 = Fa*FA
      t11 = t10*T10
      t12 = t11 * Fa
      t14 = Fb * R
      t15 = 768._FLOAT * t14
      t17 = 128._FLOAT * t5 * t1
      t21 = 256._FLOAT * t5 * R
      t22 = t5 * Fb
      t24 = 128._FLOAT * t22 * t1
      t26 = t10 * Fa
      t28 = t5*T5
      t30 = 128._FLOAT * t1 * t28
      t32 = 256._FLOAT * t22 * R
      t33 = 512._FLOAT * t5
      t34 = t28 * Fb
      t36 = 64._FLOAT * t4 * t34
      t40 = 768._FLOAT * t28 * R
      t42 = 384._FLOAT * t1 * t34
      t45 = 1536._FLOAT * t28
      t47 = 768._FLOAT * t34 * R
      t51 = exp(-0.5_FLOAT * Fa * R)
      t53 = 0.5_FLOAT * t14
      t54 = exp(-t53)
      t68 = exp(t53)
      g = (((t3 + t7 + t8) * t12 + (1536._FLOAT + t15 + t17) * t11 + &
      (-t21 - t24) * t26 + (t30 - t32 - t33 + t36) * t10 + (t40 + t42) *&
      Fa + t45 + t47) * t51 * t54 + ((t3 - t8 - t7) * t12 + (-1536._FLOAT &
      + t15 - t17) * t11 + (-t24 + t21) * t26 + (-t30 + t33 - t32 &
      + t36) * t10 + (-t40 + t42) * Fa + t47 - t45) * t51 * t68) / t1 /  &
      t12 / t34
      g=g*anorm
return
endif
end subroutine JGB_g2s2s

subroutine JGB_g1s3s(za,zb,Fa,Fb,R,g,switch,lsame)
 USE NUMBERS
 IMPLICIT DOUBLE PRECISION(A-H,O-Z)
 PARAMETER(RAD148=0.0360843918243516082) ! SQRT(1/3)/16)
logical switch,lsame
if(switch)Fb=-Fb
Anorm=SQRT((ZA**3)*(ZB**7)/7.5_FLOAT)*RAD148

if(lsame) then
  t1 = Fa * R
  t3 = exp(-0.5_FLOAT * t1)
  t4 = R ** 2
  g = -1.6_FLOAT * t3 * t4 * R * (2._FLOAT + t1) / Fa
  g=g*Anorm
else
      t3 = exp(-0.5_FLOAT * Fa * R)
      t4 = Fb ** 2
      t5 = t4 ** 2
      t6 = t5 * Fb
      t7 = t6 * Fa
      t8 = R ** 2
      t9 = Fb * R
      t10 = 0.5_FLOAT * t9
      t11 = exp(t10)
      t15 = exp(-t10)
      t16 = t8 * t15
      t19 = Fa ** 2
      t21 = t8 * R
      t22 = t21 * t15
      t25 = t19 * Fa
      t27 = t8 ** 2
      t31 = t19 ** 2
      t32 = t31 * Fa
      t33 = t8 * t32
      t45 = t4 * Fb
      t48 = t31 * t15
      t55 = t4 * t25
      t56 = t11 * R
      t59 = t15 * R
      t62 = t5 * Fa
      t73 = -6._FLOAT * t7 * t8 * t11 - 18._FLOAT * t7 * t16 - 6._FLOAT * t6 * t19 &
      * t22 - 1._FLOAT * t6 * t25 * t27 * t15 + 6._FLOAT * t33 * t11 * Fb + 18._FLOAT &
      * t33 * t15 * Fb + 6._FLOAT * t21 * t32 * t15 * t4 + t27 * t32* t15 * t45 &
      + 2._FLOAT * t48 * t45 * t21 + 12._FLOAT * t48 * t4 * t8 + 12._FLOAT * t55 * t56 &
      + 12._FLOAT * t55 * t59 + 12._FLOAT * t62 * t56 - 36._FLOAT * t62 * t59 - 12._FLOAT &
      * t5 * t19 * t16 - 2._FLOAT * t5 * t25 * t22
      t74 = t31 * t11
      t79 = t45 * t19
      t92 = R * t32
      t95 = t45 * Fa
      t100 = Fb * t25
      t111 = 12._FLOAT * t74 * t9 + 36._FLOAT * t48 * t9 + 12._FLOAT * t79 * t56 - 12._FLOAT  &
      * t79 * t59 + 48._FLOAT * t5 * t11 - 24._FLOAT * t6 * t11 * R - 24._FLOAT * t6 * t15 &
      * R - 24._FLOAT * t92 * t11 + 24._FLOAT * t95 * t11 - 24._FLOAT * t95 * t15 + 24._FLOAT &
      * t100 * t15 + 24._FLOAT * t92 * t15 - 24._FLOAT * t100 * t11 - 48._FLOAT * t5 * t15 &
      - 48._FLOAT * t74 + 48._FLOAT * t48
      g = -32._FLOAT * t3 * (t73 + t111) / t8 / t6 / t32
      g=g*Anorm
endif
return
end subroutine JGB_g1s3s

subroutine JGB_g2s3s(za,zb,Fa,Fb,R,g,switch,lsame)
 USE NUMBERS
 IMPLICIT DOUBLE PRECISION(A-H,O-Z)
PARAMETER(SETTE5=1._FLOAT/7.5_FLOAT,OVANTA6=1._FLOAT/96._FLOAT)
logical switch,lsame
Anorm=sqrt((za**5)*(zb**7)*SETTE5)*OVANTA6

if(switch) Fb=-Fb

if(lsame) then
      t1 = Fa * R
      t3 = exp(-0.5_FLOAT * t1)
      t6 = Fa*FA
      t7 = R*R
      t14 = t6*T6
      t15 = t7*T7
      g = -0.2E-8_FLOAT * R * t3 * (0.128E12_FLOAT + 0.64E11_FLOAT &
      * t1 + 0.128E11_FLOAT * t6 * t7 + 0.106666666666667E10_FLOAT * t6 * Fa * t7 &
      * R + 0.533333333333333E9_FLOAT * t14 * t15) / t14
      g=g*Anorm
else
      t3 = exp(-0.5_FLOAT * Fa * R)
      t4 = Fb*FB
      t5 = t4*T4
      t6 = Fa*FA
      t7 = t6 * Fa
      t8 = t5 * t7
      t9 = R*R
      t11 = 0.5_FLOAT * Fb * R
      t12 = exp(t11)
      t13 = t9 * t12
      t16 = t6*T6
      t17 = t16 * Fa
!      t18 = 1._FLOAT/T12
      t18 = exp(-t11)
      t21 = t5 * Fb
      t28 = t9 * t18
      t32 = t9 * R
      t33 = t32 * t18
      t36 = t5 * t4
      t38 = t9*T9
      t39 = t38 * t18
      t41 = t21 * Fa
      t42 = R * t12
      t45 = t16 * t6
      t46 = t4 * Fb
      t49 = t46 * t16
      t52 = -6._FLOAT * t8 * t13 + 120._FLOAT * t17 * t18 + 120._FLOAT * t21 * t18 &
       - 120._FLOAT * t17 * t12 - 120._FLOAT * t21 * t12 - 6._FLOAT * t8 * t28 - 2._FLOAT &
      * t5 * t16 * t33 + t36 * t7 * t39 - 48._FLOAT * t41 * t42 + t45 * t46 * t39 - 6._FLOAT * t49 * t13
      t54 = R * t18
      t60 = t46 * t6
      t63 = Fb * t16
      t66 = t5 * Fa
      t69 = t4 * t7
      t72 = t36 * t6
      t75 = t32 * t12
      t78 = Fb * t9
      t84 = Fb * t17
      t87 = -24._FLOAT * t46 * t7 * t54 - 24._FLOAT * t5 * t6 * t42 + 24._FLOAT *&
       t60 * t12 + 24._FLOAT * t63 * t18 - 24._FLOAT * t66 * t12 - 24._FLOAT * t69 &
      * t18 + 9._FLOAT * t72 * t33 + 3._FLOAT * t72 * t75 + 24._FLOAT * t78 * t45 &
      * t12 - 6._FLOAT * t49 * t28 + 48._FLOAT * t84 * t42
      t102 = t21 * t6
      t105 = t4 * t17
      t113 = t45 * t4
      t118 = 72._FLOAT * t84 * t54 + 72._FLOAT * t41 * t54 + 36._FLOAT * t78 * t45 &
      * t18 + 2._FLOAT * t46 * t17 * t33 + 24._FLOAT * t4 * t16 * t42 - 6._FLOAT   &
      * t102 * t13 - 6._FLOAT * t105 * t13 + 18._FLOAT * t105 * t28 + 2._FLOAT &
      * t21 * t7 * t33 - 3._FLOAT * t113 * t75 + 9._FLOAT * t113 * t33
      t121 = t36 * Fa
      t130 = R * t45
      t145 = 18._FLOAT * t102 * t28 + 24._FLOAT * t121 * t13 + 36._FLOAT * t121 * &
       t28 - 24._FLOAT * t60 * t18 - 24._FLOAT * t63 * t12 + 60._FLOAT * t130 * t18 &
      + 60._FLOAT * t36 * t18 * R + 24._FLOAT * t69 * t12 + 60._FLOAT * t36 * t12 * &
      R - 60._FLOAT * t130 * t12 + 24._FLOAT * t66 * t18
      g = 128._FLOAT * t3 * (t52 + t87 + t118 + t145) / t9 / t36 / t45
  g=g*Anorm
endif
return
end subroutine JGB_g2s3s

subroutine JGB_g3s3s(za,zb,Fa,Fb,R,g,SameElement)
 USE NUMBERS
 IMPLICIT DOUBLE PRECISION(A-H,O-Z)
 PARAMETER(TQUA=1._FLOAT/1440._FLOAT)
logical SameElement
Anorm=sqrt((ZA*ZB)**7)*TQUA
if(SameElement) then
      t1 = Fa * R
      t3 = exp(-0.5_FLOAT * t1)
      t5 = Fa ** 2
      t6 = t5 ** 2
      t7 = t6 * Fa
      t8 = R ** 2
      t9 = t8 ** 2
      g = -0.2E-8_FLOAT* t3 * R * (0.457142857E9_FLOAT* t7 * t9 &
      * R + 0.768E12_FLOAT* t1 + 0.1536E12_FLOAT* t5 * t8 &
      + 0.128E11_FLOAT* t5 * Fa * t8 * R + 0.914285715E9_FLOAT* t6 * t9 + 0.1536E13_FLOAT) / t7
      g=g*Anorm
return
else
      t3 = exp(-0.5_FLOAT * Fa * R)
      t4 = Fa ** 2
      t5 = t4 ** 2
      t6 = t5 * t4
      t7 = Fb * R
      t8 = 0.5_FLOAT * t7
      t9 = exp(-t8)
      t10 = t6 * t9
      t13 = Fb ** 2
      t14 = t13 * Fb
      t15 = t13 ** 2
      t16 = t15 * t14
      t17 = R ** 2
      t18 = t17 * R
      t19 = t16 * t18
      t23 = exp(t8)
      t24 = t6 * t23
      t27 = t5 * Fa
      t28 = t27 * t13
      t29 = R * t23
      t32 = t6 * t13
      t33 = t17 * t9
      t36 = t15 * Fb
      t37 = t4 * t36
      t38 = t9 * R
      t43 = t17 * t23
      t46 = t4 * Fa
      t47 = t5 * t46
      t48 = t47 * t18
      t52 = t47 * t17
      t65 = 120._FLOAT * t10 * t7 - 12._FLOAT * t19 * t4 * t9 + 120._FLOAT &
      * t24 * t7 + 24._FLOAT * t28 * t29 + 24._FLOAT * t32 * t33 + 24._FLOAT * t37 &
      * t38 - 24._FLOAT * t28 * t38 - 24._FLOAT * t32 * t43 - 12._FLOAT * t48 * t13 &
      * t23 + 60._FLOAT * t52 * t23 * Fb + 12._FLOAT * t48 * t13 * t9 + 60._FLOAT &
      * t52 * t9 * Fb - 12._FLOAT * t19 * t4 * t23
      t66 = t17 ** 2
      t67 = t16 * t66
      t74 = t27 * t14
      t77 = t6 * t14
      t78 = t18 * t23
      t81 = t46 * t15
      t86 = t27 * t15
      t89 = t5 * t36
      t90 = t18 * t9
      t97 = t46 * t36
      t104 = -1._FLOAT * t67 * t46 * t9 - 1._FLOAT * t67 * t46 * t23 - 12._FLOAT &
      * t74 * t43 + 2._FLOAT * t77 * t78 - 24._FLOAT * t81 * t29 + 24._FLOAT * t81 &
      * t38 + 2._FLOAT * t86 * t78 + 2._FLOAT * t89 * t90 - 2._FLOAT * t86 * t90 &
      + 24._FLOAT * t37 * t29 + 12._FLOAT * t97 * t33 + 2._FLOAT * t89 * t78 - 12._FLOAT * t74 * t33
      t108 = t5 * t14
      t111 = t15 * t13
      t112 = t111 * t4
      t117 = t111 * t46
      t122 = t111 * Fa
      t129 = t4 * t15
      t132 = t47 * R
      t139 = 2._FLOAT * t77 * t90 - 24._FLOAT * t108 * t38 + 24._FLOAT * t112 * t43 &
      - 24._FLOAT * t112 * t33 + 2._FLOAT * t117 * t78 - 2._FLOAT * t117 * t90 + 120._FLOAT &
      * t122 * t29 - 120._FLOAT * t122 * t38 + 12._FLOAT * t97 * t43 - 48._FLOAT * t129 &
      * t23 + 120._FLOAT * t132 * t9 - 120._FLOAT * t132 * t23 + 240._FLOAT * t111 * t23
      t140 = t47 * t66
      t145 = t16 * R
      t150 = t16 * t17
      t160 = t5 * t13
      t170 = t140 * t14 * t23 + t140 * t14 * t9 - 120._FLOAT * t145 * t9 - 24._FLOAT &
      * t108 * t29 - 60._FLOAT * t150 * Fa * t23 - 240._FLOAT * t111 * t9 - 240._FLOAT &
      * t24 + 240._FLOAT * t10 + 48._FLOAT * t129 * t9 - 48._FLOAT * t160 * t9 + 48._FLOAT &
      * t160 * t23 - 120._FLOAT * t145 * t23 - 60._FLOAT * t150 * Fa * t9
      g = -768._FLOAT * t3 * (t65 + t104 + t139 + t170) / t17 / t47 / t16
      g=g*Anorm
return
endif
end subroutine JGB_g3s3s

      subroutine lower_case(word)
      character (len=*) , intent(in out) :: word
      integer :: i,ic,nlen
      nlen = len(word)
      do i=1,nlen
      ic = ichar(word(i:i))
      if (ic >= 65 .and. ic < 90) word(i:i) = char(ic+32)
      end do
      end subroutine lower_case

subroutine charXsplit(s,wx,x)
implicit none
integer i,k,x
character(LEN=80),intent(in) :: s
character(LEN=80),intent(out) ::wx
character(LEN=80) :: w(20),a,aa
aa=adjustl(s)
i=0
do while (index(aa,' ').ne.1)
  i=i+1
  a=aa
  k=index(a,' ')
  w(i)=trim(a(:(k-1)))
  aa=adjustl(trim(a((k+1):)))
  if(i.gt.50)CALL ERRVRS(0,'charXsplit','string split error: subroutine charXsplit')
enddo
wx=w(x)
return
end subroutine

subroutine setzet(eta,za,zb)    
 USE NUMBERS
 IMPLICIT DOUBLE PRECISION(A-H,O-Z)
 PARAMETER(THIRD=1._FLOAT/3._FLOAT,MAX_PARA=36)
 DIMENSION za(MAX_PARA),zb(MAX_PARA)
 REAL(FLOAT),DIMENSION(MAX_PARA) :: ZS,ZP,ZD
data ZS /1.2000_FLOAT,1.6469_FLOAT,0.6534_FLOAT,1.0365_FLOAT,&
1.3990_FLOAT,1.7210_FLOAT,2.0348_FLOAT,2.2399_FLOAT,2.5644_FLOAT,2.8812_FLOAT,&
0.8675_FLOAT,1.1935_FLOAT,1.5143_FLOAT,1.7580_FLOAT,1.9860_FLOAT,2.1362_FLOAT,2.3617_FLOAT,2.5796_FLOAT,0.9362_FLOAT,1.2112_FLOAT,&
1.2870_FLOAT,1.3416_FLOAT,1.3570_FLOAT,1.3804_FLOAT,1.4761_FLOAT,1.5465_FLOAT,1.5650_FLOAT,1.5532_FLOAT,1.5781_FLOAT,1.7778_FLOAT,&
2.0675_FLOAT,2.2702_FLOAT,2.4546_FLOAT,2.5680_FLOAT,2.7523_FLOAT,2.9299_FLOAT/
data ZP /0.0000_FLOAT,0.0000_FLOAT,0.5305_FLOAT,0.8994_FLOAT,&
1.2685_FLOAT,1.6105_FLOAT,1.9398_FLOAT,2.0477_FLOAT,2.4022_FLOAT,2.7421_FLOAT,&
0.6148_FLOAT,0.8809_FLOAT,1.1660_FLOAT,1.4337_FLOAT,1.6755_FLOAT,1.7721_FLOAT,2.0176_FLOAT,2.2501_FLOAT,0.6914_FLOAT,0.9329_FLOAT,&
0.9828_FLOAT,1.0104_FLOAT,0.9947_FLOAT,0.9784_FLOAT,1.0641_FLOAT,1.1114_FLOAT,1.1001_FLOAT,1.0594_FLOAT,1.0527_FLOAT,1.2448_FLOAT,&
1.5073_FLOAT,1.7680_FLOAT,1.9819_FLOAT,2.0548_FLOAT,2.2652_FLOAT,2.4617_FLOAT/
data ZD /&
0.0000_FLOAT,0.0000_FLOAT,0.0000_FLOAT,0.0000_FLOAT,0.0000_FLOAT,0.0000_FLOAT,0.0000_FLOAT,0.0000_FLOAT,0.0000_FLOAT,0.0000_FLOAT,&
0.0000_FLOAT,0.0000_FLOAT,0.0000_FLOAT,0.0000_FLOAT,0.0000_FLOAT,0.0000_FLOAT,0.0000_FLOAT,0.0000_FLOAT,0.0000_FLOAT,0.0000_FLOAT,&
2.4341_FLOAT,2.6439_FLOAT,2.7809_FLOAT,2.9775_FLOAT,3.2208_FLOAT,3.4537_FLOAT,3.6023_FLOAT,3.7017_FLOAT,3.8962_FLOAT,2.0477_FLOAT,&
2.4022_FLOAT,2.7421_FLOAT,0.6148_FLOAT,0.8809_FLOAT,1.1660_FLOAT,1.4337_FLOAT/
  do i=1,MAX_PARA
select case (i)
  case(:2)
    za(i)=ZS(i)
  case(3:20,31:)
    za(i)=( ZS(i)+ZP(i) )*.5_FLOAT
  case(21:30)
    za(i)=( ZS(i)+ZP(i)+ZD(i) )*THIRD
end select 
  enddo
  za=za*eta
  zb=za
 return
end subroutine setzet
                                                                                 
      FUNCTION READAA(A,ISTART,IEND,IEND2)                                       
 USE NUMBERS
 IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      CHARACTER*(*) A                                                            
      NINE=ICHAR('9')                                                            
      IZERO=ICHAR('0')                                                           
      MINUS=ICHAR('-')                                                           
      IDOT=ICHAR('.')                                                            
      ND=ICHAR('D')                                                              
      NE=ICHAR('E')                                                              
      IBL=ICHAR(' ')                                                             
      IEND=0                                                                     
      IEND2=0                                                                    
      IDIG=0                                                                     
      C1=0._FLOAT                                                                       
      C2=0._FLOAT
      ONE=1._FLOAT
      X = 1._FLOAT
      NL=LEN(A)                                                                  
      DO 10 J=ISTART,NL-1                                                        
         N=ICHAR(A(J:J))                                                          
         M=ICHAR(A(J+1:J+1))                                                      
         IF(N.LE.NINE.AND.N.GE.IZERO .OR.N.EQ.IDOT)GOTO 20                        
         IF(N.EQ.MINUS.AND.(M.LE.NINE.AND.M.GE.IZERO.OR. M.EQ.IDOT)) GOTO 20

   10 CONTINUE
      READAA=0._FLOAT
      RETURN
   20 CONTINUE
      IEND=J
      DO 30 I=J,NL
         N=ICHAR(A(I:I))
         IF(N.LE.NINE.AND.N.GE.IZERO) THEN
            IDIG=IDIG+1
            IF (IDIG.GT.10) GOTO 60
            C1=C1*10._FLOAT+N-IZERO
         ELSEIF(N.EQ.MINUS.AND.I.EQ.J) THEN
            ONE=-1._FLOAT
         ELSEIF(N.EQ.IDOT) THEN
            GOTO 40
         ELSE
            GOTO 60
         ENDIF
   30 CONTINUE
   40 CONTINUE
      IDIG=0
      DO 50 II=I+1,NL
         N=ICHAR(A(II:II))
         IF(N.LE.NINE.AND.N.GE.IZERO) THEN
            IDIG=IDIG+1
            IF (IDIG.GT.10) GOTO 60
            C2=C2*10._FLOAT+N-IZERO
            X = X*.1_FLOAT
         ELSEIF(N.EQ.MINUS.AND.II.EQ.I) THEN
            X=-X
         ELSE
            GOTO 60
         ENDIF
   50 CONTINUE
   60 CONTINUE
      READAA= ONE * ( C1 + C2 * X)
      DO 55 J=IEND,NL
         N=ICHAR(A(J:J))
         IEND2=J
         IF(N.EQ.IBL)RETURN
   55 IF(N.EQ.ND .OR. N.EQ.NE)GOTO 57
      RETURN
   57 C1=0._FLOAT
      ONE=1._FLOAT
      DO 31 I=J+1,NL
         N=ICHAR(A(I:I))
         IEND2=I
         IF(N.EQ.IBL)GOTO 70
         IF(N.LE.NINE.AND.N.GE.IZERO) C1=C1*10._FLOAT+N-IZERO
         IF(N.EQ.MINUS)ONE=-1._FLOAT
   31 CONTINUE
   61 CONTINUE
   70 READAA=READAA*10._FLOAT**(ONE*C1)
      RETURN
      END

FUNCTION ESYM(I)
 USE NUMBERS
 IMPLICIT DOUBLE PRECISION(A-H,O-Z)
 PARAMETER(MAX_ELEM=94)
CHARACTER(LEN=2) :: ELEMNT(MAX_ELEM),ESYM 
DATA ELEMNT/'h ','he',                                           &
  'li','be','b ','c ','n ','o ','f ','ne',                       &
  'na','mg','al','si','p ','s ','cl','ar',                       &
  'k ','ca','sc','ti','v ','cr','mn','fe','co','ni','cu',        &
  'zn','ga','ge','as','se','br','kr',                            &
  'rb','sr','y ','zr','nb','mo','tc','ru','rh','pd','ag',        &
  'cd','in','sn','sb','te','i ','xe',                            &
  'cs','ba','la','ce','pr','nd','pm','sm','eu','gd','tb','dy',   &
  'ho','er','tm','yb','lu','hf','ta','w ','re','os','ir','pt',   &
  'au','hg','tl','pb','bi','po','at','rn',                       &
  'fr','ra','ac','th','pa','u ','np','pu'/                       
  ESYM=ELEMNT(I)
  RETURN
  END  

character(80) function getlevel(method)
implicit none
character(*) method
character(80) bas
integer i
  i=scan(method,'/')
  bas=trim(method((i+1):))
getlevel='basis:' // trim(bas)

select case(bas)
case('minis')
 getlevel='basis: MINIS'
case('minix')
 getlevel='basis: MINIX'
case('sv')
 getlevel='basis: SV (Ahlrichs)'
case('sv(p)')
 getlevel='basis: def2-SV(P)'
case('svp')
 getlevel='basis: def2-SVP'
case('svx')
 getlevel='basis: def2-SV(P/h,c)'
case('svp_old')
 getlevel='basis: def2-SVP (old parameters)'
case('631gd')
 getlevel='basis: 6-31G(d)'
case('lanl')
 getlevel='basis: 6-31G(d)+LANL2DZ(Sc-Zn)'
case('tz')
 getlevel='basis: def2-TZVP'
end select
end function

SUBROUTINE JGB_criteria(r_cutoff,Hlat,Itau_max)
 USE NUMBERS
 IMPLICIT DOUBLE PRECISION(A-H,O-Z)
DIMENSION Hlat(3,3),Anorm1(3),Anorm2(3),Anorm3(3),Itau_max(3)
call crossprodukt(Hlat(:,2),Hlat(:,3),Anorm1)
call crossprodukt(Hlat(:,3),Hlat(:,1),Anorm2)
call crossprodukt(Hlat(:,1),Hlat(:,2),Anorm3)
Anorm1(1:3)=Anorm1(1:3)/SQRT(ANORM1(1)*ANORM1(1)+ANORM1(2)*ANORM1(2)+ANORM1(3)*ANORM1(3))
Anorm2(1:3)=Anorm2(1:3)/SQRT(ANORM2(1)*ANORM2(1)+ANORM2(2)*ANORM2(2)+ANORM2(3)*ANORM2(3))
Anorm3(1:3)=Anorm3(1:3)/SQRT(ANORM3(1)*ANORM3(1)+ANORM3(2)*ANORM3(2)+ANORM3(3)*ANORM3(3))
cos10=SUM(Anorm1(1:3)*Hlat(1:3,1))
cos21=SUM(Anorm2(1:3)*Hlat(1:3,2))
cos32=SUM(Anorm3(1:3)*Hlat(1:3,3))
Itau_max(1)= ceiling(abs(r_cutoff/cos10))
Itau_max(2)= ceiling(abs(r_cutoff/cos21))
Itau_max(3)= ceiling(abs(r_cutoff/cos32))
END SUBROUTINE JGB_CRITERIA

SUBROUTINE crossprodukt(A,B,C)
 USE NUMBERS
 IMPLICIT DOUBLE PRECISION(A-H,O-Z)
 DIMENSION A(3),B(3),C(3)
X=A(2)*B(3)-B(2)*A(3)
Y=A(3)*B(1)-B(3)*A(1)
Z=A(1)*B(2)-B(1)*A(2)
C=(/X,Y,Z/)
END SUBROUTINE crossprodukt

subroutine abc2xyz(n,xyz,abc,Hlat)
 USE NUMBERS
 IMPLICIT DOUBLE PRECISION(A-H,O-Z)
 DIMENSION xyz(3,n),abc(3,n), Hlat(3,3),Hlatt(3,3)
xyz=0._FLOAT
do i=1,3
   do j=1,3
      Hlatt(i,j)=Hlat(j,i)
   end do
end do
!Hlatt=Hlat
do i=1,n
   do j=1,3
      xyz(j,i)=Hlatt(j,1)*abc(1,i)+Hlatt(j,2)*abc(2,i)+Hlatt(j,3)*abc(3,i)
   end do
end do
end subroutine abc2xyz

subroutine xyz2abc(n,xyz,abc,Hlat)
 USE NUMBERS
 IMPLICIT DOUBLE PRECISION(A-H,O-Z)
 PARAMETER(TOL=1E-14_FLOAT)
 DIMENSION xyz(3,n),abc(3,n),Hlat(3,3),Hlatt(3,3),olat(3,3)
olat=0._FLOAT
abc=0._FLOAT
do i=1,3
   do j=1,3
      Hlatt(i,j)=Hlat(j,i)
   end do
end do

detlat =Hlatt(1,1)*Hlatt(2,2)*Hlatt(3,3)+Hlatt(1,2)*Hlatt(2,3)*Hlatt(3,1)+Hlatt(1,3)*Hlatt(2,1)*Hlatt(3,2)
detlat=detlat-Hlatt(3,1)*Hlatt(2,2)*Hlatt(1,3)-Hlatt(3,2)*Hlatt(2,3)*Hlatt(1,1)-Hlatt(3,3)*Hlatt(2,1)*Hlatt(1,1)

if(abs(detlat).lt.TOL)CALL ERRVRS(0,'XYZ1ABC','singular cell matrix')
olat(1,1)=Hlatt(2,2)*Hlatt(3,3)-Hlatt(3,2)*Hlatt(2,3)
olat(1,2)=Hlatt(1,3)*Hlatt(3,2)-Hlatt(3,3)*Hlatt(1,2)
olat(1,3)=Hlatt(1,2)*Hlatt(2,3)-Hlatt(2,2)*Hlatt(1,3)
olat(2,1)=Hlatt(2,3)*Hlatt(3,1)-Hlatt(3,3)*Hlatt(2,1)
olat(2,2)=Hlatt(1,1)*Hlatt(3,3)-Hlatt(3,1)*Hlatt(1,3)
olat(2,3)=Hlatt(1,3)*Hlatt(2,1)-Hlatt(2,3)*Hlatt(1,1)
olat(3,1)=Hlatt(2,1)*Hlatt(3,2)-Hlatt(3,1)*Hlatt(2,2)
olat(3,2)=Hlatt(1,2)*Hlatt(3,1)-Hlatt(3,2)*Hlatt(1,1)
olat(3,3)=Hlatt(1,1)*Hlatt(2,2)-Hlatt(2,1)*Hlatt(1,2)
 DETLAT=1._FLOAT/DETLAT
do i=1,3
   do j=1,3
      olat(i,j)=olat(i,j)*detlat
   end do
end do

do i=1,n
   do j=1,3
      abc(j,i)=olat(j,1)*xyz(1,i)+olat(j,2)*xyz(2,i)+olat(j,3)*xyz(3,i)
   end do
end do
end subroutine xyz2abc

subroutine srb_egrad2(xyz,iz,Hlat,n,energy,g,cellgrad,grad,rscal,qscal,echo)
USE NUMBERS
USE PARINF_MODULE
USE PARAL1_MODULE
use gcp_module 
USE MEMORY_USE
IMPLICIT DOUBLE PRECISION(A-H,O-Z)
parameter (autokcal=627.509541_FLOAT,max_elem=94,max_para=36)
DIMENSION XYZ(3,N),G(3,N),Hlat(3,3),IZ(N),ITAU_MAX(3),pp(max_elem),co(n),Hlat_1(3,3)
real(FLOAT) :: cellgrad(3,3), stress(3,3)
REAL(FLOAT),DIMENSION(:,:),ALLOCATABLE :: r0ab
integer version,iz
logical echo,grad
character*2  esym
AUTOANG=PAR(32)

thrR=30.0_FLOAT            ! X bohr
thrE=epsilon(1._FLOAT)
 cellgrad=0.0_FLOAT
 stress=0.0_FLOAT
 co=0.0_FLOAT
energy=0.0_FLOAT
g=0.0_FLOAT

call JGB_criteria(thrR,Hlat,Itau_max)
CALL CRYALLOC(R0AB,MAX_ELEM,MAX_ELEM,'SRB_EGRAD','R0AB')
call setr0ab(max_elem,r0ab)

     ethr=0._FLOAT

do iat=1,n
   ITAU1=Itau_max(1)
   ITAU2=Itau_max(2)
   ITAU3=Itau_max(3)
   do i=-Itau1,Itau1
      do j=-Itau2,Itau2
         do k=-Itau3,Itau3
            do jat=1,n
               iselftest=0
               ! Test for equal atoms, remove selfinteraction
               if(iat.eq.jat) then
                  iselftest=iselftest+1
                  if(i.eq.0)iselftest=iselftest+1
                  if(j.eq.0)iselftest=iselftest+1
                  if(k.eq.0)iselftest=iselftest+1
               end if
               if(iselftest.eq.4)cycle
               dx=xyz(1,iat)-xyz(1,jat)+i*Hlat(1,1)+j*Hlat(1,2)+k*Hlat(1,3)
               dy=xyz(2,iat)-xyz(2,jat)+i*Hlat(2,1)+j*Hlat(2,2)+k*Hlat(2,3)
               dz=xyz(3,iat)-xyz(3,jat)+i*Hlat(3,1)+j*Hlat(3,2)+k*Hlat(3,3)
               r=sqrt(dx*dx+dy*dy+dz*dz)
               ! distance cutoff
               if(r.gt.thrR) cycle
               r0abij=r0ab(iz(iat),iz(jat))
               r0=rscal/r0ab(iz(iat),iz(jat))             
               fi=real(iz(iat))
               fj=real(iz(jat))
               ff=-(fi*fj)**0.5_FLOAT               
               ener_dum=qscal*ff*exp(-r0*r)
               !factor 1/2 from double counting? JGB: Yes.
               ener_dum=ener_dum*0.5_FLOAT
               if(abs(ener_dum).lt.Ethr) cycle
               co(iat)=co(iat)+ener_dum
               energy=energy+ener_dum               
               if(grad) then
                  rf=qscal/r
                  tmp=-ff*r0*exp(-r0*r)*rf
                  g(1,iat)=g(1,iat)+tmp*dx
                  g(2,iat)=g(2,iat)+tmp*dy
                  g(3,iat)=g(3,iat)+tmp*dz
                  stress(1,1)=stress(1,1)+dx*dx*tmp
                  stress(2,1)=stress(2,1)+dy*dx*tmp
                  stress(3,1)=stress(3,1)+dz*dx*tmp
                  
                  stress(1,2)=stress(1,2)+dx*dy*tmp
                  stress(2,2)=stress(2,2)+dy*dy*tmp
                  stress(3,2)=stress(3,2)+dz*dy*tmp
                  
                  stress(1,3)=stress(1,3)+dx*dz*tmp
                  stress(2,3)=stress(2,3)+dy*dz*tmp
                  stress(3,3)=stress(3,3)+dz*dz*tmp                  
                  
               endif               
            enddo !jat
         enddo !k
      enddo !j
   enddo !i
enddo !iat
if(echo)then
   write(IOUT,'(/2x,a5,2x,a5,4x,a15)') &
        '#','ON','SRB [kcal/mol]'
   do i=1,n
      write(IOUT,'(2x,2(i5,2x),F9.3)') i,iz(i),co(i)*AUTOKCAL
   enddo    
endif
call MINV3(Hlat,Hlat_1,Hdet)
 cellgrad=0.0_FLOAT
CALL MXMBN(stress,1,3,Hlat_1,3,1,cellgrad,1,3,3,3,3)
 cellgrad=0.5_FLOAT*cellgrad
RETURN
end subroutine srb_egrad2

