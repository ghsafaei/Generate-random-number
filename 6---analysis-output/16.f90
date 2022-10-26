
Program my_planetcircular
 IMPLICIT NONE
 REAL r0,theta,phi,vx0,vy0,vz0,a_min,a_max,pi,x0,y0,z0,G
 REAL M_star,AU,pc,M_Sun,m_p,M_J,mp_max,mp_min,M_Earth,ran3
 REAL vc,ran0,ran1,ran2,RANDOM
 INTEGER i,ISEED,idum ,seed
 INTEGER :: time,STIME
 INTEGER, DIMENSION (9) :: T
 REAL :: RANF,ecc,rp,ra,a,v_p,vpx0,vpy0,vpz0
 INTEGER N
PARAMETER(N=3)
REAL x(N),f(N)
LOGICAL check
common/constant/ x0,y0,z0,vx0,vy0,vz0,vc

!                            Random Number
integer :: values(1:8), kk
integer, dimension(:), allocatable :: seed2
call date_and_time(values=values)

 call random_seed(size=kk)
 allocate(seed2(1:kk))
 seed2(:) = values(8)
 call random_seed(put=seed2)
	!
	! Initial seed from the system time and forced to be odd
	!
	!STIME = time(%REF(0))
	CALL gmtime(TIME,T)
	ISEED = T(6)+70*(T(5)+12*(T(4)+31*(T(3)+23*(T(2)+59*T(1)))))
	IF (MOD(ISEED,2).EQ.0) ISEED = ISEED-1

	idum=iseed
	seed=iseed
	print*,seed,t,time,seed2!,	random_number(seed2)
   !initial value for 9 parameters
	m_p=-1
	x0=-1
	y0=-1
	z0=-1
	r0=-1
	vc=100.
	vx0=-1
	vy0=-1
	vz0=-1
    x(1)=100.
    x(2)=100.
    x(3)=100.
  open(1,file='mp-x-y-z-vx-vy-vz.txt')
  open(2,file='mp-x-y-z-vx-vy-vz-msun-pc.txt')
  Do i=1,2
 
    !Astronomical constant values
	 M_Sun=1.9891e30
	 M_star=1.*M_Sun
	 M_Earth=5.9736e24
	 pi=3.14159265358979323846264338327950
	 G=6.673e-11
	 AU=1.4959787066e11
	 pc=206264.806*AU
	 
     !a: Semi-major axis  (0.3 AU <a< 31 AU)
	 !r0= distance between star and planet at the start point of motion
	 a_min=0.3
	 a_max=31.
	 r0=((a_max-a_min)*RANDOM(Seed)+a_min)*AU
	 theta=2.*pi*RANDOM(Seed)
	 phi=pi*RANDOM(Seed)
	 r0=1.*AU
	 
	 !r0=sqrt(x0**2+y0**2+z0**2)
	 x0=r0*sin(theta)*cos(phi)
	 y0=r0*sin(theta)*sin(phi)
	 z0=r0*cos(theta)
  !The circular velocity at start point of motion
	 vc=sqrt (abs((G*M_star/r0)))
	 x0=x0/AU
     y0=y0/AU
     z0=z0/AU

  call newt(x,N,check)
  call funcv(N,x,f)
    if (check) then
      write(*,*) 'Convergence problems.'
    endif
     !write(*,'(1x,a5,t10,a1,t22,a1)') 'Index','x','f'

  vx0=x(1)
  vy0=x(2)
  vz0=x(3)

	!print*,G,M_star,r0,vc
	 !theta=2.*pi*RANDOM(Seed)
	 !phi=pi*RANDOM(Seed)
	 !vx0=vc*sin(theta)*cos(phi)
	 !vy0=vc*sin(theta)*sin(phi)
	 !vz0=vc*cos(theta)
	  !r=Sqrt(rx**2+ry**2+rz**2)
 !v=Sqrt(vx**2+vy**2+vz**2)
! vrx=vy*rz-vz*ry
! vry=vz*rx-vx*rz
 !vrz=vx*ry-vy*rx
! v*r = Sqrt(vrx**2+vry**2+vrz**2)
! vx*rx+vy*ry+vz*rz=0.
	!v0= sqrt(vx0**2+vy0**2+vz0**2)
	
	!eccentricity ecc:    0=<ecc<=1
	!ecc=RANDOM(Seed)
	!rp: peri center according to AU
	!a=r_p/(1.-ecc)
	!r_p=a*(1.-ecc)
	!a: Semi-major axis
	!a=r0
	!v_p=sqrt((G*M_star/a)*((1.+ecc)/(1.-ecc)))
	!theta=2.*pi*RANDOM(Seed)
	!phi=pi*RANDOM(Seed)
	!vpx0=v_p*sin(theta)*cos(phi)
	!vpy0=v_p*sin(theta)*sin(phi)
	!vpz0=v_p*cos(theta)
	!theta=pi*RANDOM(Seed)
	!r=a*(1.-ecc**2.)/(1.+ecc*cos(theta))
	!x=a*cosh(r)*cos(theta)
	!y=a*sinh(r)*sin(theta)
	!z=z
	
	 M_J=317.8*M_Earth
	 mp_max=8.*M_J
	 mp_min=0.0001*M_J
	 m_p=(mp_max-mp_min)*RANDOM(Seed)+mp_min
	 m_p=M_Earth
	write(1,*)m_p,x0,y0,z0,vx0,vy0,vz0 ,x0*vx0+y0*vy0+z0*vz0,Sqrt(vx0**2+vy0**2+vz0**2)
    write(2,*)m_p/M_Sun,x0/(206264.806),y0/(206264.806),z0/(206264.806),vx0/1000.,&
	&vy0/1000.,vz0/1000.,Sqrt(vx0**2+vy0**2+vz0**2)/1000.,vx0*x0+vy0*y0+vz0*z0,&
	&acos((vx0*x0+vy0*y0+vz0*z0)/(Sqrt(vx0**2+vy0**2+vz0**2)*r0))
	!write(1,*)m_p,x0,y0,z0,vpx0,vpy0,vpz0
	!write(1,*)'m_p=',m_p,'   ','x0=',x0,'   ','y0=',y0,'   ','z0=',z0,'   ','vx0=',vx0,'   ','vy0=',vy0,'   ','vz0=',vz0,'   ','vc=',vc,'   ','r0=',sqrt(x0**2+y0**2+z0**2),'   ','v0=',sqrt(vx0**2+vy0**2+vz0**2)
  !print*,ran3(idum),'end'
  !RANF()	
  END Do	 
	close(1) 
	close(2)
End program
!!!!!!!!!!
!*****************************************************
SUBROUTINE funcv(n,x,f)
INTEGER n
REAL x(n),f(n),vx0,vy0,vz0,vrx,vry,vrz,v0,r0,vc,x0,y0,z0
common/constant/ x0,y0,z0,vx0,vy0,vz0,vc
!vc=29870.
!print*,'***funcv'
!print*,'***funcv',x0,y0,z0,Sqrt(x0**2+y0**2+z0**2)

r0=Sqrt(x0**2+y0**2+z0**2)
 !v=Sqrt(vx**2+vy**2+vz**2)
 !x(1)=vx
 !x(2)=vy
 !x(3)=vz
 v0=Sqrt(x(1)**2+x(2)**2+x(3)**2)
!vrx=vy*rz-vz*ry
vrx=x(2)*z0-x(3)*y0
!vry=vz*rx-vx*rz
vry=x(3)*x0-x(1)*z0
!vrz=vx*ry-vy*rx
vrz=x(1)*y0-x(2)*x0

! v*r = Sqrt(vrx**2+vry**2+vrz**2)
! vx*rx+vy*ry+vz*rz=0.
	!v0= sqrt(vx0**2+vy0**2+vz0**2)
!f(1)=Acos((x(1)*x0+x(2)*y0+x(3)*z0)/(v0*r0))-(3.14/2.)	
f(1)=x(1)*x0+x(2)*y0+x(3)*z0
f(2)=v0*r0 -Sqrt(vrx**2.+vry**2.+vrz**2.)
f(3)=v0-vc
!print*,'***funcv',v0*r0,Sqrt(vrx**2.+vry**2.+vrz**2.)
!print*,'***funcv222',x(1)*x0+x(2)*y0+x(3)*z0
return
END
	





!******************************newt.for
	
SUBROUTINE newt(x,n,check)
INTEGER n,nn,NP,MAXITS
LOGICAL check
REAL x(n),fvec,TOLF,TOLMIN,TOLX,STPMX
PARAMETER (NP=40,MAXITS=200,TOLF=1.e-4,TOLMIN=1.e-6,TOLX=1.e-7,STPMX=100.)
COMMON /newtv/ fvec(NP),nn
SAVE /newtv/
!CU USES fdjac,fmin,lnsrch,lubksb,ludcmp
INTEGER i,its,j,indx(NP)
REAL d,den,f,fold,stpmax,sum,temp,test,fjac(NP,NP),g(NP),p(NP),xold(NP),fmin
EXTERNAL fmin
nn=n
f=fmin(x)
test=0.
do 11 i=1,n
if(abs(fvec(i)).gt.test)test=abs(fvec(i))
	11 continue
if(test.lt..01*TOLF)then
check=.false.
return
endif
sum=0.
do 12 i=1,n
sum=sum+x(i)**2
12 continue
stpmax=STPMX*max(sqrt(sum),float(n))
do 21 its=1,MAXITS
call fdjac(n,x,fvec,NP,fjac)
do 14 i=1,n
sum=0.
do 13 j=1,n
sum=sum+fjac(j,i)*fvec(j)
13 continue
g(i)=sum
14 continue
do 15 i=1,n
xold(i)=x(i)
15 continue
fold=f
do 16 i=1,n
p(i)=-fvec(i)
16 continue
call ludcmp(fjac,n,NP,indx,d)
call lubksb(fjac,n,NP,indx,p)
call lnsrch(n,xold,fold,g,p,x,f,stpmax,check,fmin)
test=0.
	do 17 i=1,n
if(abs(fvec(i)).gt.test)test=abs(fvec(i))
17 continue
if(test.lt.TOLF)then
check=.false.
return
endif
if(check)then
test=0.
den=max(f,.5*n)
do 18 i=1,n
temp=abs(g(i))*max(abs(x(i)),1.)/den
if(temp.gt.test)test=temp
18 continue
if(test.lt.TOLMIN)then
check=.true.
else
check=.false.
endif
return
endif
test=0.
do 19 i=1,n
temp=(abs(x(i)-xold(i)))/max(abs(x(i)),1.)
if(temp.gt.test)test=temp
19 continue
if(test.lt.TOLX)return
	21 continue
pause 'MAXITS exceeded in newt'
END

!******************************************ludcmp.for
	
SUBROUTINE ludcmp(a,n,np,indx,d)
INTEGER n,np,indx(n),NMAX
REAL d,a(np,np),TINY
PARAMETER (NMAX=500,TINY=1.0e-20)
INTEGER i,imax,j,k
REAL aamax,dum,sum,vv(NMAX)
d=1.
do 12 i=1,n
aamax=0.
do 11 j=1,n
if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
11 continue
if (aamax.eq.0.) pause 'singular matrix in ludcmp'
vv(i)=1./aamax
12 continue
do 19 j=1,n
do 14 i=1,j-1
sum=a(i,j)
do 13 k=1,i-1
sum=sum-a(i,k)*a(k,j)
13 continue
a(i,j)=sum
14 continue
aamax=0.
do 16 i=j,n
	sum=a(i,j)
do 15 k=1,j-1
sum=sum-a(i,k)*a(k,j)
15 continue
a(i,j)=sum
dum=vv(i)*abs(sum)
if (dum.ge.aamax) then
imax=i
aamax=dum
endif
16 continue
if (j.ne.imax)then
do 17 k=1,n
dum=a(imax,k)
a(imax,k)=a(j,k)
a(j,k)=dum
17 continue
d=-d
vv(imax)=vv(j)
endif
indx(j)=imax
if(a(j,j).eq.0.)a(j,j)=TINY
if(j.ne.n)then
dum=1./a(j,j)
	do 18 i=j+1,n
a(i,j)=a(i,j)*dum
18 continue
endif
19 continue
return
END

!**********************************************lubksb.for
	
SUBROUTINE lubksb(a,n,np,indx,b)
INTEGER n,np,indx(n)
REAL a(np,np),b(n)
INTEGER i,ii,j,ll
REAL sum
ii=0
do 12 i=1,n
ll=indx(i)
sum=b(ll)
b(ll)=b(i)
if (ii.ne.0)then
do 11 j=ii,i-1
sum=sum-a(i,j)*b(j)
11 continue
else if (sum.ne.0.) then
ii=i
endif
b(i)=sum
12 continue
do 14 i=n,1,-1
sum=b(i)
do 13 j=i+1,n
sum=sum-a(i,j)*b(j)
13 continue
b(i)=sum/a(i,i)
14 continue
return
END
	
!@*************************************lnsrch.for
	
SUBROUTINE lnsrch(n,xold,fold,g,p,x,f,stpmax,check,func)
INTEGER n
LOGICAL check
REAL f,fold,stpmax,g(n),p(n),x(n),xold(n),func,ALF,TOLX
PARAMETER (ALF=1.e-4,TOLX=1.e-7)
EXTERNAL func
!CU USES func
INTEGER i
REAL a,alam,alam2,alamin,b,disc,f2,fold2,rhs1,rhs2,slope,sum,temp,test,tmplam
check=.false.
sum=0.
do 11 i=1,n
sum=sum+p(i)*p(i)
11 continue
sum=sqrt(sum)
if(sum.gt.stpmax)then
do 12 i=1,n
p(i)=p(i)*stpmax/sum
12 continue
endif
slope=0.
do 13 i=1,n
slope=slope+g(i)*p(i)
	13 continue
test=0.
do 14 i=1,n
temp=abs(p(i))/max(abs(xold(i)),1.)
if(temp.gt.test)test=temp
14 continue
alamin=TOLX/test
alam=1.
1 continue
do 15 i=1,n
x(i)=xold(i)+alam*p(i)
15 continue
f=func(x)
if(alam.lt.alamin)then
do 16 i=1,n
x(i)=xold(i)
16 continue
check=.true.
return
else if(f.le.fold+ALF*alam*slope)then
return
else
if(alam.eq.1.)then
tmplam=-slope/(2.*(f-fold-slope))
else
rhs1=f-fold-alam*slope
	rhs2=f2-fold2-alam2*slope
a=(rhs1/alam**2-rhs2/alam2**2)/(alam-alam2)
b=(-alam2*rhs1/alam**2+alam*rhs2/alam2**2)/(alam-alam2)
if(a.eq.0.)then
tmplam=-slope/(2.*b)
else
disc=b*b-3.*a*slope
if(disc.lt.0.) pause 'roundoff problem in lnsrch'
tmplam=(-b+sqrt(disc))/(3.*a)
endif
if(tmplam.gt..5*alam)tmplam=.5*alam
endif
endif
alam2=alam
	f2=f
fold2=fold
alam=max(tmplam,.1*alam)
goto 1
END
!***********************************fmin.for
	
FUNCTION fmin(x)
INTEGER n,NP
REAL fmin,x(*),fvec
PARAMETER (NP=40)
COMMON /newtv/ fvec(NP),n
SAVE /newtv/
!CU USES funcv
INTEGER i
REAL sum
call funcv(n,x,fvec)
sum=0.
do 11 i=1,n
sum=sum+fvec(i)**2
11 continue
fmin=0.5*sum
return
END

!************************************fdjac.for
	
SUBROUTINE fdjac(n,x,fvec,np,df)
INTEGER n,np,NMAX
REAL df(np,np),fvec(n),x(n),EPS
PARAMETER (NMAX=40,EPS=1.e-4)
!CU USES funcv
INTEGER i,j
REAL h,temp,f(NMAX)
do 12 j=1,n
temp=x(j)
h=EPS*abs(temp)
if(h.eq.0.)h=EPS
x(j)=temp+h
h=x(j)-temp
call funcv(n,x,f)
x(j)=temp
do 11 i=1,n
df(i,j)=(f(i)-fvec(i))/h
11 continue
12 continue
return
END
!///////////////////////////////ran3.for
FUNCTION ran3(idum)
INTEGER idum
INTEGER MBIG,MSEED,MZ
!C REAL MBIG,MSEED,MZ
REAL ran3,FAC
PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1./MBIG)
!C PARAMETER (MBIG=4000000.,MSEED=1618033.,MZ=0.,FAC=1./MBIG)
INTEGER i,iff,ii,inext,inextp,k
INTEGER mj,mk,ma(55)
!C REAL mj,mk,ma(55)
SAVE iff,inext,inextp,ma
DATA iff /0/
if(idum.lt.0.or.iff.eq.0)then
iff=1
mj=MSEED-iabs(idum)
mj=mod(mj,MBIG)
ma(55)=mj
	mk=1
do 11 i=1,54
ii=mod(21*i,55)
ma(ii)=mk
mk=mj-mk
if(mk.lt.MZ)mk=mk+MBIG
mj=ma(ii)
11 continue
do 13 k=1,4
do 12 i=1,55
ma(i)=ma(i)-ma(1+mod(i+30,55))
if(ma(i).lt.MZ)ma(i)=ma(i)+MBIG
12 continue
13 continue
inext=0
inextp=31
idum=1
endif
inext=inext+1
if(inext.eq.56)inext=1
inextp=inextp+1
if(inextp.eq.56)inextp=1
mj=ma(inext)-ma(inextp)
	if(mj.lt.MZ)mj=mj+MBIG
ma(inext)=mj
ran3=mj*FAC
return
END

!/////////////////////////////ran0.for	
FUNCTION ran0(idum)
INTEGER idum,IA,IM,IQ,IR,MASK
REAL ran0,AM
PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,MASK=123459876)
INTEGER k
idum=ieor(idum,MASK)
k=idum/IQ
idum=IA*(idum-k*IQ)-IR*k
if (idum.lt.0) idum=idum+IM
ran0=AM*idum
idum=ieor(idum,MASK)
return
END

!////////////////////////ran1.for	
FUNCTION ran1(idum)
INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
REAL ran1,AM,EPS,RNMX
PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,&
&NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
INTEGER j,k,iv(NTAB),iy
SAVE iv,iy
DATA iv /NTAB*0/, iy /0/
if (idum.le.0.or.iy.eq.0) then
idum=max(-idum,1)
do 11 j=NTAB+8,1,-1
k=idum/IQ
idum=IA*(idum-k*IQ)-IR*k
if (idum.lt.0) idum=idum+IM
if (j.le.NTAB) iv(j)=idum
11 continue
iy=iv(1)
endif
k=idum/IQ
idum=IA*(idum-k*IQ)-IR*k
if (idum.lt.0) idum=idum+IM
	j=1+iy/NDIV
iy=iv(j)
iv(j)=idum
ran1=min(AM*iy,RNMX)
return
END

!/////////////////////////// ran2.for	
FUNCTION ran2(idum)
INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
REAL ran2,AM,EPS,RNMX
PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1, &
&IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791, &
&NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
INTEGER idum2,j,k,iv(NTAB),iy
SAVE iv,iy,idum2
DATA idum2/123456789/, iv/NTAB*0/, iy/0/
if (idum.le.0) then
idum=max(-idum,1)
idum2=idum
do 11 j=NTAB+8,1,-1
	k=idum/IQ1
idum=IA1*(idum-k*IQ1)-k*IR1
if (idum.lt.0) idum=idum+IM1
if (j.le.NTAB) iv(j)=idum
11 continue
iy=iv(1)
endif
k=idum/IQ1
idum=IA1*(idum-k*IQ1)-k*IR1
if (idum.lt.0) idum=idum+IM1
k=idum2/IQ2
idum2=IA2*(idum2-k*IQ2)-k*IR2
if (idum2.lt.0) idum2=idum2+IM2
j=1+iy/NDIV
iy=iv(j)-idum2
iv(j)=idum
if(iy.lt.1)iy=iy+IMM1
ran2=min(AM*iy,RNMX)
return
END

!#############
     Function RANDOM(Seed)
     Implicit None
     Real :: Random
     Integer :: Seed
     Integer :: OldSeed = 0
     Integer, Parameter :: C1 = 19423
     Integer, Parameter :: C2 = 811
     Save OldSeed

     If (OldSeed .EQ. 0) OldSeed = Seed
     OldSeed = Mod(C1 * OldSeed, C2)
     RANDOM = 1.0 * OldSeed / C2

     End Function RANDOM
	 
!!!!!!!!!!!!!	 
FUNCTION RANF() RESULT (R)
!
! Uniform random number generator x(n+1) = a*x(n) mod c with
! a=7**5 and c = 2**(31)-1.
!
!USE CSEED
IMPLICIT NONE
INTEGER :: IH,IL,IT,IA,IC,IQ,IR,ISEED
DATA IA/16807/,IC/2147483647/,IQ/127773/,IR/2836/
REAL :: R
!
IH = ISEED/IQ
IL = MOD(ISEED,IQ)
IT = IA*IL-IR*IH
IF(IT.GT.0) THEN
ISEED = IT
ELSE
ISEED = IC+IT
END IF
R = ISEED/FLOAT(IC)
END FUNCTION RANF
