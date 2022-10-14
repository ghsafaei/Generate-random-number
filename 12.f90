
Program my_planetcircular
 IMPLICIT NONE
 REAL r0,theta,phi,vx0,vy0,vz0,a_min,a_max,pi,x0,y0,z0,G
 REAL M_star,AU,pc,M_Sun,m_p,M_J,mp_max,mp_min,M_Earth,ran3
 REAL vc,ran0,ran1,ran2,RANDOM
 INTEGER i,ISEED,idum,seed
 INTEGER :: time,STIME
 INTEGER, DIMENSION (9) :: T
 REAL :: RANF,ecc,rp,ra,a,v_p,vpx0,vpy0,vpz0
	!
	! Initial seed from the system time and forced to be odd
	!
	STIME = time(%REF(0))
	CALL gmtime(STIME,T)
	ISEED = T(6)+70*(T(5)+12*(T(4)+31*(T(3)+23*(T(2)+59*T(1)))))
	IF (MOD(ISEED,2).EQ.0) ISEED = ISEED-1

	idum=iseed
	seed=iseed
   !initial value for 9 parameters
	m_p=-1
	x0=-1
	y0=-1
	z0=-1
	r0=-1
	vc=-1
	vx0=-1
	vy0=-1
	vz0=-1

  open(1,file='mp-x-y-z-vx-vy-vz.txt')
 
  Do i=1,1
 
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
	!print*,G,M_star,r0,vc
	 theta=2.*pi*RANDOM(Seed)
	 phi=pi*RANDOM(Seed)
	 vx0=vc*sin(theta)*cos(phi)
	 vy0=vc*sin(theta)*sin(phi)
	 vz0=vc*cos(theta)
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
	 !m_p=M_Earth
	write(1,*)m_p,x0,y0,z0,vx0,vy0,vz0,x0*vx0+y0*vy0+z0*vz0
	!write(1,*)m_p,x0,y0,z0,vpx0,vpy0,vpz0
	!write(1,*)'m_p=',m_p,'   ','x0=',x0,'   ','y0=',y0,'   ','z0=',z0,'   ','vx0=',vx0,'   ','vy0=',vy0,'   ','vz0=',vz0,'   ','vc=',vc,'   ','r0=',sqrt(x0**2+y0**2+z0**2),'   ','v0=',sqrt(vx0**2+vy0**2+vz0**2)
  !print*,ran3(idum),'end'
  !RANF()	
  END Do	 
	close(1) 
End program
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