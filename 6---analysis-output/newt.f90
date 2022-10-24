!****************************xnewt.for
	
PROGRAM xnewt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!C driver for routine newt
INTEGER N
PARAMETER(N=3)
INTEGER i
REAL x(N),f(N)
LOGICAL check

x(1)=-1.
x(2)=1.
x(3)=1.

call newt(x,N,check)
call funcv(N,x,f)

if (check) then
write(*,*) 'Convergence problems.'
endif

write(*,'(1x,a5,t10,a1,t22,a1)') 'Index','x','f'

do  i=1,N
write(*,'(1x,i2,2x,2f12.6)') i,x(i),f(i)
end do

END
!*****************************************************
SUBROUTINE funcv(n,x,f)
INTEGER n
REAL x(n),f(n),rx,ry,rz,vx,vy,vz,vrx,vry,vrz,v,r

rx=2.
ry=3.
rz=1.
r=Sqrt(rx**2+ry**2+rz**2)
 !v=Sqrt(vx**2+vy**2+vz**2)
 !x(1)=vx
 !x(2)=vy
 !x(3)=vz
 v=Sqrt(x(1)**2+x(2)**2+x(3)**2)
!vrx=vy*rz-vz*ry
vrx=x(2)*rz-x(3)*ry
!vry=vz*rx-vx*rz
vry=x(3)*rx-x(1)*rz
!vrz=vx*ry-vy*rx
vrz=x(1)*ry-x(2)*rx

! v*r = Sqrt(vrx**2+vry**2+vrz**2)
! vx*rx+vy*ry+vz*rz=0.
	!v0= sqrt(vx0**2+vy0**2+vz0**2)
f(1)=x(1)*rx+x(2)*ry+x(3)*rz
f(2)=v*r -Sqrt(vrx**2.+vry**2.+vrz**2.)
f(3)=v-29.8
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

