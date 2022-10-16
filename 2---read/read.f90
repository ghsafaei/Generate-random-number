Program my_reader
  Implicit none
  integer i
  INTEGER, PARAMETER :: N=10
  REAL, DIMENSION (N) ::  m0s,x0s,y0s,z0s,vx0s,vy0s,vz0s
  open(1,file='test-N10-8o5kpc-v220.PAR.NBODY')
  open(2,file='star-ms-xs-ys-zs-vxs-vys-vzs.txt')
 do i=1,N
  read(1,*)m0s(i),x0s(i),y0s(i),z0s(i),vx0s(i),vy0s(i),vz0s(i)
  print*,m0s(i),x0s(i),y0s(i),z0s(i),vx0s(i),vy0s(i),vz0s(i)
  write(2,*)m0s(i),x0s(i),y0s(i),z0s(i),vx0s(i),vy0s(i),vz0s(i)
 end do
 close(1)
 close(2)
End program
