program my_combiner
  implicit none
  integer i
  INTEGER, PARAMETER :: N=10
  REAL, DIMENSION (N) ::  m0s,x0s,y0s,z0s,vx0s,vy0s,vz0s
  REAL, DIMENSION (N) ::  m_p,x0,y0,z0,vx0,vy0,vz0
  open(1,file='star-ms-xs-ys-zs-vxs-vys-vzs.txt')
  open(2,file='mp-x-y-z-vx-vy-vz-msun-pc.txt')
  open(3,file='test-N10-8o5kpc-v220.NBODY')
 do i=1,N
  read(1,*)m0s(i),x0s(i),y0s(i),z0s(i),vx0s(i),vy0s(i),vz0s(i)
  read(2,*)m_p(i),x0(i),y0(i),z0(i),vx0(i),vy0(i),vz0(i)
 end do
 
 do i=1,2*N
  if (i.le. N) then 
    write(3,*)m0s(i),x0s(i),y0s(i),z0s(i),vx0s(i),vy0s(i),vz0s(i)
  else 	
  write(3,*)m_p(i-N),x0(i-N)+x0s(i-N),y0(i-N)+y0s(i-N),z0(i-N)+z0s(i-N),vx0(i-N)+vx0s(i-N),vy0(i-N)+vy0s(i-N),vz0(i-N)+vz0s(i-N)
  end if
 end do
  close(1)
    close(2)
	  close(3)
end program