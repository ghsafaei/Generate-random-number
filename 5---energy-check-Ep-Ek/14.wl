(* ::Package:: *)

M=1.9891 10^30;
G=6.673 10^-11;
mearth=5.9736 10^24;
r=1.4959787066 10^11;
v=29780.;
Print["Poteitial energy E_p="]
Ep=-((G M mearth)/r)
Print["Kinetic energy E_k="]
Ek=0.5 mearth v^2
Print["U = 2 k"]
Ep/Ek
Print["Centrifugal force="]
(mearth v^2)/r
Print["Gravitational force="]
(G M mearth)/r^2

rx=2;
ry=3;
rz=1;

r=\[Sqrt](rx^2+ry^2+rz^2);
vx*rx+vy*ry+vz*rz=0;
v=\[Sqrt](vx^2+vy^2+vz^2);

vrx=vy*rz-vz*ry;
vry=vz*rx-vx*rz;
vrz=vx*ry-vy*rx;





