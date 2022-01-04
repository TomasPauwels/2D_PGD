%flip the strain to be 
% [ r1 s1'
%   r2's2 
%   r1's1 + r2s1']

clear all; close all; clc;

%GEOMETRY
L=1;
H=1;
nelx=1;
nely=1;
[coorX,coorY]= meshgrid(0:(L/nelx):L,0:(H/nely):H); 
Nodes = [coorX(:),coorY(:)];

%VARIABLES
E = 1e6;
nu = 0.3;
error_iter=1;
error = 1;
iter=0;
mode_iter=0;
max_iter=25;
tau = 1e-3;
Aprt = 0;
TOL = 1e-3;

rho = 1;
g =10;
b = [0; -1];
%boundary conditions

BCx = [1 2];
BCy = [];

CCx=setxor([1:numel(coorX)],BCx);
CCy=setxor([1:numel(coorY)],BCy);

%Pre-allocation
Fhis=zeros(numel(coorX),1);
Ghis=zeros(numel(coorY),1);
F=zeros(numel(coorX),1);
G=zeros(numel(coorY),1);
qs = zeros(numel(coorY),1);
qr = rand(numel(coorX),1);
qs(BCx) = 0;
qr(BCy) = 0;

qr_old = qr;
qs_old = qs;

%Enrichment
%while error_iter >TOL && mode_iter<max_iter   
    mode_iter=mode_iter +1;
    
while abs(error) > tau && iter<max_iter
     iter = iter+1;
%% CALCULATE s ASSUMING r IS KNOWN
%Stiffness terms    
        ky = zeros(numel(coorY));
        gauss = [-1/sqrt(3) 1/sqrt(3)]; 
        Jx = L/nelx/2;                      %Jacobian x-direction
        Jy = H/nely/2;                      %Jacobian y-direction
for j = 1:size(gauss,2)

            R = [0                0.5*((1-gauss(j))*qr(1) + (1+gauss(j))*qr(3))   0                 0;
                 0                0                                               (qr(4)-qr(2))/2   0;
                 (qr(3)-qr(1))/2  0                                               0                 0.5*((1-gauss(j))*qr(2) + (1+gauss(j))*qr(4))];   
           
            C = E/(1-nu^2)*[1  nu 0;
                            nu 1  0;
                            0  0  (1-nu)/2];
            ky = ky + R'*C*R*Jy;
end


Kx = zeros(numel(coorX));
for j = 1:size(gauss,2)
    
   Bx = [ 0.5*(1-gauss(j)) 0                0.5*(1+gauss(j)) 0;
         -0.5              0                0.5              0
          0                0.5*(1-gauss(j)) 0                0.5*(1+gauss(j));
          0               -0.5              0                0.5];
   
    Kx = Kx + Bx'*ky*Bx *Jx;
end

%force terms
fy = zeros(2,1);
for j = 1:size(gauss,2)
      r = [0.5*((1-gauss(j))*qr(1) + (1+gauss(j))*qr(3));
           0.5*((1-gauss(j))*qr(2) + (1+gauss(j))*qr(4))];
     
  fy = fy + r.*b*Jx;
end
Fx = zeros(numel(coorX),1);
for j = 1:size(gauss,2)
    Nx = [0.5*(1-gauss(j))  0                   0.5*(1+gauss(j))  0;
          0                0.5*(1-gauss(j))    0                  0.5*(1+gauss(j))];

    Fx = Fx + Nx'*fy*Jx;
end
      
%residuals (negeren)
res=zeros(numel(coorX),1);
for j=1:size(gauss,2)
      R = [0                0.5*((1-gauss(j))*qr(1) + (1+gauss(j))*qr(3))   0                 0;
           0                0                                               (qr(4)-qr(2))/2   0;
           (qr(3)-qr(1))/2  0                                               0                 0.5*((1-gauss(j))*qr(2) + (1+gauss(j))*qr(4))];   
           
     C = E/(1-nu^2)*[1  nu 0;
                     nu 1  0;
                     0  0  (1-nu)/2];
                        
                        
      eta_u=[ (G(3)-G(1))/2                              *    (0.5*((1-gauss(j))*F(1) + (1+gauss(j))*F(3)));
              (0.5*((1-gauss(j))*G(2) + (1+gauss(j))*G(4)))*    (F(4)-F(2))/2;
              (0.5*((1-gauss(j))*G(1) + (1+gauss(j))*G(3)))*    (F(3)-F(1))/2                             + (G(4)-G(2))/2 * (0.5*((1-gauss(j))*F(2) + (1+gauss(j))*F(4)))];
              
      res = res + R'*C*eta_u * Jy;
end

Res=zeros(numel(coorY),1);
for j=1:size(gauss,2)
     Bx = [ 0.5*(1-gauss(j))  0                 0.5*(1+gauss(j))  0;
           -0.5               0                 0.5               0
            0                 0.5*(1-gauss(j))  0                 0.5*(1+gauss(j));
            0                -0.5               0                 0.5];
      
      Res= Res + Bx'*res *Jx;
end

                        
qs(CCx) = Kx(CCx,CCx)\(Fx(CCx)-Res(CCx))

%% Calculate r assuming s is known
kx = zeros(numel(coorY));
for j = 1:size(gauss,2)
    S = [0                0.5*((1-gauss(j))*qs(1) + (1+gauss(j))*qs(3))   0                 0;
         0                0                                               (qs(4)-qs(2))/2   0;
         (qs(3)-qs(1))/2  0                                               0                 0.5*((1-gauss(j))*qs(2) + (1+gauss(j))*qs(4))];   
    
    C = E/(1-nu^2)*[1  nu 0;
                    nu 1  0;
                    0  0  (1-nu)/2];
    kx = kx + S'*C*S *Jx;
end

Ky = zeros(numel(coorY));
for j = 1:size(gauss,2)
   By = [ 0.5*(1-gauss(j))            0                0.5*(1+gauss(j)) 0;
         -0.5                         0                0.5              0
          0                           0.5*(1-gauss(j)) 0                0.5*(1+gauss(j));
          0                          -0.5              0                0.5];
   Ky = Ky + By'*kx*By *Jy;
end

%force terms
fx = zeros(2,1);
for j = 1:size(gauss,2)
      s = [0.5*((1-gauss(j))*qs(1) + (1+gauss(j))*qs(3));
           0.5*((1-gauss(j))*qs(2) + (1+gauss(j))*qs(4))];    
      fx = fx + s.*b*Jx;
end
Fy = zeros(numel(coorX),1);
for j = 1:size(gauss,2)
    Ny = [0.5*(1-gauss(j))  0                   0.5*(1+gauss(j))  0;
          0                0.5*(1-gauss(j))    0                  0.5*(1+gauss(j))];

    Fy = Fy + Ny'*fx*Jy;
end
      

%residuals
res=zeros(numel(coorX),1);
for j=1:size(gauss,2)
      S = [0                0.5*((1-gauss(j))*qr(1) + (1+gauss(j))*qs(3))   0                 0;
           0                0                                               (qs(4)-qs(2))/2   0;
           (qs(3)-qs(1))/2  0                                               0                 0.5*((1-gauss(j))*qs(2) + (1+gauss(j))*qs(4))];   
          
      C = E/(1-nu^2)*[1  nu 0;
                      nu 1  0;
                      0  0  (1-nu)/2];
                        
                        
      eta_u=[ (G(3)-G(1))/2                              *    (0.5*((1-gauss(j))*F(1) + (1+gauss(j))*F(3)));
              (0.5*((1-gauss(j))*G(2) + (1+gauss(j))*G(4)))*    (F(4)-F(2))/2;
              (0.5*((1-gauss(j))*G(1) + (1+gauss(j))*G(3)))*    (F(3)-F(1))/2                             + (G(4)-G(2))/2 * (0.5*((1-gauss(j))*F(2) + (1+gauss(j))*F(4)))];
              
      res = res + R'*C*eta_u * Jy;
end

Res=zeros(numel(coorX),1);
for j=1:size(gauss,2)
    By= [ 0.5*(1-gauss(j))  0                0.5*(1+gauss(j)) 0;
         -0.5               0                0.5              0
          0                 0.5*(1-gauss(j)) 0                0.5*(1+gauss(j));
          0                -0.5              0                0.5];
     Res = Res + By'*res *Jx;
end

qr(CCy) = Ky(CCy,CCy)\(Fy(CCy)-Res(CCy))
qr = qr./norm(qr);

%% Error 
error = max(abs(sum(qr_old-qr)),abs(sum(qs_old-qs)))
qr_old = qr;
qs_old = qs;
end

