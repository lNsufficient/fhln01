function Me = flw2i4m(ex,ey,x)
% Me=flw2i4m(ex,ey,x)
%-------------------------------------------------------------
% PURPOSE
%  Compute the quantity: Ce=x*int(N^T*N)dA with full integration, i.e. 3x3
%  gauss points
%
% INPUT:  ex,ey;       Element coordinates
%	
%	  x
%
% OUTPUT: Theta :      Matix 4 x 4
%-------------------------------------------------------------

w1 = [0.55555555555555555555,0.8888888888888888,0.55555555555555555555,0.55555555555555555555,0.8888888888888888,0.55555555555555555555,0.55555555555555555555,0.8888888888888888,0.55555555555555555555];
w2 = [0.55555555555555555555,0.55555555555555555555,0.55555555555555555555,0.8888888888888888,0.8888888888888888,0.8888888888888888,0.55555555555555555555,0.55555555555555555555,0.55555555555555555555];
gp1 = -0.774596669241483;
gp2 = 0;
gp3 = -gp1;

gp_m = [gp1 gp1;
    gp2 gp1;
    gp3 gp1;
    gp1 gp2;
    gp2 gp2;
    gp3 gp2;
    gp1 gp3;
    gp2 gp3;
    gp3 gp3;];

Me = zeros(4);
for gp=1:9
     
     xsi = gp_m(gp,1);
     eta = gp_m(gp,2);

     dN1dxsi = 0.25*(eta-1);
     dN1deta = 0.25*(xsi-1);
     dN2dxsi = -0.25*(eta-1);
     dN2deta = -0.25*(xsi+1);
     dN3dxsi = 0.25*(eta+1);
     dN3deta = 0.25*(xsi+1);
     dN4dxsi = -0.25*(eta+1);
     dN4deta = -0.25*(xsi-1);
     
     dNdxsi = [dN1dxsi, dN2dxsi, dN3dxsi, dN4dxsi];
     dNdeta = [dN1deta, dN2deta, dN3deta, dN4deta];
     ax = [ex(1); ex(2); ex(3); ex(4)];
     ay = [ey(1); ey(2); ey(3); ey(4)];
    
    dxdxsi = dNdxsi*ax; dxdeta= dNdeta*ax;
    dydxsi = dNdxsi*ay; dydeta = dNdeta*ay;

    J = [dxdxsi, dxdeta; dydxsi, dydeta];
    detJ = det(J); %ska vara lika med Areaan av elementen*4
     
   N1 = 1/4*(xsi-1)*(eta-1);
   N2 = -1/4*(xsi+1)*(eta-1);
   N3 = 1/4*(xsi+1)*(eta+1);
   N4 = -1/4*(xsi-1)*(eta+1);
    
   N = [N1, N2,N3,N4];
     
   Me = Me + N'*N*detJ*x*w1(gp)*w2(gp);
    
end
