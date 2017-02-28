function fe=getdgdrhotilde_el(ex,ey,ep,D,ed_u,ed_rho,q)
%
%-------------------------------------------------------------
% PURPOSE
%   Computes the sensitivity of the compliance with respect to a
%   continues density field (rhotilde)
%
% INPUT:  ex = [x1 x2 x3 x4]  element coordinates
%         ey = [y1 y2 y3 y4]
%                             
%         ep =[ptype t ir]    element property 
%                               ptype: analysis type
%                               ir: integration rule
%                               t : thickness
%
%         D                   constitutive matrix
%
%     ed_u:             Element displacement vector = [u1, u2, ... u8];
%
%     ed_rho:           Element density vector = [rho1, rho2, rho3, rho4]
%
%     q:                Exponent in density field.
%
%
%
% OUTPUT: fe : Nodal sensitivites (4 x 1)
%-------------------------------------------------------------

% LAST MODIFIED: M Ristinmaa  1995-10-25
%                Eric Borgqvist  2014-03-04
% Copyright (c)  Division of Structural Mechanics and
%                Department of Solid Mechanics.
%                Lund Institute of Technology
%-------------------------------------------------------------
  ptype=ep(1); t=ep(2);  ir=ep(3);  ngp=ir*ir;

%--------- gauss points --------------------------------------
  if ir==1
    g1=0.0; w1=2.0;
    gp=[ g1 g1 ];  w=[ w1 w1 ];
  elseif ir==2
    g1=0.577350269189626; w1=1;
    gp(:,1)=[-g1; g1;-g1; g1];  gp(:,2)=[-g1;-g1; g1; g1];
    w(:,1)=[ w1; w1; w1; w1];   w(:,2)=[ w1; w1; w1; w1];
  elseif ir==3
    g1=0.774596699241483; g2=0.;
    w1=0.555555555555555; w2=0.888888888888888;
    gp(:,1)=[-g1;-g2; g1;-g1; g2; g1;-g1; g2; g1];
    gp(:,2)=[-g1;-g1;-g1; g2; g2; g2; g1; g1; g1];
    w(:,1)=[ w1; w2; w1; w1; w2; w1; w1; w2; w1];
    w(:,2)=[ w1; w1; w1; w2; w2; w2; w1; w1; w1];
  else
    disp('Used number of integration points not implemented');
    return
  end
  wp=w(:,1).*w(:,2);
  xsi=gp(:,1);  eta=gp(:,2);  r2=ngp*2;

%--------- shape functions -----------------------------------
  N(:,1)=(1-xsi).*(1-eta)/4;  N(:,2)=(1+xsi).*(1-eta)/4;
  N(:,3)=(1+xsi).*(1+eta)/4;  N(:,4)=(1-xsi).*(1+eta)/4;

  dNr(1:2:r2,1)=-(1-eta)/4;     dNr(1:2:r2,2)= (1-eta)/4;
  dNr(1:2:r2,3)= (1+eta)/4;     dNr(1:2:r2,4)=-(1+eta)/4;
  dNr(2:2:r2+1,1)=-(1-xsi)/4;   dNr(2:2:r2+1,2)=-(1+xsi)/4;
  dNr(2:2:r2+1,3)= (1+xsi)/4;   dNr(2:2:r2+1,4)= (1-xsi)/4;


  fe=zeros(4,1);
  JT=dNr*[ex;ey]';

%--------- plane stress --------------------------------------
if ptype==1
  
  colD=size(D,2);
  if colD>3
    Cm=inv(D);
    Dm=inv(Cm([1 2 4],[1 2 4]));
  else
    Dm=D;
  end
    
  for i=1:ngp
      indx=[ 2*i-1; 2*i ];
      detJ=det(JT(indx,:));
      if detJ<10*eps
        disp('Jacobideterminant equal or less than zero!')
      end
      JTinv=inv(JT(indx,:));
      dNx=JTinv*dNr(indx,:);

      B(1,1:2:8-1)=dNx(1,:);
      B(2,2:2:8)  =dNx(2,:);
      B(3,1:2:8-1)=dNx(2,:);
      B(3,2:2:8)  =dNx(1,:);

      N2(1,1:2:8-1)=N(i,:);
      N2(2,2:2:8)  =N(i,:);

      N1 = N(i,:);
      rhogp = N1*ed_rho';
      
      fe=fe+N1'*q*rhogp^(q-1)*(ed_u*B'*Dm*B*ed_u')*detJ*wp(i)*t;
  end
%--------- plane strain --------------------------------------
elseif ptype==2
  
  colD=size(D,2);
  if colD>3
    Dm=D([1 2 4],[1 2 4]);
  else
    Dm=D;
  end
    
  for i=1:ngp
      indx=[ 2*i-1; 2*i ];
      detJ=det(JT(indx,:));
      if detJ<10*eps
        disp('Jacobideterminant equal or less than zero!')
      end
      JTinv=inv(JT(indx,:));
      dNx=JTinv*dNr(indx,:);

      B(1,1:2:8-1)=dNx(1,:);
      B(2,2:2:8)  =dNx(2,:);
      B(3,1:2:8-1)=dNx(2,:);
      B(3,2:2:8)  =dNx(1,:);

      N2(1,1:2:8-1)=N(i,:);
      N2(2,2:2:8)  =N(i,:);

      N1 = N(i,:);
      rhogp = N1*ed_rho';
      
      fe=fe+N1'*q*rhogp^(q-1)*(ed_u*B'*Dm*B*ed_u')*detJ*wp(i)*t;
  end
  
else
   error('Error ! Check first argument, ptype=1 or 2 allowed')
   return
end
%--------------------------end--------------------------------
