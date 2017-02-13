function xstar=xstar(lambda,parameter1,parameter2)

tmp=lambda*parameter1*parameter2;

if tmp<0
   xstar=0.5;
else
   xstar=1;
end
  


