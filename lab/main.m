

parameter1=1;
parameter2=2;
lambda_min=.1;
lambda_max=100;


lambdastar = fzero(@(lambda) dphidlambda(lambda,parameter1,parameter2),[lambda_min lambda_max])
