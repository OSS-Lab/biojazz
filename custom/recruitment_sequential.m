function score = recruitment_sequential(t,A,X,Y,AXY) % A scoring function to encourage sequential recruitment
    A0 = A(1);
    X0 = X(1);
    Y0 = Y(1);

    
    Efficacy = AXY(end)/min(A0,X0,Y0);
    Aeff = AXY(end)/A0;
    Xeff = AXY(end)/X0;
    Yeff = AXY(end)/Y0;
    aggeff = 3*A(end)*X(end)*Y(end)/(A0+X0+Y0);
    %score = (Efficacy + aggeff)/2
    %score = Aeff
    %score = Xeff
    %score = Yeff
    %     score = aggeff
    k = 0.5 % importance scalar for intermediate AX
    score = (AXY(end)*AX(end)^k)^(1/(k+1))