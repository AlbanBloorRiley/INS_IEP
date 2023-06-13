function [x,NIter,flag,xs,Fprev] = Deflated_GaussNewton_Type1_Method(Evaluate,xint,constants,tol,MaxIter,Y,p,linesearch_type,c1,alphaint,tau,minalpha,theta,f,f_0,J_fformula)

%Initial Values
xs = xint;  NIter = 0;  hprev = inf(length(xint),1); xprev = xint;
Fprev=Evaluate(xint,constants);

while true
    NIter = NIter +1;   alpha = alphaint;
    [M,gradM] = deflation(Y,xprev,p);%Calculate deflation operators
    [~,Rx,Jx] = Evaluate(xprev,constants);%Calculate residual, Jacobian of R

    h =  GN_Step(Jx,Rx,M,gradM);     %Calculate Gauss-Newton Step
    x = xprev- alpha*h;
    linesearch_1     %Apply linesearch

    xs = [xs,x];    %Saves x values for plotting

    Fprev = Evaluate(x,constants)      %Calculate objective function


    [tf,flag] = IsMin(Fprev,h,hprev,Jx'*Rx,NIter,tol,MaxIter); %Stopping criteria

    if tf
        return
    end
    xprev = x;
    hprev = h;
end
end

%Gauss--Newton Step
function [h,J] = GN_Step(J,R,M,J_M)
h = lsqminnorm(J'*J,J'*R);
h = h/(1+(1/M)*J_M*h);
end
