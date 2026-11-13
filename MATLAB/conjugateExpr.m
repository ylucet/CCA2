% conjugate F*(s) = sup(s.x-f(x))  restricted to edge
% edge : edge
% f : f
% make appropriate substitutions to clean this

function conj = conjugateExpr(edge,f,x,y)
    %edge
    %f

    lambda = sym('lambda');
    infs = sym('inf');

    s1 = sym('s1');
    s2 = sym('s2');

    %edge
    %f
    dedge = sym.empty();
    dedge(1) = diff(edge,x);
    dedge(2) = diff(edge,y);

    df = sym.empty();
    df(1) = diff(f,x);
    df(2) = diff(f,y);
    
    eq1 = s1 -  df(1) - lambda*dedge(1);
    eq2 = s2 -  df(2) - lambda*dedge(2);
    eq3 = simplifyFraction(eq1*dedge(2)-eq2*dedge(1));
    eq3 = simplify(eq3);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    xyl = solve([eq1,eq2,edge],[x,y,lambda]);
    
    %size(xyl.x)
    if size(xyl.x,1) > 1
        xyl.x = xyl.x(2);
        xyl.y = xyl.y(2);
        xyl.lambda = xyl.lambda(2);
    end



     % xy = solve([eq3,edge],[x,y]);
    
    conj = infs;
    
    if isempty(xyl)
       return
    elseif isempty(xyl.x) | isempty(xyl.y) | isempty(xyl.lambda)
       return
    end
    
    %conj = s1*xy.x + s2*xy.y - subs(f,[x,y],[xy.x,xy.y]);
    conj = s1*xyl.x + s2*xyl.y - subs(f,[x,y],[xyl.x,xyl.y]);
    conj = simplifyFraction(conj);
    %[co, va] = coeffs(conj,[s1,s2])
    conj = subs(conj,[s1,s2],[x,y]);

    %[n,d] = numden(conj)
   
    %[c,t] = coeffs(n,[x,y])
    %[c,t] = coeffs(d,[x,y])
    
    
end