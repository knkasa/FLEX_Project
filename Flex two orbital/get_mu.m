function out = get_mu(S,P,WN,ek,muLim,xifill,norb,beta,useSymmetry)
%GET_MU Find the chemical potential with bisection method.
%   The bisection method is guaranteed to converge if DOS(w) >= 0.

%erreps = 1e-4;
erreps = 1e-8;

fill = 2*norb*fermi(xifill,beta);       

%set the initial values of the chemical potential.
muL = muLim(1);
muR = muLim(2);
muM = (muR+muL)/2;

fillL = get_filling(S,P,WN,ek,muL,xifill,norb,beta,useSymmetry);
fillR = get_filling(S,P,WN,ek,muR,xifill,norb,beta,useSymmetry);

%In some rare cases, the self-energy X is so large that fill is not in
%[fillL fillR], so we expand the interval [muL muR] below.
iter = 0;
Wband = muR - muL;
while ( (fill-fillL)*(fill-fillR)>0 )
    iter = iter + 1;
    if (fill > fillR)
        muR = muR + Wband*0.5;
        fillR = get_filling(S,P,WN,ek,muR,xifill,norb,beta,useSymmetry);
    elseif (fill < fillL)
        muL = muL - Wband*0.5;
        fillL = get_filling(S,P,WN,ek,muL,xifill,norb,beta,useSymmetry);
    end
    if (iter > 4)
        error('Self-energy too large, cannot find the chemical potential!')
    end
end
muLim = [muL muR];

%find mu using Matlab fzero method (super-linear convergence) -------------
%fzero uses the algorithm that was originated by T. Dekker, which uses a
%combination of bisection, secant, and inverse quadratic interpolation
%methods. See Press et al., Numerical Recipes, Chapter 9.3 Van
%Wijngaarden-Dekker-Brent Method.
options = optimset;
options.TolX = erreps;
out = fzero(@(xmu) get_filling(S,P,WN,ek,xmu,xifill,norb,beta,useSymmetry)-fill,muLim,options);
