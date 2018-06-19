function out = get_filling(S,P,WN,ek,mu,xifill,norb,beta,useSymmetry)

%S and P should be 4-dimensional array and ek 3-dimensional array. G only
%includes all the diagonal matrix elements in orbital space. 

[nk1 nk2 nw nel] = size(S);
if mod(nw,2) ~= 0 || nw ~= numel(WN)
    error('Wrong number of frequency grid!')
end
Nk(2) = nk2/2;
Nk(1) = nk1/2;
nkpt = nk1*nk2;
numwi = nw/2;
G = solve_dyson(S,P,WN,ek,mu,norb,Nk,numwi,useSymmetry,0);

dim = size(S);
if dim(4)==2
  G = G(:,:,:,[1 1]);     %(1,1) is equivalent to (2,2)
elseif dim(4)==3
  G = G(:,:,:,[1 3]);     %(1,1) ~= (2,2)
end


%  my code for Nelm=3
G = real(G);
fill = 2*norb*fermi(xifill,beta);  
for nn = 1:length(WN)
     wn = WN(nn);
     fill = fill + 2/nkpt/beta*( sum(sum(G(:,:,nn,1))) + sum(sum(G(:,:,nn,2)))  ) ...
         + (2*norb/beta)*xifill/(wn^2 + xifill^2);
end



%{
%  Yan's code for Nelm=2
fill = 2*norb*fermi(xifill,beta);       %for two spins
G = reshape(permute(G,[1 2 4 3]),nkpt*norb,[]);
G = real(G(:,1:numwi));
for nn = 1:numwi
    wn = WN(nn);
    fill = fill + 4*sum(G(:,nn))/(nkpt*beta) ...
        + (4*norb/beta)*xifill/(wn^2 + xifill^2);
end
%}

        
out = fill;
