function [S P] = solve_dysonbwd(G,F,WN,ek,mu,norb,Nk,Nc,useSymmetry,varargin)
%SOLVE_DYSONBWD Solve Dyson's equation backward to find normal (S) and
%anomalous (P) self-energy from normal (G) and anomalous (F) Green's
%function.


useFreqSymmetry = 1;

% useFreqSymmetry: =1 use frequency symmetry G(k,-i*w_n) = conj[G(k,i*w_n)]
% and F(k,-i*w_n) = conj[F(k,i*w_n)].

if norb~=2
    error('Currently, only two-orbital case is supported!')
end

if isempty(varargin)
    solveCase = 1;  %default: calculate both G and F
else
    solveCase = 0;  %calculate only G
end


siz_G = size(G);
siz_e = size(ek);
if ~isequal(siz_G(1:3),[2*Nk(1) 2*Nk(2) 2*Nc])
    error(' Input array G has a wrong size.')
end
if ~isequal(siz_e(1:2),[2*Nk(1) 2*Nk(2)])
    error(' Input array ek has a wrong size.')
end
if size(WN,1)~=1
    WN = reshape(WN,1,[]);
end
nel = siz_G(4);     %number of elements of the matrix in orbital space stored
if ~isequal(siz_e(3),nel)
    error(' Input array ek or G has a wrong size.')
end

if useSymmetry
    G = G(1:Nk(1)+1,1:Nk(2)+1,:,:);
    F = F(1:Nk(1)+1,1:Nk(2)+1,:,:);
    ek = ek(1:Nk(1)+1,1:Nk(2)+1,:);
    nk1 = Nk(1)+1;
    nk2 = Nk(2)+1;
else
    nk1 = 2*Nk(1);
    nk2 = 2*Nk(2);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if useFreqSymmetry
    G = G(:,:,1:Nc,:);
    F = F(:,:,1:Nc,:);
    WN = WN(1:Nc);
    nwpt = Nc;
else
    nwpt = 2*Nc;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nkpt = nk1*nk2;
npt = nkpt*nwpt;
norbsq = norb*norb;

dim = size(G);
if dim(4)==2
  idxmap = [1 2 2 1];       %(1,1) is equivalent to (2,2)
elseif dim(4)==3
  idxmap = [1 2 2 3];       %(1,1) is not equivalent to (2,2)
end

G = G(:,:,:,idxmap);
F = F(:,:,:,idxmap);
G = reshape(G,npt,norbsq);   %now G becomes [nk1*nk2*nw norbsq]
F = reshape(F,npt,norbsq);
%F is Hermitian
F(:,3) = conj(F(:,3));
F(:,[1 4]) = real(F(:,[1 4]));

ek = ek(:,:,idxmap);
iG0(1:npt,1:norbsq) = 0;
for ss = 1:4
    if ss==1 || ss==4
        iG0(:,ss) = complex(reshape(repmat(mu - ek(:,:,ss), [1, 1, nwpt]),npt,1), ...
            reshape(repmat(WN, [nkpt, 1]),npt,1));
    else
        iG0(:,ss) = reshape(repmat(-ek(:,:,ss), [1, 1, nwpt]),npt,1);
    end
end

Mtmp = mul2(inv2(G),F);
S = iG0 - inv2(G + mul2(F,conj(Mtmp)));
if solveCase == 1
    P = mul2(Mtmp,conj(S-iG0));
    %P is Hermitian
    P(:,[1 4]) = real(P(:,[1 4]));
    if dim(4)==2
      P = reshape(P(:,[1 2]),[nk1, nk2, nwpt, nel]);
    elseif dim(4)==3
      P = reshape(P(:,[1 2 4]),[nk1, nk2, nwpt, nel]);
    end
else
    P = [];
end

if dim(4)==2
    S = reshape(S(:,[1 2]),[nk1, nk2, nwpt, nel]);
elseif dim(4)==3
    S = reshape(S(:,[1 2 4]),[nk1, nk2, nwpt, nel]);
end



%Change back to (2*Nk(1), 2*Nk(2), 2*Nc, nel)-array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if useFreqSymmetry
    S = S(:,:,[1:end, end:-1:1],:);
    S(:,:,(Nc+1):end,:) = conj(S(:,:,(Nc+1):end,:));
    if solveCase == 1
        P = P(:,:,[1:end, end:-1:1],:);
        P(:,:,(Nc+1):end,:) = conj(P(:,:,(Nc+1):end,:));
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if useSymmetry
    S = S([1:end, end-1:-1:2],[1:end, end-1:-1:2],:,:);
    if solveCase == 1
        P = P([1:end, end-1:-1:2],[1:end, end-1:-1:2],:,:);
    end
end
