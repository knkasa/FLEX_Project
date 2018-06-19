function [G F] = solve_dyson(S,P,WN,ek,mu,norb,Nk,Nc,useSymmetry,varargin)
%SOLVE_DYSON Solve Dyson's equation to find normal (G) and anomalous (F)
%Green's function from normal (S) and anomalous (P) self-energy.


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


siz_S = size(S);
siz_e = size(ek);
if ~isequal(siz_S(1:3),[2*Nk(1) 2*Nk(2) 2*Nc])
    error(' Input array S has a wrong size.')
end
if ~isequal(siz_e(1:2),[2*Nk(1) 2*Nk(2)])
    error(' Input array ek has a wrong size.')
end
if size(WN,1)~=1
    WN = reshape(WN,1,[]);
end
nel = siz_S(4);     %number of elements of the matrix in orbital space stored
if ~isequal(siz_e(3),nel)
    error(' Input array ek or S has a wrong size.')
end

if useSymmetry
    S = S(1:Nk(1)+1,1:Nk(2)+1,:,:);
    P = P(1:Nk(1)+1,1:Nk(2)+1,:,:);
    ek = ek(1:Nk(1)+1,1:Nk(2)+1,:);
    nk1 = Nk(1)+1;
    nk2 = Nk(2)+1;
else
    nk1 = 2*Nk(1);
    nk2 = 2*Nk(2);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if useFreqSymmetry
    S = S(:,:,1:Nc,:);
    P = P(:,:,1:Nc,:);
    WN = WN(1:Nc);
    nwpt = Nc;
else
    nwpt = 2*Nc;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nkpt = nk1*nk2;
npt = nkpt*nwpt;     % nk*numwn  into single column
norbsq = norb*norb;  % #orb*#orb=4  4x4 matrix for chi


dim = size(S);
if dim(4)==2
    idxmap = [1 2 2 1];      %eps(1,1) is equivalent to eps(2,2)
elseif dim(4)==3
    idxmap = [1 2 2 3];      %eps(1,1) is not equivalent to eps(2,2)
end

% Note: Self-energy in order S(1,1) S(1,2) S(2,1) S(2,2) 
S = S(:,:,:,idxmap);  % making it 4 component vector
P = P(:,:,:,idxmap);
S = reshape(S,npt,norbsq);   %now S becomes [nk1*nk2*nw, norbsq]
P = reshape(P,npt,norbsq);
%P needs to be Hermitian
P(:,3) = conj(P(:,3));
P(:,[1 4]) = real(P(:,[1 4]));

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

% Getting Green's function
Mtmp = mul2(inv2(S-iG0),P);
G = inv2(iG0 - S - mul2(P,conj(Mtmp)));
if solveCase == 1
    F = mul2(Mtmp,conj(G));
    %F is Hermitian
    F(:,[1 4]) = real(F(:,[1 4]));
    if dim(4)==2
      F = reshape(F(:,[1 2]),[nk1, nk2, nwpt, nel]);   %if eps(1,1)=eps(2,2)
    elseif dim(4)==3
      F = reshape(F(:,[1 2 4]),[nk1, nk2, nwpt, nel]); %if eps(1,1)~=eps(2,2)
    end    
else
    F = [];
end

if dim(4)==2
    G = reshape(G(:,[1 2]),[nk1, nk2, nwpt, nel]);
elseif dim(4)==3
    G = reshape(G(:,[1 2 4]),[nk1, nk2, nwpt, nel]);
end
 


%Change back to (2*Nk(1), 2*Nk(2), 2*Nc, nel)-array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if useFreqSymmetry
    G = G(:,:,[1:end, end:-1:1],:);
    G(:,:,(Nc+1):end,:) = conj(G(:,:,(Nc+1):end,:));
    if solveCase == 1
        F = F(:,:,[1:end, end:-1:1],:);
        F(:,:,(Nc+1):end,:) = conj(F(:,:,(Nc+1):end,:));
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if useSymmetry
    G = G([1:end, end-1:-1:2],[1:end, end-1:-1:2],:,:);
    if solveCase == 1
        F = F([1:end, end-1:-1:2],[1:end, end-1:-1:2],:,:);
    end
end
