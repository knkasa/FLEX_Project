function [Sut Put mu] = interpSelfE(WNin,varargin)
%INTERPSELFE Interpolate the self-energy on different Matsubara frequency
%grid WN = pi*(2*(-numwi:numwi-1)+1)/beta.

flySelfE = varargin{1};
if flySelfE == 1
    S = varargin{2};
    P = varargin{3};
    WN = varargin{4};
    mu = varargin{5};
else
    dataFileName = varargin{2};
    load(dataFileName,'S','P','WN','mu','-mat')
end
[nk1 nk2 nw nel] = size(S);        %nel>=1
nwin = numel(WNin);
if mod(nw,2)~=0 || mod(nwin,2)~=0
    error('Number of Matsubara frequencies must be even!')
end
if nw~=numel(WN)
    error('Wrong grid for self-energy data!')
end

S = permute(S,[3 1 2 4]);          %now S becomes [nw nk1 nk2 nel]
P = permute(P,[3 1 2 4]);
Sut(1:nwin,1:nk1,1:nk2,1:nel) = 0;
Put(1:nwin,1:nk1,1:nk2,1:nel) = 0;

%1D interpolation with spline-extrapolated values
Sut = interp1(WN,S,WNin,'spline','extrap');
Put = interp1(WN,P,WNin,'spline','extrap');

%Replace spline-extrapolated values by 'nearest' values
idxL = find(WNin<WN(1));
idxR = find(WNin>WN(end));
if ~isempty(idxL)
    idxL = idxL(end);
    Sut(1:idxL,:,:,:) = Sut(repmat(idxL+1,1,idxL),:,:,:);
    Put(1:idxL,:,:,:) = Put(repmat(idxL+1,1,idxL),:,:,:);
end
if ~isempty(idxR)
    idxR = idxR(1);
    Sut(idxR:nwin,:,:,:) = Sut(repmat(idxR-1,1,nwin-idxR+1),:,:,:);
    Put(idxR:nwin,:,:,:) = Put(repmat(idxR-1,1,nwin-idxR+1),:,:,:);
end

Sut = permute(Sut,[2 3 1 4]);      %now Sut becomes [nk1 nk2 nwin nel]
Put = permute(Put,[2 3 1 4]);