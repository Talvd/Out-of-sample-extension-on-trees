function [ idx, dsts] = RANN( X, NN )
% 
% wrap and run the RANN function of Andrei Osipov wi
%
% Nir Sharon, June 2013


if (NN<1)||(mod(NN,1)~=0)
    error('Wrong nearest neighbour parameter');
end

n = size(X,2);
if NN==1
    idx = 1:n;
    dsts= zeros(1,n);
    return;
end

% the rann doesnot return the nodes themselfs in its list
NN = NN-1;

% parameters initializing
idx=int64(zeros(NN,n));
dsts=zeros(NN,n);
numit=int32(5);
stats=zeros(1,100);
iisuper=int32(1);
iistat=int32(0);

% main run
rann64(X, idx, dsts, numit, stats, iisuper, iistat);

idx  = [ 1:n        ; double(idx) ];
dsts = [ zeros(1,n) ; dsts        ];

end

