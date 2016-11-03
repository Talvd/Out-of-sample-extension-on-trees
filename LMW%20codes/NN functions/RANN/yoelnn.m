m=900;
n=60000;
X=rand(m,n);
nn=5;
idx=int64(zeros(nn,n));
dsts=zeros(nn,n);
numit=int32(2);
stats=zeros(1,100);
iisuper=int32(1);
iistat=int32(0);

rann64(X, idx, dsts, numit, stats, iisuper, iistat);