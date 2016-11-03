function [ out ] = mu( data,k,j,i )
% Computing mu(i,j) due to data and k parameter
start = 1+(i-1)*k*(2^j);
stop = i*k*(2^j);
out = (data(start)+data(stop))/2;
end

