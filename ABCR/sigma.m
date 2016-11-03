function [ out ] = sigma( data,k,j,i )
% Computing sigma(i,j) due to data and k parameter
start = 1+(i-1)*k*(2^j);
stop = i*k*(2^j);
out = (data(stop)-data(start))/2;
end
