function [ val ] =  Set_Epsilon_V3(index, distance, test, percent)
%---------------------------------------------------------------
% Automatic set for the epsilon parameter:
%  We measure the total "energy" in the affinity matrix and decide the
%  epsilon due to the linear part of the function 
%             energy = energy(log epsilon)

% Note: this version uses relaxed optimization problem to improve the run time

% Input
%  index,distance - are the info for the nearest neighbors of each point.
%  NN - Number of nearest neighbours to consider in the affinity
%       construction matrix
%
% Output
%  Val - the recomended epsilon
%---------------------------------------------------------------
% Nir Sharon 29-10-12 and June 2013
%---------------------------------------------------------------
defualt_percent = .65;
if nargin<4
    percent = defualt_percent;
end

if nargin<3
    test = 0;
end

if isempty(percent)
        percent = defualt_percent;
else if (percent<0)||(percent>1)
          percent = defualt_percent;
    end
end

%------- get the epsilon behaviour ---------------

exp_lim  = 15; 
exp_step = .5;
number_of_graphs = exp_lim*2+1;

if isempty(nonzeros(distance))
    val = 1;
    return
end

md = median(nonzeros(distance));

%[W,~] = Graph_Matrix(index,distance,1,0);
% assert(~isnan(ep),'Wrong epsilon choice');
% ep is in [0,1], we map this segment to [-exp_lim,0] which is the segment of
% lower exponent of our epsilon

start = -exp_lim; %   floor(exp_lim*(md-1));
x = zeros(number_of_graphs,1);
y = zeros(number_of_graphs,1);

% reconstructing many graphs with different epsi values
for j = start:(start+2*exp_lim);
    current = j-start+1;
    epsi = md*10^(j*exp_step);

  %  epsi = md*jump(current)      %10^(j*exp_step);
    [AffinityMatrix, ~] = Graph_Matrix(index,distance,epsi,0);
    x(current) = epsi;
    y(current) = sum(sum(AffinityMatrix));
end

thr    = 0.05;
delta = y(end) - y(1);

starts = find(y>(thr*delta+y(1)),1,'first');
ends = find(y>((1-thr)*delta+y(1)),1,'first');
rang = starts:ends;
logx = log(x);

% % L_1 minimazation
% tol=1.0e-14;
% maxiter=1000;
% A      =[ones(numel(rang),1) logx(rang)];
% brob = IRLS(A,y(rang),1,maxiter,tol);
% b = brob(1);
% a = brob(2);

% Matlab's LS
p = polyfit(logx(rang),y(rang),1);
a = p(1);
b = p(2);

if a<10*eps 
    x_0 = logx(starts);
    x_1 = logx(ends);
    y_0 =  y(starts);
    y_1 = y(ends);
else
    x_0 = (y(starts)-b)/a;
    x_1 = (y(ends)-b)/a;
    y_0 =  y(starts);
    y_1 = y(ends);
end

dist = norm([x_0-x_1,y_0-y_1 ])*percent;
add     = dist/sqrt(1+a^2);
val   = exp(x_0 + add);




if test
    figure
    loglog(x,y,'.-');    
    hold on;
    scatter(x(rang),y(rang),'r')
    valx =[val,val];
    valy=[y(starts),y(ends)];
    plot(valx,valy,'k--');
    plot(exp([x_0,x_1]),[y_0,y_1],'k')
    legend('full values','Choosen range','Choosen value');
end

end

