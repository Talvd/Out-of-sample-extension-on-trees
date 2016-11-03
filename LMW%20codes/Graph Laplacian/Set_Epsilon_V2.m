function [ val ] =  Set_Epsilon_V2(index, distance, test, percent)
%---------------------------------------------------------------
% Automatic set for the epsilon parameter:
%  We measure the total "energy" in the affinity matrix and decide the
%  epsilon due to the linear part of the function 
%             energy = energy(log epsilon)
%
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

if nargin<4
    percent = .7;
end

if nargin<3
    test = 0;
end

if isempty(percent)
        percent = .7;
else if (percent<0)||(percent>1)
          percent = .7;
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
%z = zeros(number_of_graphs,1);

for j = start:(start+2*exp_lim);
    current = j-start+1;
    epsi = md*10^(j*exp_step);

  %  epsi = md*jump(current)      %10^(j*exp_step);
    [AffinityMatrix, ~] = Graph_Matrix(index,distance,epsi,0);
    x(current) = epsi;
    y(current) = sum(sum(AffinityMatrix));
%     if j>start
%         z(current) = (y(current)-y(current-1))/log(exp_step); %(log(x(current))-x(current-1));
%     end
end

%------- get the linear part ---------------
min = inf;
j_s = 1;
i_s = number_of_graphs;
f = zeros(number_of_graphs,1);
logx = log(x);
for j=1:number_of_graphs-1
    for i=j+1:number_of_graphs
        % the constant part
        f(1:(j-1)) = y(1)  ;
        f(i+1:end) = y(end);
        % the linear part
        p      = polyfit(logx(j:i),y(j:i),1);
        f(j:i) = polyval(p,logx(j:i));
        e      = norm(f-y);
        % update the minimum
        if e<min
            min = e;
            j_s = j;
            i_s = i;
        end
    end
end

%---------------- summarize ------------------
ind = floor( (1-percent)*j_s + percent* i_s);
val = x(ind);

if test
    figure
    loglog(x,y,'.-');    
    hold on;
    scatter(x(j_s:i_s),y(j_s:i_s),'r')
    scatter(x(ind),y(ind),'filled','k');
     legend('full values','sub domain','Choosen value');
end

end

