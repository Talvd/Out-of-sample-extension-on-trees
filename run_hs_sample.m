
load('C:\Users\tal\Dropbox\data\data.mat');
warning('off','all')

NN = 5;
D=4;

run_LMW = 1;
run_LS = 1;
run_KNN =1;
run_MDSE =1;
new_data=1;

R1 = 40:2:89;
R2 = 100:149;  
R3 = 41:2:90;

% R1 = 40:2:89;
% R2 = 50:100;  
% R3 = 41:2:90;

s=size(R1,2)*size(R2,2);
points=zeros(D,s);
for i=1:D
    im=imgs{1}(R1,R2,i);
    points(i,:) =im(:)';
end;

im_f=imgs{1}(R1,R2,10);
f = im_f(:);

rs = size(R3,2)*size(R2,2);
new_points=zeros(D,rs);
for i=1:D
    im=imgs{1}(R3,R2,i);
    new_points(i,:) =im(:)';
end;
im_f=imgs{1}(R3,R2,10);
f_new_org = im_f(:);


%LMW
if (run_LMW)
    k_params = [ 5 10 100];
%     k_params = 1;

    for k=k_params
        f_new =  run_LMW_extension(f, points, new_points, k, D, [], new_data, NN);
        figure_title =['LMW k=' num2str(k) ' hs'];       
        figure_and_error_hs(points, new_points, f, f_new_org, f_new, figure_title, 0);
    end
end

% LS
if (run_LS)
    f_new = [];
    for p=new_points   
        idx= knnsearch(points',p','k',NN);
        v=quaddric_fit_extension(points(:,idx),p,f(idx,:),1);   
        f_new =[f_new; v(end)];   
    end
    figure_and_error_hs(points, new_points, f, f_new_org, f_new, 'LS ', 0);
end

%KNN
if (run_KNN)
    f_new = [];
    for p=new_points    
        idx= knnsearch(points',p','k',NN);    
        M = mean(f(idx));    
        f_new = [f_new; M];       
    end
    figure_and_error_hs(points, new_points, f, f_new_org, f_new, 'KNN ', 0);
end

%MDSE
if (run_MDSE)
    data_mse = [points new_points];
    d = squareform(pdist(data_mse'));

    N=rs;
    in_dists = d(1:N,1:N);
    out_dists = d(N+1:end,1:N);

    [AF S EF errors] = distancesMultiScale(f, in_dists, out_dists, 0.002 ,1);

    f_new_mse = EF(:,size(EF,2)); 
    figure_and_error_hs(points, new_points, f, f_new_org, f_new_mse, 'MDSE ', 0);
end
