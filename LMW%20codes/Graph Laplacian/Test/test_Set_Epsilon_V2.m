function [ test ] = test_Set_Epsilon_V2( func_NN )
% Script name: test_Set_Epsilon_V2
%
% test the function "Set_Epsilon_V2"
%
% Nir Sharon, June 2013

N    = 100;
NN   = 10;

%first example - Equidistant on the circle
teta = linspace(0,2*pi,N);
X    = [cos(teta) ; sin(teta)];

[index,distance] = func_NN(X,NN);

Set_Epsilon_V2(index,distance,1 );

%second example - random on the circle
d = 2;
X = rand(d,N);
X = X*diag((sum(X.*X)).^(-.5));
[index,distance] = func_NN(X,NN);
Set_Epsilon_V2(index,distance,1 );

%second example - random on R^d
d = 4;
X = rand(d,N);
[index,distance] = func_NN(X,NN);
test = Set_Epsilon_V2(index,distance,1 );



test = test~=0;

end

