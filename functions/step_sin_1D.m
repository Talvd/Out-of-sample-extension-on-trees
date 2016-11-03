function [ val ] = step_sin_1D( x )

ind_1 = x >=0 & x<= 1;
ind_2 = x >= 1 & x <= 2;

val = zeros(size(x)); 

val(ind_1) = sin(4*x(ind_1));
val(ind_2) = -sin(4*x(ind_2));

end
