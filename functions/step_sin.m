function [ val ] = step_sin( x )

ind_1 = x <=0;
ind_2 = x >= 0;

val = zeros(size(x)); 

val(ind_1) = sin(4*x(ind_1));
val(ind_2) = -sin(4*x(ind_2));

end

