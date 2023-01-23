function [euler] = R2(a)
    euler = [cos(a) 0 -sin(a); 0 1 0; sin(a) 0 cos(a)];
end