function [euler] = R1(a)
    euler = [1 0 0; 0 cos(a) sin(a); 0 -sin(a) cos(a)];
end