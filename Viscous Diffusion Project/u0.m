% Fully developed Poiseuille flow through rectangular duct.
% Will be the initial velocity profile.
% Equation 3-48 of the textbook, non-dimensionalized.
% Assumes a = 1.

function thing = u0(x, y, b)
N = 100; % total number of terms
% If too many terms are used, the cosh() will be undefined.

sum = 0;
for it = 1:2:2*N-1
    neg = (-1)^((it-1)/2);
    C = cos(it*pi*x/2)/(it^3);
    H = 1 - cosh(it*pi*y/2)/cosh(it*pi*b/2);
    sum = sum + neg*C*H;
end
thing = sum;
end


