function determine_x(p)
global M;
global N;
global X;
X = zeros(M,N);
for i = 1:numel(p)
    X(i,p(i)) = 1;
end