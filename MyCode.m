global M
global N
M = 10;
N = 6;
PHR = 0.0005 + (0.00010 - 0.0005) .* rand(1,N);
E = 15 + (25 - 15) .* rand(M,N);
CBW = 1 + (4 - 1) .* rand(N,N);
CR = 0 + (25 - 0) .* rand(M,M);
PL = 100 + (200 - 100) .* rand(N,N);
CHR = 0.00015 + (0.00030 - 0.00015) .* rand(N,N);
Mem = 100 + (200 - 100) .* rand(1,N);
mem = 1 + (50 - 1) .* rand(1,M);
L = 100 + (200 - 100) .* rand(1,N);
l = 1 + (50 - 1) .* rand(1,M);

dim = M;
popsize = 3;

p = randi(N, [popsize,dim]);            % eq. 1
oldp = p;


index = randperm(numel(oldp));  %   for eq. 4
oldp = oldp(index);         % eq. 4 

determine_x(oldp);

F = 3 .* normrnd(0,1);
Mutant = p + F .* (oldp - p);           % eq. 5