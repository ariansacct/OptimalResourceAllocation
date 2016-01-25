function ObjVal = penalty2( x,mydata )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
global M;
global N;
global cost;
M = 10;
N = 3;
% PHR = 0.00005 + (0.00010 - 0.00005) .* rand(1,N);
% E = 15 + (25 - 15) .* rand(M,N);
% CBW = 1 + (4 - 1) .* rand(N,N);
% CR = 0 + (25 - 0) .* rand(M,M);
% PL = 100 + (200 - 100) .* rand(N,N);
% CHR = 0.00015 + (0.00030 - 0.00015) .* rand(N,N);
% Mem = 100 + (200 - 100) .* rand(1,N);
% mem = 1 + (50 - 1) .* rand(1,M);
% L = 100 + (200 - 100) .* rand(1,N);
% l = 1 + (50 - 1) .* rand(1,M);

global PHR;
global E;
global CBW;
global CR;
global PL;
global CHR;
global Mem;
global mem;
global L;
global l;

eta = 0.01;    %   this should be set appropriately


%   mohasebeye cost
C = 0;
C1 = 0;
C2 = 0;
for k=1:N
    for i=1:M
        C1 = C1 + PHR(1,k) .* x(i,k) .* CR(i,k);
    end
end

for p = 1:N-1
    for q = p+1:N
        for i = 1:M-1
            for j = i+1:M
                C2 = C2 + CHR(p,q) .* x(i,p) .* x(j,q) .* CR(i,j) / CBW(p,q);
            end
        end
    end
end

C = C1 + C2;


%   mohasebeye penaltie memory
Pm = 0;
for k = 1:N
    temp = 0;
    for i = 1:M
        temp = temp + mem(1,i) .* x(i,k) - Mem(1,k);
    end
    Pm = Pm + max(0, temp);
end

Pl = 0;
for K = 1:N
    temp = 0;
    for i = 1:M
        temp = temp + l(1,i) .* x(i,k) - L(1,k);
    end
    Pl = Pl + max(0, temp);
end

Pc = 0;
for p = 1:N
    for q = p+1:N
        
        temp = 0;
        for i = 1:M
            for j = 1:M
                temp = temp + x(i,p) .* x(j,q) .* CR(i,j);
            end
        end
        maxsecondarg = temp - PL(p,q);
        
        
        
        
        Pc = Pc + max(0, maxsecondarg);
        
    end
end

cost = C;
ObjVal = C + eta .* (Pm + Pl + Pc);
%ObjVal = C;



end

