function ObjVal = penalty( x,mydata )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
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

eta = 1;    %   this should be set appropriately

% syms i j k p q
% C = symsum( symsum( PHR(1,k) .* x(i,k) .* CR(i,k) , 1 , M ) , k , 1 , N ) + symsum( symsum( symsum( symsum( CHR(p,q) .* x(i,p) .* x(j,q) .* CR(i,j) / CBW(p,q) , j , i+1 , M ) , i , 1 , M-1 ) , q , p+1 , N ) , p , 1 , N-1 );
% Pm = symsum( max( 0, sumsum( mem(1,i) .* x(i,k) - Mem(1,k) , i , 1 , M ) , k , 1 , N ) );
% Pl = symsum( max( 0, sumsum( l(1,i) .* x(i,k) - L(1,k) , i , 1 , M ) , k , 1 , N ) );
% Pc = symsum( symsum( max( 0, symsum( symsum( x(i,p) .* x(j,q) .* CR(i,j) , j , 1 , M ) , i , 1 , M ) - PL(p,q) ) , q , p+1 , N ) , p , 1 , N );

%   lets replace i by s and j by t
syms s t k p q
C = symsum( symsum( PHR(1,k) .* x(s,k) .* CR(s,k) , 1 , M ) , k , 1 , N ) + symsum( symsum( symsum( symsum( CHR(p,q) .* x(s,p) .* x(t,q) .* CR(s,t) / CBW(p,q) , t , s+1 , M ) , s , 1 , M-1 ) , q , p+1 , N ) , p , 1 , N-1 );
Pm = symsum( max( 0, sumsum( mem(1,s) .* x(s,k) - Mem(1,k) , s , 1 , M ) , k , 1 , N ) );
Pl = symsum( max( 0, sumsum( l(1,s) .* x(s,k) - L(1,k) , s , 1 , M ) , k , 1 , N ) );
Pc = symsum( symsum( max( 0, symsum( symsum( x(s,p) .* x(t,q) .* CR(s,t) , t , 1 , M ) , s , 1 , M ) - PL(p,q) ) , q , p+1 , N ) , p , 1 , N );

ObjVal = C + eta .* (Pm + Pl + Pc);

end

