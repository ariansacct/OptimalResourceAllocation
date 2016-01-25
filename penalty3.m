function ObjVal = penalty3( x,mydata )
%UNTITLED2 Summary of this function goes here
%   esme function be eshtebah penalty ast; esmesh bayad totalcost bashe
%   Detailed explanation goes here
global M;
global N;
global cost;
% M = 10;
% N = 3;
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
global L;   %   same as total processor available
global l;   %   same as processor demand
global hwhr;
global vmmhr;
global X;

global ST;
global st;

eta = 0.01;    %   this should be set appropriately


%   mohasebeye cost
value=0;
for j=1:N
%     disp('salam goosi');
% findomega(j)
% disp('salam toole');
	value = value + findomega(j).*( hwhr(j)+vmmhr(j) );
end

%   mohasebeye penaltie memory
Pm = 0;
for k = 1:N
    temp = 0;
    for i = 1:M
        temp = temp + mem(1,i).* X(i,k) - Mem(1,k);
    end
    Pm = Pm + max(0, temp);
end

%   mohasebeye penaltie processor (load)
Pl = 0;
for K = 1:N
    temp = 0;
    for i = 1:M
        temp = temp + l(1,i).* X(i,k) - L(1,k);
    end
    Pl = Pl + max(0, temp);
end

% Pc = 0;
% for p = 1:N
%     for q = p+1:N
%         
%         temp = 0;
%         for i = 1:M
%             for j = 1:M
%                 temp = temp + X(i,p) .* X(j,q) .* CR(i,j);
%             end
%         end
%         maxsecondarg = temp - PL(p,q);
%         
%         
%         
%         
%         Pc = Pc + max(0, maxsecondarg);
%         
%     end
% end



%   mohasebeye penaltie storage
Ps = 0;
for K = 1:N
    temp = 0;
    for i = 1:M
        temp = temp + st(1,i).* X(i,k) - ST(1,k);
    end
    Ps = Ps + max(0, temp);
end






ObjVal = value + eta .* (Pm + Pl + Ps);
cost = value;
%ObjVal = C;

end

function omega=findomega(j)
global X;
% X=[1 0;1 0;1 0];
omega=max( max( 0,X(:,j) ) );
end