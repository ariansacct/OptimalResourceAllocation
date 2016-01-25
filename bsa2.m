function bsa2(fnc,mydata,popsize,dim,DIM_RATE,low,up,epoch)

global M;
global N;
global cost;
global X;
M = 3;
N = 2;
dim = M;
%popsize = 3;   we determine popsize in the call to bsa()

low = 1;
up = N;

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

global hwhr;
global vmmhr;

global ST;
global st;

PHR = 0.00005 + (0.00010 - 0.00005) .* rand(1,N);
E = 15 + (25 - 15) .* rand(M,N);
CBW = 1 + (4 - 1) .* rand(N,N);
CR = 0 + (25 - 0) .* rand(M,M);
PL = 100 + (200 - 100) .* rand(N,N);
CHR = 0.00015 + (0.00030 - 0.00015) .* rand(N,N);

% Mem = 100 + (200 - 100) .* rand(1,N);
% mem = 1 + (50 - 1) .* rand(1,M);
Mem = [1000 1800];
mem=[300 500 1500];

L = 100 + (200 - 100) .* rand(1,N);
l = 1 + (50 - 1) .* rand(1,M);

ST = 100 + (200 - 100) .* rand(1,N);    %   FELAN GOZASHTIMESH MESLE MEMORY; IN BAYAD AVAZ SHE
st = 1 + (50 - 1) .* rand(1,M);         %   FELAN GOZASHTIMESH MESLE MEMORY; IN BAYAD AVAZ SHE

hwhr = 0.00005 + (0.00010 - 0.00005) .* rand(1,N);  %   felan gozashtimesh mesle processor
vmmhr = 0.00005 + (0.00010 - 0.00005) .* rand(1,N); %   felan gozashtimesh mesle processor


%INITIALIZATION
if numel(low)==1, low=low*ones(1,dim); up=up*ones(1,dim); end % this line must be adapted to your problem
pop=GeneratePopulation(popsize,dim,low,up); % see Eq.1 in [1]
determine_x(pop);
fitnesspop=feval(fnc,pop,mydata);
historical_pop=GeneratePopulation(popsize,dim,low,up); % see Eq.2 in [1]

% historical_pop  is swarm-memory of BSA as mentioned in [1].

% ------------------------------------------------------------------------------------------

for epk=1:epoch
    %SELECTION-I
    if rand<rand, historical_pop=pop; end  % see Eq.3 in [1]
    historical_pop=historical_pop(randperm(popsize),:); % see Eq.4 in [1]
    determine_x(pop);
    F=get_scale_factor; % see Eq.5 in [1], you can other F generation strategies
    map=zeros(popsize,dim); % see Algorithm-2 in [1]
    if rand<rand,
        for i=1:popsize,  u=randperm(dim); map(i,u(1:ceil(DIM_RATE*rand*dim)))=1; end
    else
        for i=1:popsize,  map(i,randi(dim))=1; end
    end
    % RECOMBINATION (MUTATION+CROSSOVER)
    offsprings=round( pop+(map.*F).*(historical_pop-pop) );   % see Eq.5 in [1]    must be integer
    offsprings=BoundaryControl(offsprings,low,up); % see Algorithm-3 in [1]
    determine_x(offsprings);
    
    % SELECTON-II
    fitnessoffsprings=feval(fnc,offsprings,mydata);
    ind=fitnessoffsprings<fitnesspop;
%     disp('salam salam salam');
%     fitnessoffsprings==fitnesspop
%     disp('khodafez khodafez khodafez');
    fitnesspop(ind)=fitnessoffsprings(ind);
    pop(ind,:)=offsprings(ind,:);
    [globalminimum,ind]=min(fitnesspop);
    globalminimizer=pop(ind,:);
    % EXPORT SOLUTIONS
    assignin('base','globalminimizer',globalminimizer);
    assignin('base','globalminimum',globalminimum);
%         fprintf('BSA|%5.0f -----> %9.16f\n',epk,globalminimum);
    fprintf('BSA|%5.0f -----> %9.16f\n',epk,cost);
    
end
return

function pop=GeneratePopulation(popsize,dim,low,up)
pop=ones(popsize,dim);
for i=1:popsize
    for j=1:dim
        pop(i,j)=round( rand*(up(j)-low(j))+low(j) );   %   must be integer     floor GHALATE!
    end
end
return

function pop=BoundaryControl(pop,low,up)
[popsize,dim]=size(pop);
for i=1:popsize
    for j=1:dim
        k=rand<rand; % you can change boundary-control strategy
        if pop(i,j)<low(j), if k, pop(i,j)=low(j) ;  else pop(i,j)=round( rand*(up(j)-low(j))+low(j) );  end, end    % must be integer
        if pop(i,j)>up(j),  if k, pop(i,j)=up(j);  else pop(i,j)=round( rand*(up(j)-low(j))+low(j) );  end, end    % must be integer
    end
end
return



function F=get_scale_factor % you can change generation strategy of scale-factor,F
%F=3*randn; % STANDARD brownian-walk
%F=4*randn;  % brownian-walk
% F=4*randg;  % brownian-walk
%F=lognrnd(rand,5*rand);  % brownian-walk
% F=1/normrnd(0,5);        % pseudo-stable walk (levy-like)
F=1./gamrnd(1,0.5);      % pseudo-stable walk (levy-like, simulates inverse gamma distribution; levy-distiribution)
% F=2.*rand^2;    % male man!
return