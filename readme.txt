Backtracking Search Algorithm and its source code courtesy of Pinar Civicioglu of Erciyes University, Turkey.


Sample execution:
FncNumber=1; % FncNumber={1=griewank, 2=rastrigin, 3=rosenbrock, 4=FoxHoles, 5=ackley}
settingOfBenchmarkFnc(FncNumber); bsa3(fnc,[],30,dim,1,low,up,1e6)


BSA|    1 -----> 520.5563288824233700
BSA|    2 -----> 520.5563288824233700
BSA|    3 -----> 520.5563288824233700
...
BSA| 3226 -----> 0.0000000000000001
BSA| 3227 -----> 0.0000000000000000
