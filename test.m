addpath ServiceRoutines
% Driver script for solving the 1D(spherically symmetric reduction) Generalized harmonic Einstein equations
Globals1D;

% Order of polymomials used for approximation 
N = 7;

% Generate simple mesh
[Nv, VX, K, EToV] = MeshGen1D(1.5, 2.5, 1);

% Initialize solver and construct grid and metric
StartUp1D;

% Set initial conditions
A = power(x,-2);
B = power(x,-1);
A_exact = A;
B_exact = B;
time = 0;                                                                                                
FinalTime = 10000.0; 
                                                                                                         
% Runge-Kutta residual storage                                                                           
res_A = zeros(Np,K);                                                                                      
res_B = zeros(Np,K);                                                                                      
                                                                                                         
% compute time step size                                                                                 
xmin = min(abs(x(1,:)-x(2,:)));                                                                          
dt = 0.8/(2*N+1)*xmin;                                                            
Nsteps = ceil(FinalTime/dt); dt = FinalTime/Nsteps;                                                      
time_seq = [];
rhs_A_seq = [];
rhs_B_seq = [];
                                                                                                         
for tstep=1:Nsteps                                                                                   
    for INTRK = 1:5                                                                            
        timelocal = time + rk4c(INTRK)*dt;                                                           
        [rhs_A, rhs_B] = RHS_test(A, B, A, B)
        rhs_A(vmapO) = 0.5*rhs_A(vmapO) + 0.5*rhs_B(vmapO);
        rhs_B(vmapO) = 0.5*rhs_A(vmapO) + 0.5*rhs_B(vmapO);

        res_A = rk4a(INTRK)*res_A + dt*rhs_A;                                                           
        res_B = rk4a(INTRK)*res_B + dt*rhs_B;                                                           

        A = A+rk4b(INTRK)*res_A;                                                                      
        B = B+rk4b(INTRK)*res_B;                                                                      
    end;                                                                                             

    % Increment time                                                                                 
    time = time+dt;                                                                                  
    time_seq = [time_seq, time];
    rhs_A_seq = [rhs_A_seq, max(max(abs(rhs_A)))];
    rhs_B_seq = [rhs_B_seq, max(max(abs(rhs_B)))];
    if (mod(tstep, 100) == 0)
        figure(1); plot(x, A-A_exact); title(['t = ', num2str(time)]); drawnow; pause(.1);
        figure(2); plot(x, B-B_exact); title(['t = ', num2str(time)]); drawnow; pause(.1);
        figure(3); plot(x, rhs_A); title(['rhs\_A, t = ', num2str(time)]); drawnow; pause(.1);
        figure(4); plot(x, rhs_B); title(['rhs\_B, t = ', num2str(time)]); drawnow; pause(.1);
        figure(5); semilogy(time_seq, rhs_A_seq); title(['max of rhs\_A with time']); drawnow; pause(.1);
        figure(6); semilogy(time_seq, rhs_B_seq); title(['max of rhs\_B with time']); drawnow; pause(.1);
    end;
end;  

