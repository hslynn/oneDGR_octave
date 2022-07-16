addpath ServiceRoutines
% Driver script for solving the 1D(spherically symmetric reduction) Generalized harmonic Einstein equations
Globals1D;

% Order of polymomials used for approximation 
N = 1;

% Generate simple mesh
[Nv, VX, K, EToV] = MeshGen1D(1.5, 2.5, 10);

% Initialize solver and construct grid and metric
StartUp1D;

%filter
WITHFILTER = 0;
eta = 0.3;
s = 12;
Nc = floor(eta*N);
F = Filter1D(N, Nc, s);

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
xmin = min(abs(VX(1:end-1)-VX(2:end)));                                                                          
dt = 0.5*xmin/(2*N+1);                                                            
Nsteps = ceil(FinalTime/dt); dt = FinalTime/Nsteps;                                                      
time_seq = [];
rhs_A_seq = [];
rhs_B_seq = [];
err_A_seq = [];
err_B_seq = [];

res_A = zeros(Np,K);                                                                                      
res_B = zeros(Np,K);                                                                                      
                                                                                                         
for tstep=1:Nsteps                                                                                   
    for INTRK = 1:5                                                                            
        timelocal = time + rk4c(INTRK)*dt;                                                           
        [rhs_A, rhs_B] = RHS_test(A, B, A_exact, B_exact);
        res_A = rk4a(INTRK)*res_A + dt*rhs_A;                                                           
        res_B = rk4a(INTRK)*res_B + dt*rhs_B;                                                           

        A = A+rk4b(INTRK)*res_A;                                                                      
        B = B+rk4b(INTRK)*res_B;                                                                      
        if(WITHFILTER)
            A = F*A;
            B = F*B;
        end;
    end;                                                                                             

    % Increment time                                                                                 
    time = time+dt;                                                                                  
    time_seq = [time_seq, time];
    rhs_A_seq = [rhs_A_seq, L2norm(rhs_A)];
    rhs_B_seq = [rhs_B_seq, L2norm(rhs_B)];
    err_A_seq = [err_A_seq, L2norm(A-A_exact)];
    err_B_seq = [err_B_seq, L2norm(B-B_exact)];
    if (mod(tstep, 100) == 0)
        figure(1); 
        subplot(2,4,1);
        plot(x, A-A_exact); title(['error of A, t = ', num2str(time)]); drawnow;
        subplot(2,4,2);
        plot(x, B-B_exact); title(['error of B, t = ', num2str(time)]); drawnow;
        subplot(2,4,3);
        plot(x, rhs_A); title(['rhs\_A, t = ', num2str(time)]); drawnow;
        subplot(2,4,4);
        plot(x, rhs_B); title(['rhs\_B, t = ', num2str(time)]); drawnow;
        subplot(2,4,5);
        semilogy(time_seq, rhs_A_seq); title(['L2 of rhs\_A with time']); drawnow;
        subplot(2,4,6);
        semilogy(time_seq, rhs_B_seq); title(['L2 of rhs\_B with time']); drawnow;
        subplot(2,4,7);
        semilogy(time_seq, err_A_seq); title(['L2 of error of A with time']); drawnow;
        subplot(2,4,8);
        semilogy(time_seq, err_B_seq); title(['L2 of error of B with time']); drawnow;
    end;
end;  
