addpath ServiceRoutines
% Driver script for solving the 1D(spherically symmetric reduction) Generalized harmonic Einstein equations
Globals1D;

% Order of polymomials used for approximation 
N = 1;

% Generate simple mesh
[Nv, VX, K, EToV] = MeshGen1D(1.99, 2.01, 1);

% Initialize solver and construct grid and metric
StartUp1D;

% Set initial conditions
%init_min;
init_schwarzschild_kerr_schild;
%init_PG;
%init_test;
time = 0;                                                                                                
FinalTime = 10000.0; 
time_seq = [];
rhs_g11_seq = [];
rhs_Pi11_seq = [];
rhs_Phi11_seq = [];
C0_seq = [];
C1_seq = [];
Cr11_seq = [];
                                                                                                         
% Runge-Kutta residual storage                                                                           
res_g00 = zeros(Np,K);                                                                                      
res_g01 = zeros(Np,K);                                                                                      
res_g11 = zeros(Np,K);                                                                                      
res_Pi00 = zeros(Np,K);                                                                                      
res_Pi01 = zeros(Np,K);                                                                                      
res_Pi11 = zeros(Np,K);                                                                                      
res_Phi00 = zeros(Np,K);                                                                                      
res_Phi01 = zeros(Np,K);                                                                                      
res_Phi11 = zeros(Np,K);                                                                                      
                                                                                                         
% compute time step size                                                                                 
xmin = min(abs(x(1,:)-x(2,:)));                                                                          
dt = 0.9/(2*N+1)*xmin;                                                            
Nsteps = ceil(FinalTime/dt); dt = FinalTime/Nsteps;                                                      

[rhs_g00, rhs_g01, rhs_g11, rhs_Pi00, rhs_Pi01, rhs_Pi11, rhs_Phi00, rhs_Phi01, rhs_Phi11] = compute_RHS
                                                                                                         
for tstep=1:Nsteps                                                                                   
    for INTRK = 1:5                                                                            
        timelocal = time + rk4c(INTRK)*dt;                                                           
        [rhs_g00, rhs_g01, rhs_g11, rhs_Pi00, rhs_Pi01, rhs_Pi11, rhs_Phi00, rhs_Phi01, rhs_Phi11] = compute_RHS;

        b = sqrt(1/g11(vmapO));
        rhs_Pi00(vmapO) = paragamma2/2*rhs_g00(vmapO) + 1/2*rhs_Pi00(vmapO) ...
                + b/2*rhs_Phi00(vmapO) ;
        rhs_Pi01(vmapO) = paragamma2/2*rhs_g01(vmapO) + 1/2*rhs_Pi01(vmapO) ...
                + b/2*rhs_Phi01(vmapO) ;
        rhs_Pi11(vmapO) = paragamma2/2*rhs_g11(vmapO) + 1/2*rhs_Pi11(vmapO) ...
                + b/2*rhs_Phi11(vmapO) ;
        rhs_Phi00(vmapO) = -paragamma2/2/b*rhs_g00(vmapO) ...
                + 1/2/b*rhs_Pi00(vmapO) + 1/2*rhs_Phi00(vmapO) ;
        rhs_Phi01(vmapO) = -paragamma2/2/b*rhs_g01(vmapO) ...
                + 1/2/b*rhs_Pi01(vmapO) + 1/2*rhs_Phi01(vmapO) ;
        rhs_Phi11(vmapO) = -paragamma2/2/b*rhs_g11(vmapO) ...
                + 1/2/b*rhs_Pi11(vmapO) + 1/2*rhs_Phi11(vmapO) ;

        %rhs_Pi00(vmapO) = paragamma2*rhs_g00(vmapO);
        %rhs_Pi01(vmapO) = paragamma2*rhs_g01(vmapO);
        %rhs_Pi11(vmapO) = paragamma2*rhs_g11(vmapO);
        %rhs_Phi00(vmapO) = 0;
        %rhs_Phi01(vmapO) = 0;
        %rhs_Phi11(vmapO) = 0;


        res_g00 = rk4a(INTRK)*res_g00 + dt*rhs_g00;                                                           
        res_g01 = rk4a(INTRK)*res_g01 + dt*rhs_g01;                                                           
        res_g11 = rk4a(INTRK)*res_g11 + dt*rhs_g11;                                                           
        res_Pi00 = rk4a(INTRK)*res_Pi00 + dt*rhs_Pi00;                                                           
        res_Pi01 = rk4a(INTRK)*res_Pi01 + dt*rhs_Pi01;                                                           
        res_Pi11 = rk4a(INTRK)*res_Pi11 + dt*rhs_Pi11;                                                           
        res_Phi00 = rk4a(INTRK)*res_Phi00 + dt*rhs_Phi00;                                                           
        res_Phi01 = rk4a(INTRK)*res_Phi01 + dt*rhs_Phi01;                                                           
        res_Phi11 = rk4a(INTRK)*res_Phi11 + dt*rhs_Phi11;                                                           

        g00 = g00+rk4b(INTRK)*res_g00;                                                                      
        g01 = g01+rk4b(INTRK)*res_g01;                                                                      
        g11 = g11+rk4b(INTRK)*res_g11;                                                                      
        Pi00 = Pi00+rk4b(INTRK)*res_Pi00;                                                                      
        Pi01 = Pi01+rk4b(INTRK)*res_Pi01;                                                                      
        Pi11 = Pi11+rk4b(INTRK)*res_Pi11;                                                                      
        Phi00 = Phi00+rk4b(INTRK)*res_Phi00;                                                                      
        Phi01 = Phi01+rk4b(INTRK)*res_Phi01;                                                                      
        Phi11 = Phi11+rk4b(INTRK)*res_Phi11;                                                                      
    end;                                                                                             

    % limiter
    %g00 = SlopeLimitN(g00);
    %g01 = SlopeLimitN(g01);
    %g11 = SlopeLimitN(g11);
    %Pi00 = SlopeLimitN(Pi00);
    %Pi01 = SlopeLimitN(Pi01);
    %Pi11 = SlopeLimitN(Pi11);
    %Phi00 = SlopeLimitN(Phi00);
    %Phi01 = SlopeLimitN(Phi01);
    %Phi11 = SlopeLimitN(Phi11);

    %filter
    %F = Filter1D(N, 0, 50);
    %g00 = F*g00;
    %g01 = F*g01;
    %g11 = F*g11;
    %Pi00 = F*Pi00;
    %Pi01 = F*Pi01;
    %Pi11 = F*Pi11;
    %Phi00 = F*Phi00;
    %Phi01 = F*Phi01;
    %Phi11 = F*Phi11;

    % Increment time                                                                                 
    time = time+dt;                                                                                  
    time_seq = [time_seq, time];
    rhs_g11_seq = [rhs_g11_seq, max(max(abs(rhs_g11)))];
    rhs_Pi11_seq = [rhs_Pi11_seq, max(max(abs(rhs_Pi11)))];
    rhs_Phi11_seq = [rhs_Phi11_seq, max(max(abs(rhs_Phi11)))];

    C0_seq = [C0_seq, max(max(abs(C0)))];
    C1_seq = [C1_seq, max(max(abs(C1)))];
    Cr11_seq = [Cr11_seq, max(max(abs(Cr11)))];
    if (mod(tstep, 1000) == 0)
        %figure(1); plot(x, g00-g00_exact); title(['Error of g00, t = ', num2str(time)]); drawnow; pause(.1);
        %figure(2); plot(x, g01-g01_exact); title(['Error of g01, t = ', num2str(time)]); drawnow; pause(.1);
        figure(3); plot(x, g11-g11_exact); title(['Errof of g11, t = ', num2str(time)]); drawnow; pause(.1);
        %figure(4); plot(x, Pi00-Pi00_exact); title(['t = ', num2str(time)]); drawnow; pause(.1);
        %figure(5); plot(x, Pi01-Pi01_exact); title(['t = ', num2str(time)]); drawnow; pause(.1);
        figure(6); plot(x, Pi11-Pi11_exact); title(['Error of Pi11, t = ', num2str(time)]); drawnow; pause(.1);
        %figure(7); plot(x, Phi00-Phi00_exact); title(['t = ', num2str(time)]); drawnow; pause(.1);
        %figure(8); plot(x, Phi01-Phi01_exact); title(['t = ', num2str(time)]); drawnow; pause(.1);
        figure(9); plot(x, Phi11-Phi11_exact); title(['Error of Phi11, t = ', num2str(time)]); drawnow; pause(.1);
        figure(22); plot(x, g11.*g00 - g01.*g01); title(['Det(g), t = ', num2str(time)]); drawnow; pause(.1);

        %figure(10); plot(x, rhs_g11); title(['rhs\_g11, t = ', num2str(time)]); drawnow; pause(.1);
        %figure(11); plot(x, rhs_Pi11); title(['rhs\_Pi11, t = ', num2str(time)]); drawnow; pause(.1);
        %figure(12); plot(x, rhs_Phi11); title(['rhs\_Phi11, t = ', num2str(time)]); drawnow; pause(.1);

        figure(13); semilogy(time_seq, C0_seq); title(['max of C0 with time']); drawnow; pause(.1);
        figure(14); semilogy(time_seq, C1_seq); title(['max of C1 with time']); drawnow; pause(.1);
        figure(15); semilogy(time_seq, Cr11_seq); title(['max of Cr11 with time']); drawnow; pause(.1);

        figure(16); semilogy(time_seq, rhs_g11_seq); title(['max of rhs\_g11 with time']); drawnow; pause(.1);
        figure(17); semilogy(time_seq, rhs_Pi11_seq); title(['max of rhs\_Pi11 with time']); drawnow; pause(.1);
        figure(18); semilogy(time_seq, rhs_Phi11_seq); title(['max of rhs\_Phi11 with time']); drawnow; pause(.1);

    end;
end;  

