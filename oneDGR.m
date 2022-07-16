addpath ServiceRoutines
% Driver script for solving the 1D(spherically symmetric reduction) Generalized harmonic Einstein equations
Globals1D;
GlobalsGR;

% Order of polymomials used for approximation
N = 2; 

% Generate simple mesh
inB = 1.8;
outB = 10;
meshNum = 100;
[Nv, VX, K, EToV] = MeshGen1D(inB, outB, meshNum);

fprintf("N=%d\n", N);
fprintf("mesh numbers=%d\n", meshNum);

% Initialize solver and construct grid and metric
StartUp1D;

% Limiter and filter
WITHLIMITER = 0;
WITHFILTER = 0;
if (WITHFILTER)
    eta = 0.5;
    s = 32; 
    Nc = floor(eta*N);
    F = Filter1D(Nc, s);
end

% compute time step size                                                                                 
time = 0;                                                                                                
FinalTime = 1000.0; 
xmin = min(abs(VX(1:end-1)-VX(2:end)));                                                                          
dt = 0.8/(2*N+1)*xmin;      
Nsteps = FinalTime/dt;

% Set initial conditions
%init_min;
%init_iso;
%init_PG_out;
%init_PG_in;
%init_ks_in;
%init_ks_out;
init_min_random;

%Set bdry condition
DIRICHLET = 1;
FREEZING = 2;
bdry_type = 1; 

%Evolve S or psi
USE_S = 0;
USE_psi = 0;

%lists to store results
time_seq = [];
rhs_g11_seq = [];
rhs_Pi11_seq = [];
rhs_Phi11_seq = [];

err_g11_seq = [];
err_Pi11_seq = [];
err_Phi11_seq = [];

C0_seq = [];
C1_seq = [];

Cr00_seq = [];
Cr01_seq = [];
Cr11_seq = [];

AH_seq = [];
                                                                                                         
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
res_S = zeros(Np,K);                                                                                                    
res_Pi_S = zeros(Np,K);
res_Phi_S = zeros(Np,K);
res_psi = zeros(Np,K);
res_Pi_psi = zeros(Np,K);
res_Phi_psi = zeros(Np,K);

tic;
for tstep=1:Nsteps                                                                                   
    for INTRK = 1:5                                                                            
        timelocal = time + rk4c(INTRK)*dt;                                                           
        compute_RHS;
        rhs_g00(1) = 0.0;
        rhs_g01(1) = 0.0;
        rhs_g11(1) = 0.0;
        rhs_Pi00(1) = 0.0;
        rhs_Pi01(1) = 0.0;
        rhs_Pi11(1) = 0.0;
        rhs_Phi00(1) = 0.0;
        rhs_Phi01(1) = 0.0;
        rhs_Phi11(1) = 0.0;
        rhs_S(1) = 0.0;
        rhs_Pi_S(1) = 0.0;
        rhs_Phi_S(1) = 0.0;


        rhs_g00(end) = 0.0;
        rhs_g01(end) = 0.0;
        rhs_g11(end) = 0.0;
        rhs_Pi00(end) = 0.0;
        rhs_Pi01(end) = 0.0;
        rhs_Pi11(end) = 0.0;
        rhs_Phi00(end) = 0.0;
        rhs_Phi01(end) = 0.0;
        rhs_Phi11(end) = 0.0;
        rhs_S(end) = 0.0;
        rhs_Pi_S(end) = 0.0;
        rhs_Phi_S(end) = 0.0;


        if (bdry_type == FREEZING)
            apply_freezing;
        end

        res_g00 = rk4a(INTRK)*res_g00 + dt*rhs_g00;                                                           
        res_g01 = rk4a(INTRK)*res_g01 + dt*rhs_g01;                                                           
        res_g11 = rk4a(INTRK)*res_g11 + dt*rhs_g11;                                                           
        res_Pi00 = rk4a(INTRK)*res_Pi00 + dt*rhs_Pi00;                                                           
        res_Pi01 = rk4a(INTRK)*res_Pi01 + dt*rhs_Pi01;                                                           
        res_Pi11 = rk4a(INTRK)*res_Pi11 + dt*rhs_Pi11;                                                           
        res_Phi00 = rk4a(INTRK)*res_Phi00 + dt*rhs_Phi00;                                                           
        res_Phi01 = rk4a(INTRK)*res_Phi01 + dt*rhs_Phi01;                                                           
        res_Phi11 = rk4a(INTRK)*res_Phi11 + dt*rhs_Phi11;                                                           
        if (USE_S)
            res_S = rk4a(INTRK)*res_S + dt*rhs_S;                                                           
            res_Pi_S = rk4a(INTRK)*res_Pi_S + dt*rhs_Pi_S;                                                           
            res_Phi_S = rk4a(INTRK)*res_Phi_S + dt*rhs_Phi_S;                                                           
        end
        if (USE_psi)
            res_psi = rk4a(INTRK)*res_psi + dt*rhs_psi;                                                           
            res_Pi_psi = rk4a(INTRK)*res_Pi_psi + dt*rhs_Pi_psi;                                                           
            res_Phi_psi = rk4a(INTRK)*res_Phi_psi + dt*rhs_Phi_psi;                                                           
        end

        g00 = g00+rk4b(INTRK)*res_g00;                                                                      
        g01 = g01+rk4b(INTRK)*res_g01;                                                                      
        g11 = g11+rk4b(INTRK)*res_g11;                                                                      
        Pi00 = Pi00+rk4b(INTRK)*res_Pi00;                                                                      
        Pi01 = Pi01+rk4b(INTRK)*res_Pi01;                                                                      
        Pi11 = Pi11+rk4b(INTRK)*res_Pi11;                                                                      
        Phi00 = Phi00+rk4b(INTRK)*res_Phi00;                                                                      
        Phi01 = Phi01+rk4b(INTRK)*res_Phi01;                                                                      
        Phi11 = Phi11+rk4b(INTRK)*res_Phi11;                                                                      
        if (USE_S)
            S = S+rk4b(INTRK)*res_S;                                                                      
            Pi_S = Pi_S+rk4b(INTRK)*res_Pi_S;                                                                      
            Phi_S = Phi_S+rk4b(INTRK)*res_Phi_S;                                                                      
        end
        if (USE_psi)
            psi = psi+rk4b(INTRK)*res_psi;                                                                      
            Pi_psi = Pi_psi+rk4b(INTRK)*res_Pi_psi;                                                                      
            Phi_psi = Phi_psi+rk4b(INTRK)*res_Phi_psi;                                                                      
        end

        
    end                                                                                           
    %filter
    if (WITHFILTER)
        g00 = F*g00;
        g01 = F*g01;
        g11 = F*g11;
        Pi00 = F*Pi00;
        Pi01 = F*Pi01;
        Pi11 = F*Pi11;
        Phi00 = F*Phi00;
        Phi01 = F*Phi01;
        Phi11 = F*Phi11;
        %S = F*S;
        %Pi_S = F*Pi_S;
        %Phi_S = F*Phi_S;
        %psi = F*psi;
        %Pi_psi = F*Pi_psi;
        %Phi_psi = F*Phi_psi;
    end

    % limiter
    if (WITHLIMITER)
        g00 = SlopeLimit1(g00);
        g01 = SlopeLimit1(g01);
        g11 = SlopeLimit1(g11);
        Pi00 = SlopeLimit1(Pi00);
        Pi01 = SlopeLimit1(Pi01);
        Pi11 = SlopeLimit1(Pi11);
        Phi00 = SlopeLimit1(Phi00);
        Phi01 = SlopeLimit1(Phi01);
        Phi11 = SlopeLimit1(Phi11);
        S = SlopeLimit1(S);
        Pi_S = SlopeLimit1(Pi_S);
        Phi_S = SlopeLimit1(Phi_S);
        psi = SlopeLimit1(psi);
        Pi_psi = SlopeLimit1(Pi_psi);
        Phi_psi = SlopeLimit1(Phi_psi);
    end

    % Increment time                                                                                 
    time = time+dt;                                                                                  
    time_seq = [time_seq, time];
    rhs_g11_seq = [rhs_g11_seq, L2norm(rhs_g11)];
    rhs_Pi11_seq = [rhs_Pi11_seq, L2norm(rhs_Pi11)];
    rhs_Phi11_seq = [rhs_Phi11_seq, L2norm(rhs_Phi11)];

    err_g11_seq = [err_g11_seq, L2norm(g11-g11_exact)];
    err_Pi11_seq = [err_Pi11_seq, L2norm(Pi11-Pi11_exact)];
    err_Phi11_seq = [err_Phi11_seq, L2norm(Phi11-Phi11_exact)];

    C0_seq = [C0_seq, L2norm(C0)];
    C1_seq = [C1_seq, L2norm(C1)];

    Cr00_seq = [Cr00_seq, L2norm(Cr00)];
    Cr01_seq = [Cr01_seq, L2norm(Cr01)];
    Cr11_seq = [Cr11_seq, L2norm(Cr11)];
    for (i=1:size(AH_indicator, 2)-1)
        if (AH_indicator(i) <= 0.0 && AH_indicator(i+1) > 0.0)
            AH_seq = [AH_seq, (x(i+1)*AH_indicator(i)-x(i)*AH_indicator(i+1))/(AH_indicator(i)-AH_indicator(i+1)) - 2.0];
            break;
        end
    end

    if (mod(time, 20) < dt)
       subplot(2,4,1);
       plot(x, g11-g11_exact, 'r'); title(['Error of g11, Pi11, Phi11, t = ', num2str(time)]); hold on;
       plot(x, Pi11-Pi11_exact, 'b'); hold on;
       plot(x, Phi11-Phi11_exact, 'y'); hold off;
       subplot(2,4,2);
       semilogy(time_seq, err_g11_seq, 'r'); title(['L2 norm of error with time']); hold on;
       semilogy(time_seq, err_Pi11_seq, 'b'); hold on; 
       semilogy(time_seq, err_Phi11_seq, 'y'); hold off;
       subplot(2,4,3);
       plot(x, C0, 'r'); title(['C0 and C1, t = ', num2str(time)]); hold on;
       plot(x, C1, 'b'); hold off;
       subplot(2,4,4);
       semilogy(time_seq, C0_seq, 'r'); title(['L2 norm of C0 and C1 with time']); hold on;
       semilogy(time_seq, C1_seq, 'b'); hold off;
       subplot(2,4,5);
       plot(x, g11, 'r'); title(['g11, Pi11, Phi11, t = ', num2str(time)]); hold on;
       plot(x, Pi11, 'b'); hold on;
       plot(x, Phi11, 'y'); hold off;
       %plot(x, Cr00, 'r'); title(['Cr, t= ', num2str(time)]); hold on;
       %plot(x, Cr01, 'b'); hold on;
       %plot(x, Cr11, 'y'); hold off;
       subplot(2,4,6);
       semilogy(time_seq, Cr00_seq, 'r'); title(['L2 norm of Cr with time']); hold on;
       semilogy(time_seq, Cr01_seq, 'b'); hold on;
       semilogy(time_seq, Cr11_seq, 'y'); hold off;

       %plot(x, rhs_g11, 'r'); title(['rhs, t= ', num2str(time)]); hold on;
       %plot(x, rhs_Pi11, 'b'); hold on;
       %plot(x, rhs_Phi11, 'y'); hold off;

       subplot(2,4,7);
       semilogy(time_seq, rhs_g11_seq, 'r'); title(['L2 norm of rhs with time']); hold on;
       semilogy(time_seq, rhs_Pi11_seq, 'b'); hold on;
       semilogy(time_seq, rhs_Phi11_seq, 'y'); hold off;

       subplot(2,4,8);
       plot(x, 0.5*sqrt(2).*(1./power(g11, 0.5) - shift./lapse)); title(['Expansion, t = ', num2str(time)]); drawnow;
       %plot(time_seq, AH_seq, 'r.'); title(['AH position drift with time']); drawnow;

    end
    if (time >= 10 && time - dt <= 10)
        %save(["isoN2mesh", num2str(meshNum),".data"]);
        %fprintf("N=%d\n", N);
        %fprintf("mesh numbers=%d\n", meshNum);

        %fprintf("L2 norm of error of g11 when t=100:  %.2e\n", err_g11_seq(end));
        %fprintf("L2 norm of error of Pi11 when t=100:  %.2e\n", err_Pi11_seq(end));
        %fprintf("L2 norm of error of Phi11 when t=100:  %.2e\n", err_Phi11_seq(end));

        fprintf("Time costs: %.6f\n", toc);
        break;
    end
end
