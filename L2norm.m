function [norm] = L2norm(u);
    Globals1D;
    norm = sqrt(sum(sum((invV*u).*(invV*u)).*abs(VX(2:end)-VX(1:end-1)))/2);

