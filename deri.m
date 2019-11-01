function [pPlus, pMinus] = deri(u, u_left, u_right)

Globals1D;

duP = zeros(Nfp*Nfaces,K);
duM = zeros(Nfp*Nfaces,K);
duM(:) = (u(vmapM)-u(vmapP)).*(nx(:)-abs(nx(:)))/2;
duP(:) = (u(vmapM)-u(vmapP)).*(nx(:)+abs(nx(:)))/2;

duM(mapI) = (u(vmapI)-u_left(vmapI)).*(nx(mapI)-abs(nx(mapI)))/2;
duM(mapO) = 0;

duP(mapI) = 0;
duP(mapO) = (u(vmapO)-u_right(vmapO)).*(nx(mapO)+abs(nx(mapO)))/2;

pPlus = rx.*(Dr*u) - LIFT*(Fscale.*(duP));
pMinus = rx.*(Dr*u) - LIFT*(Fscale.*(duM));
return

