function [rho] = waterDensity(T,S)
%% WDENSITY Density of water with supplied Temp and Salinity
% rho = WDENSITY(Temperature, Salinity) is the density of water with the
% given temperature and salinity
% T is in °C
% S is in Pratical Salinity Scale units (dimensionless)
% T and S should be the same size, unless S is a scalar
% rho is in grams/Liter
%
% <<<--- Effective range: 0-40°C, 0.5-43 Salinity --->>>
%
%----Author: Jordan S Read 2011 ----
% Editied by:
% 1/30/2010 - Luke Winslow <lawinslow@gmail.com>

% three options: use all UNESCO, use all Martin & McCutcheon, or element
% wise for both

MM = false; UN = false; elm = false;
if eq(nargin,1);
    S = 0;
end


Trng = [0 40];
Srng = [0.5 43];

% check to see if all 
if all(lt(S,Srng(1)))
    MM = true;      % use Martin & McCutcheon
elseif ~(sum(T<Trng(1))||sum(T>Trng(2))) && ~(sum(S<Srng(1))||sum(S>Srng(2)))
    UN = true;      % use UNESCO
else
    elm = true;     % element-wise choice between M&M and UNESCO
end


%% use methods:
if MM
    % <<equation provided by:
    %°° Martin, J.L., McCutcheon, S.C., 1999. Hydrodynamics and Transport °°
    %°° for Water Quality Modeling. Lewis Publications, Boca              °°
    %°° Raton, FL, 794pp.  >>
    if ~eq(length(T),1)
        dens = @(T)(1000*(1-(T+288.9414)*(T-3.9863)^2/(508929.2*(T+68.12963))));
        rho  = arrayfun(dens,T);
    else
        rho = (1000*(1-(T+288.9414)*(T-3.9863)^2/(508929.2*(T+68.12963))));
    end
elseif UN
    % <<equations provided by:
    %°° Millero, F.J., Poisson, A., 1981. International one-atmosphere    °°
    %°° equation of state of seawater. UNESCO Technical Papers in Marine  °°
    %°° Science. No. 36     >>
    % --eqn (1):
    rho_0 = 999.842594+6.793952e-2*T-9.095290e-3*T.^2+...
        1.001685e-4*T.^3-1.120083e-6*T.^4+6.536335e-9*T.^5;
    % --eqn (2):
    A = 8.24493e-1-4.0899e-3*T+7.6438e-5*T.^2-8.2467e-7*T.^3+5.3875e-9*T.^4;
    % --eqn (3):
    B = -5.72466e-3+1.0227e-4*T-1.6546e-6*T.^2;
    % --eqn (4):
    C = 4.8314e-4;
    % --eqn (5):
    rho = rho_0+A.*S+B.*S.^(3/2)+C.*S.^2;
else
    rho = T*NaN;
    for j = 1:length(T)
        rho(j) = waterDensity(T(j),S(j));
    end
end

end
