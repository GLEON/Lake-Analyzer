function [ uStar ] = uStar( wndSpeed,wndHeight,averageEpiDense )
%----Author: Jordan S Read 2009 ----
% **** updated 19 Feb 2010 ****

rhoAir = 1.2;
vonK = 0.4; %von Karman k

%equation (1) provided by:
%같 Hicks, B.B., 1972. A procedure for the formulation of bulk transfer
%같 coefficients over water bodies of different sizes. Boundary-Layer
%같 Meterology 3: 201-213

%equation (2) provided by:
%같 Amorocho, J., DeVries, J.J., 1980. A new evaluation of the wind 같
%같 stress coefficient over water surfaces. Journal of Geophysical  같
%같 Research 85: 433-442.

%equation (3) provided by:
%같 Fischer, H.B., List, E.J., Koh, R.C.Y., Imberger, J., Brooks, N.H.,
%같 1979. Mixing in inland and coastal waters. Academic Press.

%equation (4) provided by:
%같 Imberger, J., 1985. The diurnal mixed layer. Limnology and Oceanography
%같 30: 737-770.

% ~ eqn (1) ~
if wndSpeed < 5
    CD = 0.001;
else
    CD = 0.0015;
end

% ~ eqn (2) ~
%If the windspeed is not measured at a height of 10 meters, corrected by
%(eqn 21) from Amorocho and DeVries, 1980
if wndHeight ~= 10
    wndSpeed = wndSpeed/(1-sqrt(CD)/vonK*log(10/wndHeight));
end


% ~ eqn (3) ~
tau = CD*rhoAir*wndSpeed^2;

% ~ eqn (4) ~
uStar = sqrt(tau/averageEpiDense);

end

