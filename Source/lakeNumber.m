function Ln = lakeNumber(bthA,bthD,uSt,St,metaT,metaB,rhoHyp)
%----Author: Jordan S Read 2009 ----
% updated on 2 April 2010
% updated on 25 January 2012

%Calculates the lake number of a system using the
%following equation:
%
%   Ln = (g*St*(1-(ht/hm)))/(p0*(uStar^2)*(A0^1.5)*(1-hv/hm)).

%
%References:
%   -Imberger, Jorg, and John C. Patterson. "Physical Limnology." Advances in 
%    Applied Mechanics 27 (1990): 314-317.

g = 9.81;
dz = 0.1;



% if bathymetry has negative values, intepolate to 0
if lt(min(bthD),0)
    useI = ge(bthD,0);
    if ~eq(bthD,0)
        depT = [0 bthD(useI)];
    else
        depT = bthD(useI);
    end
    bthA = interp1(bthD,bthA,depT);
    bthD = depT;
end

[Zo,Io] = min(bthD);
Ao = bthA(Io);
if eq(Ao,0)
    error('Surface area cannot be zero, check *.bth file')
end

%interpolates the bathymetry data
layerD = Zo:dz:max(bthD);
layerA  = interp1(bthD,bthA,layerD);


%find depth to the center of volume
Zv = layerD.*layerA*dz;                    
Zcv = sum(Zv)/sum(layerA)/dz;   %should only need to do this once per  
                                %lake analyzer run...move out.
St_uC = St*Ao/g;
% Calculates the Lake Number according to the formula provided
Ln = g*St_uC*(metaT+metaB)/(2*rhoHyp*uSt^2*Ao^(3/2)*Zcv);

end
