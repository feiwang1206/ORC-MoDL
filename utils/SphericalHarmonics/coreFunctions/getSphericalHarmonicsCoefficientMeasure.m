function [bc] = getSphericalHarmonicsCoefficientMeasure(B,degreeMax,orderMax,rk,mode)


% Same as getSphericalHarmonicsCoefficientMeasure4 but with another
% normalization

maxSphericalHarmonicsOrder = orderMax;
% Phi = rectangle method - 2m max+1 points - order
% Theta = Gauss method - lmax+1 points - degree
nbrPointsPhi = 2*orderMax+1;
nbrPointsTheta = degreeMax+1;

[gaussLegendreAbscissas , gaussLegendreWeight] = gauss(nbrPointsTheta);
% Phi - Alex recommand to use the left rectangle technic, to match
% with his code
startPhi = 0;
stopPhi = 2*pi;
stepPhi = (stopPhi-startPhi)/nbrPointsPhi;
AbsissasPhi = startPhi:stepPhi:stopPhi-stepPhi;
WeightPhi = AbsissasPhi(2)-AbsissasPhi(1);

% theta
AbsissasTheta = gaussLegendreAbscissas;
WeightTheta = gaussLegendreWeight;


for l=0:degreeMax
    
    maxOrder = min(l,orderMax);
    for m=-l:l
        sumBx_bc = 0;
        sumBy_bc = 0;
        sumBz_bc = 0;
        sumBx_bs = 0;
        sumBy_bs = 0;
        sumBz_bs = 0;
        point=1;

        %scaling = (2*l+1)/(4*pi);
        for i=1:nbrPointsPhi % Phi index
%             if rk(point,5) ~= AbsissasPhi(i)
%                 disp('error, Phi for the field calculation and for the integration are not the same')
%             end
            for j=1:nbrPointsTheta % Theta index
%                 if rk(point,6) ~= acos(AbsissasTheta(j))
%                    disp('error, Theta for the field calculation and for the Gauss Quadrature are not the same')
%                 end
                if strcmp(mode,'norm')
                    K = sqrt((2*l+1)*factorial(l-abs(m))/(factorial(l+abs(m))*4*pi)); %mathematical normalization
                elseif strcmp(mode,'sch')
                    K = (2*l+1)/(4*pi)*sqrt(factorial(l-abs(m))/(factorial(l+abs(m)))); % schmidt quasi normal
                else
                    disp('error, the normalization is not recognized')
                    K = 0;
                end
                S = legendre(l,AbsissasTheta(j)); %and not cos(theta), thank to super integral trick from Alex
                
                if m>0
                    Yc = sqrt(2)*K*S(abs(m)+1)*cos(abs(m)*AbsissasPhi(i));
                    Ys = sqrt(2)*K*S(abs(m)+1)*sin(abs(m)*AbsissasPhi(i));
                elseif m<0
                    Yc = sqrt(2)*K*S(abs(m)+1)*sin(abs(m)*AbsissasPhi(i));
                    Ys = sqrt(2)*K*S(abs(m)+1)*cos(abs(m)*AbsissasPhi(i));
                else
                    Yc = K*S(m+1)*cos(m*AbsissasPhi(i));
                    Ys = 0; % This one is just impossible (i.e. m<0 and m=0)
                end
                sumBx_bc = sumBx_bc + WeightPhi*WeightTheta(j)*B(point)*Yc;
                
                point = point+1;
            end
        end

        bc(1).coefficient(l+1,m+l+1) = sumBx_bc;
    end
end

