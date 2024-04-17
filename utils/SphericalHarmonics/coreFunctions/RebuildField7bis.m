function [ B ] = RebuildField7bis(bc,bs,rhoReference,rk,mode )
% Calculate the fields at given position
% the list rk gives each position
% with the spherical harmonic equations
numberDregree = max(size(bc(1).coefficient,1),size(bs(1).coefficient,1))-1;
numberOrder = max(size(bc(1).coefficient,2),size(bs(1).coefficient,2))-1;

if size(rk,2) <6
    for i=1:size(rk,1)
        x = rk(i,1);
        y = rk(i,2);
        z = rk(i,3);
        rk(i,4) = sqrt(x^2+y^2+z^2); % rho
        rk(i,6) = acos(z/rk(i,4)); % theta
        rk(i,5) = atan2(y,x); % phi
    end
end

B = zeros(3,size(rk,1));
for i=1:size(rk,1)
    rho = rk(i,4);
    phi = rk(i,5);
    theta = rk(i,6);
        for l=0:numberDregree
            scaling=(rho/rhoReference)^(l);
            L = legendre(l,cos(theta));
            %calculate the minimum between l and maxM
            maxOrder = min(l,numberOrder);
            for m=-l:l
                if strcmp(mode,'norm')
                    K = sqrt((2*l+1)*factorial(l-abs(m))/(factorial(l+abs(m))*4*pi)); %mathematical normalization
                elseif strcmp(mode,'sch')
                    K = sqrt(factorial(l-abs(m))/(factorial(l+abs(m)))); % schmidt quasi normal without the factor (2*l+1)/(4*pi)
                else
                    disp('error, the normalization is not recognized')
                    K = 0;
                end
                if m>0
                    Yc = sqrt(2)*K*L(abs(m)+1)*cos(abs(m)*phi);
                    Ys = sqrt(2)*K*L(abs(m)+1)*sin(abs(m)*phi);
                elseif m<0
                    Yc = sqrt(2)*K*L(abs(m)+1)*sin(abs(m)*phi);
                    Ys = sqrt(2)*K*L(abs(m)+1)*cos(abs(m)*phi);
                else
                    Yc = K*L(-m+1)*cos(m*phi);
                    Ys = 0; % This one is just impossible (i.e. m<0 and m=0)
                end
                B(1,i) = B(1,i) + bc(1).coefficient(l+1,m+l+1)*scaling*Yc;
                B(2,i) = B(2,i) + bc(2).coefficient(l+1,m+l+1)*scaling*Yc;
                B(3,i) = B(3,i) + bc(3).coefficient(l+1,m+l+1)*scaling*Yc;

                B(1,i) = B(1,i) + bs(1).coefficient(l+1,m+l+1)*scaling*Ys;
                B(2,i) = B(2,i) + bs(2).coefficient(l+1,m+l+1)*scaling*Ys;
                B(3,i) = B(3,i) + bs(3).coefficient(l+1,m+l+1)*scaling*Ys;
            end
        end
end
%l order; m degree; B is superpositon of the l order and m degree SH,
%numberOrder is the max order, so do numberDegree.
%together. B(1,:) is Bx field, B(2,:) is By field, and B(3,:) is Bz field.


