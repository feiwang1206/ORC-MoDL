function [ B ] = RebuildField7(bc,rhoReference,x,y,z,mode )
% Calculate the fields on a grid
% define by the vector x,y and z
numberDregree = max(size(bc(1).coefficient,1))-1;
numberOrder = max(size(bc(1).coefficient,2))-1;

%activate the paralle function
matlabVersion = version;
matlabVersion = str2num(matlabVersion(1:3));
if matlabVersion < 8.2
    [TF,~] = license('checkout', 'Distrib_Computing_Toolbox');
    if TF
        schd = findResource('scheduler', 'configuration', 'local');
        numWorkers = schd.ClusterSize;
    end

    if matlabpool('size') == 0  && TF && numWorkers >1
        % checking to see if the pool is already open and of we have the licence
        % and at least 2 cores
        matlabpool open
    end
elseif matlabVersion >= 8.2 
    poolobj = gcp('nocreate'); % If no pool, do not create new one.
	if isempty(poolobj)
% 		parpool;
    end
end

B = zeros(3,size(x,2),size(y,2),size(z,2));
Bx = zeros(size(x,2),size(y,2),size(z,2));
By = zeros(size(x,2),size(y,2),size(z,2));
Bz = zeros(size(x,2),size(y,2),size(z,2));
% parfor i=1:size(x,2)
for i=1:size(x,2)
    Tempx = zeros(size(y,2),size(z,2));
    Tempy = zeros(size(y,2),size(z,2));
    Tempz = zeros(size(y,2),size(z,2));
    for j=1:size(y,2)
        for k=1:size(z,2)
            rho = sqrt(x(i)^2+y(j)^2+z(k)^2);
            theta = acos(z(k)/rho);
            %To correct the phase stuff
            phi = atan2(y(j),x(i));
%             if rho<rhoReference
                for l=0:numberDregree
                    scaling=(rho/rhoReference)^(l);
                    L = legendre(l,cos(theta));
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
                                Yc = K*L(m+1)*cos(m*phi);
                                Ys = 0; % This one is just impossible (i.e. m<0 and m=0)
                            end
                            Tempx(j,k) = Tempx(j,k) + bc(1).coefficient(l+1,m+l+1)*scaling*Yc;

                    end
                end

%             end
        end
    end
    Bx(i,:,:) = Tempx;
end
B = Bx;

end

