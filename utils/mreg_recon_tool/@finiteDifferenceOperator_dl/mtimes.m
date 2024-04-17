function Q = mtimes(A,B)

if strcmp(class(A),'finiteDifferenceOperator_dl')    
    if A.adjoint==1
        if isvector(B)
            Q = zeros(size(B));
            Q(end+1) = 0;
            Q(1) = -B(1);
            Q(2:end-1) = B(1:end-1)-B(2:end);
            Q(end) = B(end);
            
        else
            s = [1:length(size(B))];
            s(A.direction) = [];
            s = [A.direction, s];
            B = permute(B, s);

            ddots = '';
            for k=1:length(size(B))-1
                ddots = [',:', ddots];
            end

            E1 = eval(['B(1', ddots, ');']);
            Eend = eval(['B(end', ddots, ');']);

%             Q = cat(1, -E1, -diff(B,1,1), Eend);
            Q = cat(1, -E1, B(1:end-1,:,:)-B(2:end,:,:), Eend);
            Q = ipermute(Q, s);
        end

    else
        if isvector(B)
            Q = B(2:end)-B(1:end-1);
        else
            if A.direction==1
                Q=B(2:end,:,:)-B(1:end-1,:,:);
            elseif A.direction==2
                Q=B(:,2:end,:)-B(:,1:end-1,:);
            elseif A.direction==3
                Q=B(:,:,2:end)-B(:,:,1:end-1);
            end
        end
        
    end

% now B is the operator and A is the vector
else
    Q = mtimes(B',A')';
    
end