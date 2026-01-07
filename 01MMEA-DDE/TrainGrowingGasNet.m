function [net,genFlag] = TrainGrowingGasNet(temp1,net,params,Problem,genFlag)
    N = params.N;                        
    MaxIt = params.MaxIt;                
    L = params.L;                       
    epsilon_b = params.epsilon_b;        
    epsilon_n = params.epsilon_n;        
    alpha = params.alpha;               
    delta = params.delta;                
    T = params.T;
    C = net.C;
    w = net.w;
    t = net.t;
    E = net.E;
    nx = net.nx;
    ageSumBefore = net.ageSumBefore;
    flag = net.flag;
    
    for i = 1:size(w,1)
        neighbor = find(C(i,:)==1);
        ageSum(i,:) = sum(t(i,neighbor,:),2); 
        if ageSum(i,:) == ageSumBefore(i,:)
            flag(i,:) = flag(i,:) + 1;
        end
    end

    ageSumBefore = ageSum; 
    maxN = 1.5; 

    gen = ceil(Problem.FE/Problem.N);
    maxgen = ceil(Problem.maxFE/Problem.N);

    if gen <= round(0.9*maxgen)
        maxIter = 1;
        maxPZ = maxN;
        if size(w,1) == round(maxN*N)
            [~,rankFlag] = sort(flag,'descend');
            r = rankFlag(1:round(maxN*N)-N);
            C(r, :) = [];
            C(:, r) = [];
            t(r, :) = [];
            t(:, r) = [];
            w(r, :) = [];
            E(r) = [];
            ageSumBefore(r,:) = [];
            flag(r,:) = [];
            flag = zeros(N,1);
        end
    else
        if size(w,1) < round(maxN*N) && isempty(genFlag)
            maxPZ = maxN;
            maxIter = 1;
        else 
            maxPZ = 1;
            maxIter = 0;
        end
        if size(w,1) == round(maxN*N)
            maxPZ = 1;
            maxIter = 0;
            genFlag = gen;
        end
    end

    if isempty(genFlag)
        for iter = 1:maxIter
            for kk = 1:size(temp1,1)
                nx = nx + 1;
                x = temp1(kk,:);
                d = pdist2(x, w);
                [~, SortOrder] = sort(d);
                s1 = SortOrder(1);
                s2 = SortOrder(2);

                t(s1, :) = t(s1, :) + 1; 
                t(:, s1) = t(:, s1) + 1; 

                E(s1) = E(s1) + d(s1)^2;

                w(s1,:) = w(s1,:) + epsilon_b*(x-w(s1,:)); 
                Ns1 = find(C(s1,:)==1);  
                for j=Ns1                
                    w(j,:) = w(j,:) + epsilon_n*(x-w(j,:));
                end

                C(s1,s2) = 1;           
                C(s2,s1) = 1;          
                t(s1,s2) = 0;            
                t(s2,s1) = 0;         

                C(t>T) = 0;             
                nNeighbor = sum(C);     
                AloneNodes = (nNeighbor==0); 
                if ~isempty(find(AloneNodes == true))  
                    % AloneNodes          
                end
                C(AloneNodes, :) = [];   
                C(:, AloneNodes) = [];   
                t(AloneNodes, :) = [];   
                t(:, AloneNodes) = [];   
                w(AloneNodes, :) = [];   
                E(AloneNodes) = [];     
                ageSumBefore(AloneNodes,:) = []; 
                flag(AloneNodes,:) = []; 

                % Add New Nodes
                if mod(nx, L) == 0 && size(w,1) < round(maxPZ*N) 
                    [~, q] = max(E); 
                    [~, f] = max(C(:,q).*E); 
                    r = size(w,1) + 1; 
                    w(r,:) = (w(q,:) + w(f,:))/2; 
                    C(q,f) = 0;  
                    C(f,q) = 0; 
                    C(q,r) = 1; 
                    C(r,q) = 1; 
                    C(r,f) = 1;  
                    C(f,r) = 1; 
                    t(r,:) = 0; 
                    t(:, r) = 0;
                    E(q) = alpha*E(q);
                    E(f) = alpha*E(f);
                    E(r) = E(q);
                    ageSumBefore(r,:) = 0; 
                    flag(r,:) = 0;
                    
                    if mod(nx, 2*L) == 0 && size(w,1) < round(maxPZ*N)
                        for i = 1:size(w,1) 
                            edgeSum(i) = sum(C(i,:)) + sum(C(:,i));  
                        end
                        [~, q] = min(edgeSum);  
                        D = pdist2(w,w,'cityblock');  
                        D(logical(eye(length(D)))) = inf;  
                        D = D(q,:);         
                        D(C(q,:)==1) = inf;  
                        D(C(:,q)==1) = inf;  
                        [~, f] = min(D);     
                        r = size(w,1) + 1;   
                        w(r,:) = (w(q,:) + w(f,:))/2; 
                        C(q,f) = 0;          
                        C(f,q) = 0;          
                        C(q,r) = 1;          
                        C(r,q) = 1;         
                        C(r,f) = 1;         
                        C(f,r) = 1;          
                        t(r,:) = 0;          
                        t(:, r) = 0;        
                        E(q) = alpha*E(q);   
                        E(f) = alpha*E(f);   
                        E(r) = E(q);        
                        ageSumBefore(r,:) = 0;  
                        flag(r,:) = 0;       
                    end
                end

                % Decrease Errors
                E = delta*E;         
            end
        end
        net.w = w;                   
        net.E = E;                   
        net.C = C;                  
        net.t = t;                 
        net.nx = nx;                
        net.ageSumBefore = ageSumBefore; 
        net.flag = flag;        
    end
end