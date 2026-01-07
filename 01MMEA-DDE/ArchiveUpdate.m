function [Population,CrowdDis,net,genFlag,flag_Ini_Arc] = ArchiveUpdate(Population,N,eps,st,Problem,params,net,genFlag,flag_Ini_Arc)
    n = length(Population);
    if eps~=1 && st<0.5
       eps = 2*(1-eps)/(2*st+1)+2*eps-1;
    end
    %% Select global Pareto front
    [FrontNo,MaxFNo] = NDSort(Population.objs,n);
    next = FrontNo==1;
    first_pf = Population(next);
    new_pop = first_pf;
    remain_pop = Population(~next);

    V=0.2*prod(max(Population.decs)-min(Population.decs)).^(1./size(Population.decs,2));

    while ~isempty(remain_pop)
        %% Delete close solutions
        Dec_dist = min(pdist2(new_pop.decs,remain_pop.decs));
        index = Dec_dist<V;
        remain_pop(index) = [];
        if isempty(remain_pop)
            break;
        end

        %% Select remaining solutions
        [FrontNo,MaxFNo] = NDSort(remain_pop.objs,length(remain_pop));
        pick_pop = remain_pop(FrontNo==1);
        [nF,~] = NDSort([pick_pop.objs .* (1-eps); first_pf.objs],length(pick_pop)+length(first_pf));
        nF = nF(1:length(pick_pop));

        maxnF = max(nF);

        if maxnF>1
            new_pop = [new_pop pick_pop(nF==1)];
            remain_pop = remain_pop(FrontNo~=1);
            break;
        else
            new_pop = [new_pop pick_pop];
            remain_pop = remain_pop(FrontNo~=1);
        end
    end
    Population = new_pop;
    
    %% Truncation
    if length(Population) > N
         Del = Truncation(Population.objs, Population.decs, length(Population) - N );
         Population = Population(~Del);
    end
    CrowdDis = Crowding(Population.decs);
   

    if flag_Ini_Arc &&length(Population) > 2 && isempty(find(isnan(Population.decs)==true)) && st<1 && isempty(genFlag)
        [net,genFlag] = TrainGrowingGasNet(Population.decs,net,params,Problem,genFlag);
    end
    flag_Ini_Arc=1;
end

function Del = Truncation(PopObj, PopDec, K)
    %% Truncation
    Distance_Pop = pdist2(PopObj, PopObj);
    Distance_Dec = pdist2(PopDec, PopDec);
    
    Distance_Pop_norm = (Distance_Pop - min(Distance_Pop(:))) / (max(Distance_Pop(:)) - min(Distance_Pop(:)) + eps);
    Distance_Pop_norm(logical(eye(length(Distance_Pop_norm)))) = inf;
    Distance_Dec_norm = (Distance_Dec - min(Distance_Dec(:))) / (max(Distance_Dec(:)) - min(Distance_Dec(:)) + eps);
    Distance_Dec_norm(logical(eye(length(Distance_Dec_norm)))) = inf;
    
    D =  Distance_Pop_norm + Distance_Dec_norm;
    
    Del = false(1, size(PopObj, 1));
    while sum(Del) < K
        Remain = find(~Del);
        Temp = sort(D(Remain, Remain), 2);
        [~, Rank] = sortrows(Temp);
        Del(Remain(Rank(1))) = true;
    end
end

function [CrowdDis]=Crowding(Pop)
    [N, ~]=size(Pop);
    K=N-1;
    Z = min(Pop,[],1);
    Zmax = max(Pop,[],1);
    pop=(Pop-repmat(Z,N,1))./repmat(Zmax-Z,N,1);
    distance=pdist2(pop,pop);
    [value,~]=sort(distance,2);
    CrowdDis=K./sum(1./value(:,2:N),2);
end