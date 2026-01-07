function [Population,CrowdDis] = EnvironmentalSelection(Population,N)
    n = length(Population);
    dist = pdist2(Population.decs,Population.decs);
    V = 0.2 * prod(max(Population.decs) - min(Population.decs)) .^ (1./size(Population.decs,2));
    DominationX = zeros(n);
    for i = 1 : n
        for j = i + 1 : n
            if dist(i,j) > V
                continue
            end
            L1 = Population(i).objs < Population(j).objs;
            L2 = Population(i).objs > Population(j).objs;
            if all(L1 | (~L2))
                DominationX(i,j) = 0;
                DominationX(j,i) = 1;
            elseif all(L2 | (~L1))
                DominationX(i,j) = 1;
                DominationX(j,i) = 0;
            end
        end
    end
    LocalC = zeros(1,n);
    for i = 1 : n
        tmp = dist(i,:);
        index = tmp < V;
        LocalC(i) = (sum(DominationX(i,index))) ./ sum(index);
    end
    unique_LocalC = unique(LocalC);
    selected = false(1,n);
    selected_count = 0;

    for val = unique_LocalC
        idx = find(LocalC == val);
        remaining = N - selected_count;
        if length(idx) + selected_count <= N
            selected(idx) = true;
            selected_count = selected_count + length(idx);
        else
            temp_pop = Population(idx);
            K = length(idx) - remaining;
            Del = Truncation(temp_pop.decs,temp_pop.objs, K);
            selected(idx(~Del)) = true;
            selected_count = selected_count + sum(~Del);
            if selected_count >= N
                break;
            end
        end
    end

    Population = Population(selected);
    CrowdDis = Crowding(Population.decs);

end


function Del = Truncation(PopObj, PopDec, K)
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
