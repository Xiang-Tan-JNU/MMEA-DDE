function [Population,Fitness] = EnvironmentalSelectionSup(Population,N,net)
    V = net.w;
    C = net.C;
    Distance = pdist2(V,V);
    Distance = Distance.*~C;
    d = zeros(1,length(V));
    for i = 1:length(V)
        Dis = Distance(i,:);
        Dis(Dis==0) = [];
        d(i) = min(Dis);
    end
    theta = max(d);
    Fitness = CalFitnessSup(Population.decs,V);

    Next = Fitness < theta;
    if sum(Next) < N
        [~,Rank] = sort(Fitness);
        Next(Rank(1:N)) = true;
    elseif sum(Next) > N
        Del  = Truncation(Population(Next).decs,sum(Next)-N);
        Temp = find(Next);
        Next(Temp(Del)) = false;
    end
    Population = Population(Next);
    Fitness    = Fitness(Next);
    [Fitness,rank] = sort(Fitness);
    Population = Population(rank);
end

function Del = Truncation(PopDec,K)
    D = pdist2(PopDec,PopDec);
    D(logical(eye(length(D)))) = inf;
    Del = false(1,size(PopDec,1));
    while sum(Del) < K
        Remain   = find(~Del);
        Temp     = sort(D(Remain,Remain),2);
        [~,Rank] = sortrows(Temp);
        Del(Remain(Rank(1))) = true;
    end
end

function D_Dec = CalFitnessSup(PopDec,V)   
    Distance = pdist2(PopDec,V);
    Distance = sort(Distance,2);
    D_Dec = Distance(:,1);
end