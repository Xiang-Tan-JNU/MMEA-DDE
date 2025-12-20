function D_Dec = CalFitnessSup(PopDec,V)
    %% Calculate D(i)    
    Distance = pdist2(PopDec,V);
    Distance = sort(Distance,2);
    D_Dec = Distance(:,1);
end