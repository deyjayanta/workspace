function [x y] = random_data(data,res)
    dt = [data res];
    s = size(dt);
    dt = datasample(dt,s(1),'Replace',false);
    
    x = dt(:,1:s(2)-1);
    y = dt(:,s(2));
end