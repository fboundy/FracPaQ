function Limits = getRotatedMapLimits(traces, angle)
    if angle > 45
        angle = angle - 90;
    end;
    if angle < -45
        angle = angle + 90;
    end;
    theta = deg2rad(angle);    
    NodesXY = [];

%    f=figure;
    for i = 1:length(traces)
        for j = 1:traces(i).nNodes
            NodesXY = [NodesXY [traces(i).Node(j).x; traces(i).Node(j).y]];
        end 
    end 
    
    xMin = min(NodesXY(1,:)) ;
    xMax = max(NodesXY(1,:)) ;
    yMin = min(NodesXY(2,:)) ;
    yMax = max(NodesXY(2,:)) ;
    
    
    disp(['Xmin: ',num2str(xMin,'%.0f')]);
    disp(size(NodesXY,2));
    
    xMid = (xMin + xMax) / 2;
    yMid = (yMin + yMax) / 2;
    
    ShiftXY = [repelem(xMid, size(NodesXY,2));repelem(yMid, size(NodesXY,2))];
    S = NodesXY - ShiftXY;
    
    R = [cos(theta) -sin(theta) ; sin(theta) cos(theta)];
    
    NodesIJ = R*S;
    
%    plot(S(1,:),S(2,:),'g',NodesIJ(1,:), NodesIJ(2,:),'r');
    iMin = min(NodesIJ(1,:));
    iMax = max(NodesIJ(1,:));
    jMin = min(NodesIJ(2,:));
    jMax = max(NodesIJ(2,:));
    
    v00 = [iMin; jMin];
    v01 = [iMin; jMax];
    v10 = [iMax; jMin];
    v11 = [iMax; jMax];
    V = [v00 v01 v10 v11];  
%    plot(S(1,:),S(2,:),'g',NodesIJ(1,:), NodesIJ(2,:),'r',V(1,:),V(2,:),'b');

    
    m = [xMid; yMid];
    M = [m m m m];       

    R = [cos(-theta) -sin(-theta) ; sin(-theta) cos(-theta)];

    Z = R * V + M;

    disp(num2str(Z,'%12.0f'));
    
    
%    f=figure;
    XAx = [Z(1,1) Z(1,3);  Z(2,1), Z(2,3)];
    YAx = [Z(1,1) Z(1,2);  Z(2,1), Z(2,2)];
%    plot(NodesXY(1,:),NodesXY(2,:),'g',XAx(1,:), XAx(2,:),'r',YAx(1,:), YAx(2,:),'b');

    Limits = Z;
end
