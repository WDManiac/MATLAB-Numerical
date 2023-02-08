function [V,C] = LloydOptimizePeri(Ms, tar, err, L)
% This function optimize an existing Voronoi Tesselation to a more uniform
% state, by which means, decrease its polydispersity(CA).
% Ms --- Initial Centroids
% tar --- Target CA
% err --- Error Margin
% L --- Size of Periodic Box

%% First get initial CA
Len = length(Ms)/9;
CArec = [];

[V,C] = voronoin(Ms); % Generate Voronoi Tesselation for given points

% Plot and show the initial Lattice
figure,
for j = 4*Len+1 : 5*Len
    hold on, patch([V(C{j},1)', V(C{j}(1),1)],[V(C{j},2)', V(C{j}(1),2)],'y')
end
axis equal

% Calculate Volume of Initial Lattice to get CA
for j = 1:Len
    Siz(j) = length(C{4*Len + j});
    % For each cell, go through its vertices and compute area with
    % cross product
    vol = 0;
    for ss = 1:Siz(j)-1
        vol = vol + V(C{4*Len + j}(ss),1)*V(C{4*Len + j}(ss+1),2) - V(C{4*Len + j}(ss),2)*V(C{4*Len + j}(ss+1),1);
    end
    vol = vol + V(C{4*Len + j}(Siz(j)),1)*V(C{4*Len + j}(1),2) - V(C{4*Len + j}(Siz(j)),2)*V(C{4*Len + j}(1),1);
    
    % Check if the result is possitive or negative
    if vol>0
        Vol(j) = vol/2;
        Sign(j) = 1;
    else
        Vol(j) = -vol/2;
        Sign(j) = -1;
    end
end
% Calculate CA, save picture
[para,~] = gamfit(Vol);
CA = sqrt(para(1) * para(2)^2);
disp(['Polydispersity: ',num2str(CA)])
hold on, title(['Iter Step: 1, c_{A} = ',num2str(CA)])
set(gca,'fontsize',14)
% saveas(gcf,['LloydIterStep1'],'jpg')

% Record CA
CArec = [CArec, CA];


%% Then start Lloyd Iteration
% The idea is to use centroids of previous lattice as seeds for new
% lattice!
if CA > tar
    i = 1;
    while abs(CA - tar) >= err && CA > tar
        
        i = i+1; % Update Iteration Step No.
        
        % Get New Centroids Locations
        XM = []; YM = [];
        for j = 1:length(C)
            Vx{j} = V(C{j},1);
            Vy{j} = V(C{j},2);
            
            x(j) =  mean(Vx{j});
            y(j) =  mean(Vy{j});
            
            if x(j)>0 && x(j)<L && y(j)>0 && y(j)<L
                XM = [XM,x(j)];
                YM = [YM,y(j)];
            end
        end
        
        % Reorganize
        Xcord = []; Ycord = [];
        Xcord = [Xcord, XM - L]; Ycord = [Ycord, YM + L];
        Xcord = [Xcord, XM]; Ycord = [Ycord, YM + L];
        Xcord = [Xcord, XM + L]; Ycord = [Ycord, YM + L];
        Xcord = [Xcord, XM - L]; Ycord = [Ycord, YM];
        Xcord = [Xcord, XM]; Ycord = [Ycord, YM];
        Xcord = [Xcord, XM + L]; Ycord = [Ycord, YM];
        Xcord = [Xcord, XM - L]; Ycord = [Ycord, YM - L];
        Xcord = [Xcord, XM]; Ycord = [Ycord, YM - L];
        Xcord = [Xcord, XM + L]; Ycord = [Ycord, YM - L];
        
        % Generate New Lattice!
        Ms = [Xcord; Ycord]';
        [V,C] = voronoin(Ms);
        % Show New Lattice!
        figure,
        for j = 4*Len+1 : 5*Len
            hold on, patch([V(C{j},1)', V(C{j}(1),1)],[V(C{j},2)', V(C{j}(1),2)],'y')
        end
        axis equal
        
        % For New Lattice, Compute CA for iteration
        Vol = [];
        Sign = [];
        for j = 1:Len
            Siz(j) = length(C{4*Len + j});
            % For each cell, go through its vertices and compute area with
            % cross product
            vol = 0;
            for ss = 1:Siz(j)-1
                vol = vol + V(C{4*Len + j}(ss),1)*V(C{4*Len + j}(ss+1),2) - V(C{4*Len + j}(ss),2)*V(C{4*Len + j}(ss+1),1);
            end
            vol = vol + V(C{4*Len + j}(Siz(j)),1)*V(C{4*Len + j}(1),2) - V(C{4*Len + j}(Siz(j)),2)*V(C{4*Len + j}(1),1);
            
            % Check if the result is possitive or negative
            if vol>0
                Vol(j) = vol/2;
                Sign(j) = 1;
            else
                Vol(j) = -vol/2;
                Sign(j) = -1;
            end
        end
        
        % Calculate CA and Save Figure
        [para,~] = gamfit(Vol);
        CA = sqrt(para(1) * para(2)^2);
        disp(['Polydispersity: ',num2str(CA)])
        hold on, title(['Iter Step: ',num2str(i),', c_{A} = ',num2str(CA)])
        set(gca,'fontsize',14)
        saveas(gcf,['LloydIterStep',num2str(i)],'jpg')
        
        CArec = [CArec, CA]; % Record CA
    end
    
    % Record the whole process
    figure,
    plot(1:i, CArec,'ro-','linewidth',1.5)
    set(gca,'fontsize',14)
    title(['Lloyd Algorithm Optimized c_{A}, with ',num2str(Len),' cells'])
    grid on,
    ylabel 'c_{A}'
    xlabel 'Iter. Step #'
    saveas(gcf,'OptCurve','jpg')
    saveas(gcf,'OptCurve','fig')
    
    if abs(CA - tar) >= err
        disp('CanT Converge!')
    end
else
    disp('Initial Polydispersity is smaller than target already, please specify another value!')
end
