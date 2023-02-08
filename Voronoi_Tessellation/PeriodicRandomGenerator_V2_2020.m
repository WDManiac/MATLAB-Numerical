%% Automatic Script for Surface Evolver Periodic Hexagon Lattice Generation (V1.2)
% Author: Ran Li
% Using "voronoin" function to generate Periodic Voronoi Tesselation to represent
% cell aggregate.

% This code is modified to generate square Periodic Cellular Lattice with CVT method 

clear all
clc
close all

cp = 0.045; % Quenched Disorder!
gamma = 0.22; % adhesion term
tar = 0.24; % Target CA
err = 0.005; % CA error Margin

%% Generate Centroids (Poisson Process)
% 1> Parent Process
% lambda = 100;                       % density
% M = 0;
% U = unifrnd(0,1);

% while U >= exp(-lambda)           % condition
%    U = U*unifrnd(0,1);
%    M=M+1;
% end 
% 
% if M < 1
%    M = 1;
% end

Num = 33;
M = Num^2;

a = 0;b = Num;      %Take [0,100]*[0,100]*[0,100] as distribution region
c = 0;d = Num;
A = zeros(1,M);
B = zeros(1,M);
for i = 1:M 
    U1 = unifrnd(0,1);
    A(i) = (b-a)*U1;
    U2 = unifrnd(0,1);
    B(i) = (d-c)*U2;
end

X = []; Y = [];

% 2> Child Process for j=1:M
% for j=1:M
%     n = Num*Num/M;
%     r = 10;
%     u1 = zeros(1,n);
%     u2 = zeros(1,n);
%     R = zeros(1,n);
%     x = zeros(1,n);
%     y = zeros(1,n);
%     theta = zeros(1,n);
% 
%     for i = 1:n
%         u1(i) = unifrnd(0,1);
%     end
% 
%     R = r * sqrt(u1);
%     R = sort(R);
% 
%     for i = 1:n
%         u2(i) = unifrnd(0,1);
%     end
% 
%     theta = 2*pi*u2;
% 
%     for i = 1:n
%         x(i) = A(j) + R(i) * cos(theta(i));
%         y(i) = B(j) + R(i) * sin(theta(i));
%         
%         if x(i)<0 
%             x(i) = x(i) + Num;
%         elseif x(i)>Num
%             x(i) = x(i) - Num;
%         end
%         
%         if y(i)<0 
%             y(i) = y(i) + Num;
%         elseif y(i)>Num
%             y(i) = y(i) - Num;
%         end
%         
% %         if x(i)>0 && x(i)<Num && y(i)>0 && y(i)<Num
%             X = [X,x(i)];
%             Y = [Y,y(i)];
% %         end
%     end
%     
% end

X = A; Y = B;

L = length(X);

Xcord = []; Ycord = [];
Xcord = [Xcord, X - Num]; Ycord = [Ycord, Y + Num];
Xcord = [Xcord, X]; Ycord = [Ycord, Y + Num];
Xcord = [Xcord, X + Num]; Ycord = [Ycord, Y + Num];
Xcord = [Xcord, X - Num]; Ycord = [Ycord, Y];
Xcord = [Xcord, X]; Ycord = [Ycord, Y];
Xcord = [Xcord, X + Num]; Ycord = [Ycord, Y];
Xcord = [Xcord, X - Num]; Ycord = [Ycord, Y - Num];
Xcord = [Xcord, X]; Ycord = [Ycord, Y - Num];
Xcord = [Xcord, X + Num]; Ycord = [Ycord, Y - Num];

Ms = [Xcord; Ycord]';

%% Generate Centroids (Ordered Rectanguler)
% Num=20;% Size of periodic box
% Cord = meshgrid(0:Num);
% lambda = 0.9; % Coefficient for randomness:
% %     multiplying normally distributed randomness
% % range = 9; % Radius of the selected region
% % rad = 15;
% 
% B = 0.06; % ratio of margin
% 
% % % Tension Info
% % beta = 1;
% % gamma1 = 1;
% % gamma2 = 1.67;
% 
% Xm=Cord;
% Ym=Cord';
% 
% for i=1:Num+1
%     if mod(i,2)==0
%         %     Xm(i,:) = Xm(i,:)+0.5; % Rotate 90 degrees
%         Ym(:,i) = Ym(:,i)+0.5;
%     end
%     
%     % Perturb each node according to normal distribution
%     seed = randn(Num+1,1);
%     Xm(:,i) = Xm(:,i) + lambda*cos(seed);
%     Ym(:,i) = Ym(:,i) + lambda*sin(seed);
%     
%     for j = 1:Num+1
%         if Xm(j,i) <= 0 || Xm(j,i) >= Num ||...
%                 Ym(j,i) <= 0 || Ym(j,i) >= Num
%             Xm(j,i) = 0;
%             Ym(j,i) = 0;
%         end
%     end
% end
% 
% % Reshape the centroids matrix into array
% Xref = reshape(Xm,(Num+1)^2,1);
% Yref = reshape(Ym,(Num+1)^2,1);
% % Pick-out centroids outside of periodic box
% X = nonzeros(Xref);
% Y = nonzeros(Yref);
% 
% L = length(X);
% 
% % Replicate the periodic box on 8 neighboring positions
% Ms(1:L,1) = X - Num; Ms(1:L,2) = Y + Num;
% Ms(L+1:2*L,1) = X; Ms(L+1:2*L,2) = Y + Num;
% Ms(2*L+1:3*L,1) = X + Num; Ms(2*L+1:3*L,2) = Y + Num;
% Ms(3*L+1:4*L,1) = X - Num; Ms(3*L+1:4*L,2) = Y;
% Ms(4*L+1:5*L,1) = X; Ms(4*L+1:5*L,2) = Y; % This part is what we wanted
% Ms(5*L+1:6*L,1) = X + Num; Ms(5*L+1:6*L,2) = Y;
% Ms(6*L+1:7*L,1) = X - Num; Ms(6*L+1:7*L,2) = Y - Num;
% Ms(7*L+1:8*L,1) = X; Ms(7*L+1:8*L,2) = Y - Num;
% Ms(8*L+1:9*L,1) = X + Num; Ms(8*L+1:9*L,2) = Y - Num;

%% Voronoi Diagram Generation

[V,C] = LloydOptimizePeri(Ms, tar, err, Num);
% [V,C]=voronoin(Ms);

V = ceil(V*100)/100;

figure,
voronoi(Ms(:,1),Ms(:,2)); % Show the Voronoi Tesselation
for i = 4*L+1:5*L
    hold on, patch(V(C{i,1},1),V(C{i,1},2),'y')
end
hold on, plot([0 Num Num 0 0],[0 0 Num Num 0],'r-','LineWidth',1.5)
hold on, plot(X,Y,'co','LineWidth',3)
title('Random Lattice Generated')
set(gca,'FontSize',14)
axis equal
grid on
saveas(gcf,'PeriodicRandomLattice.fig')
saveas(gcf,'PeriodicRandomLattice.jpg')

figure,
for i = 4*L+1:5*L
    hold on, patch(V(C{i},1),V(C{i},2),'y')
    hold on, plot(V(C{i},1),V(C{i},2),'b*')
end
hold on, plot([0 Num Num 0 0],[0 0 Num Num 0],'r-','LineWidth',1.5)
hold on, plot(X,Y,'co','LineWidth',3)
axis equal
grid on
set(gca,'FontSize',14)

% Save cells within the periodic bx
Cell = C(4*L+1:5*L);
Vflag = zeros(length(V),1);
for i = 4*L+1:5*L
    Vflag(C{i}) = 1; % Mark vertices which are within the box!
end

%% Organize orientation of each cell
for k = 1:length(Cell)
    Vt = V(Cell{k},:);
    Vt = Vt - mean(Vt);
    [~,ID] = ShapeOrder(Vt);
    Cell{k} = Cell{k}(ID);
end


%% Clear useless vertices
[m,~,~] = find(Vflag == 1);
for i = 1:length(V)
    if Vflag(i) ~= 1 % Remove unused boxing
        V(i,:) = 0;
    else % Remove repeative vertices
        ids = find(ismember(V,V(i,:),'rows') == 1);
        
        ls = length(ids);
        if ls ~= 1 % Repeative Vertices needs replacement
            Vref = ids(1);
            for ii = 2:ls
                Vrep = ids(ii);
                
                for k = 1:length(Cell)
                    out = find(Cell{k} == Vrep);
                    if ~isempty(out)
                        Cell{k}(out) = Vref;
                    end
                end
                % Clear the replaced vertex
                V(Vrep,:) = [0 0];
            end    
            
        end
    end
    
end


%% Generating Edges
n=1; sumE = 0;
Siz = zeros(length(Cell),1);
for k = 1:length(Cell)
    Km = length(Cell{k});% Number of nodes/edges in cell number k
    for K = 1:Km-1
        edge(sumE+K,1:2) = [Cell{k}(K),Cell{k}(K+1)];
    end
    edge(sumE+Km,1:2) = [Cell{k}(Km),Cell{k}(1)];
    sumE = sumE+Km;
    
    Siz(k) = length(Cell{k});
end

ne=1:length(edge);


% Calculate Cell Volume(Cross-Section Area)
for i = 1:length(Cell)
    Siz(i) = length(Cell{i});
        % For each cell, go through its vertices and compute area with
        % cross product
        vol = 0;
        for ss = 1:Siz(i)-1
            vol = vol + V(Cell{i}(ss),1)*V(Cell{i}(ss+1),2) - V(Cell{i}(ss),2)*V(Cell{i}(ss+1),1);
        end
        vol = vol + V(Cell{i}(Siz(i)),1)*V(Cell{i}(1),2) - V(Cell{i}(Siz(i)),2)*V(Cell{i}(1),1);
        
        % Check if the result is possitive or negative
        if vol>0
            Vol(i) = vol/2;
            Sign(i) = 1;
        else
            Vol(i) = -vol/2;
            Sign(i) = -1;
        end
end

% Compute standard deviation of area (Gamma Dist.)
% CA = std(Vol);
Coeff = gamfit(Vol);
CA = sqrt(Coeff(1)*Coeff(2)^2);
hold on, title(['Periodic Box, c_{A} = ',num2str(CA)])
saveas(gcf,'PeriodicBox.fig')
saveas(gcf,'PeriodicBox.jpg')

% ================= Include new volume from gamma distribution ================= 
% testV = randg(tar^-2, [length(Vol),1]);
% testV = testV/mean(testV);
% testV(end) = Num^2 - sum(testV(1:end-1));
% [~,Vid] = sort(Vol,'ascend');
% [~,VTid] = sort(testV,'ascend');
% VolTar(Vid) = testV(VTid);
% % hist(VolTar)
% Coeff = gamfit(VolTar);
% CA = sqrt(Coeff(1)*Coeff(2)^2);
% hold on, title(['Periodic Box, c_{A} = ',num2str(CA)])
% saveas(gcf,'PeriodicBox.fig')
% saveas(gcf,'PeriodicBox.jpg')
% 
% Vol = VolTar;
% ================= Include new volume from gamma distribution ================= 

% Compute target perimeter of each cell
Ptar = 2*sqrt(pi*Vol); % target perimeter of each cell
PtarAvg = mean(Ptar); % margin length

%% Merging Edges
ind = 1;
valence = ones(length(edge),1);
for M = 1:length(edge)
    Buff = edge(M,1:2);
    for N = 1:M-1
        if Buff == edge(N,1:2)
            ne(M) = N;
            edge(M,1:2) = [0,0];
            valence(M) = 2;
            valence(N) = 2;
        elseif Buff == flip(edge(N,1:2))
            ne(M) = -N;
            edge(M,1:2) = [0,0];
            valence(M) = 2;
            valence(N) = 2;
        else
            %             edge(ind) = e(M);
            %             ind = ind + 1;
        end
    end
end

%% Pair-up Vertices
Vnew = V;                   % Altered value --- buff: zeros(length(V),2)
Vpair = zeros(length(V),2); % Pair of changed vertices
Vc = zeros(length(V),1);    % Mark of vertex warp
Vtar = zeros(1,2);          % Target vertex coordinates
for i = 1:length(V)
    if Vflag(i) % We worry only those vertices within periodic box
        
        if V(i,1)>=Num && V(i,2)<Num && V(i,2)>0% Right
            Vtar(1) = V(i,1) - Num; % Target X
            Vtar(2) = V(i,2);       % Target Y
            Vc(i) = 1;              % Mark this vertex as warped
            
            % Rounding
            Vtar(1) = round(Vtar(1)*100)/100;
            
            % Check if the target vertex is within the box
            [Mark,Tar] = ismember(Vtar,Vnew,'rows');
            
            if Mark == 0 % If there's NO corresponding vertex inside
                Vpair(i,:) = [i,i]; % Pair with itself
                Vnew(i,:) = Vtar; % Update actual coordinate
            else % If there IS corresponding vertex inside
                Vpair(i,:) = [Tar,i]; % Pair with the other
                Vnew(i,:) = [0,0]; % Update actual coordinate to 0
            end
            
        elseif V(i,1)<0 && V(i,2)<Num && V(i,2)>0 % Left
            Vtar(1) = V(i,1) + Num; % Target X
            Vtar(2) = V(i,2);       % Target Y
            Vc(i) = 3;              % Mark this vertex as warped
            
            % Check if the target vertex is within the box
            [Mark,Tar] = ismember(Vtar,Vnew,'rows');
            
            if Mark == 0 % If there's NO corresponding vertex inside
                Vpair(i,:) = [i,i]; % Pair with itself
                Vnew(i,:) = Vtar; % Update actual coordinate
            else % If there IS corresponding vertex inside
                Vpair(i,:) = [Tar,i]; % Pair with the other
                Vnew(i,:) = [0,0]; % Update actual coordinate to 0
            end
            
        elseif V(i,2)>=Num && V(i,1)<Num && V(i,1)>0 % Top
            Vtar(1) = V(i,1);       % Target X
            Vtar(2) = V(i,2) - Num; % Target Y
            Vc(i) = 2;              % Mark this vertex as warped
            
            % Rounding
            Vtar(2) = round(Vtar(2)*100)/100;
            
            % Check if the target vertex is within the box
            [Mark,Tar] = ismember(Vtar,Vnew,'rows');
            
            if Mark == 0 % If there's NO corresponding vertex inside
                Vpair(i,:) = [i,i]; % Pair with itself
                Vnew(i,:) = Vtar; % Update actual coordinate
            else % If there IS corresponding vertex inside
                Vpair(i,:) = [Tar,i]; % Pair with the other
                Vnew(i,:) = [0,0]; % Update actual coordinate to 0
            end
            
        elseif V(i,2)<0 && V(i,1)<Num && V(i,1)>0 % Bottom
            Vtar(1) = V(i,1);       % Target X
            Vtar(2) = V(i,2) + Num; % Target Y
            Vc(i) = 4;              % Mark this vertex as warped
            
            % Check if the target vertex is within the box
            [Mark,Tar] = ismember(Vtar,Vnew,'rows');
            
            if Mark == 0 % If there's NO corresponding vertex inside
                Vpair(i,:) = [i,i]; % Pair with itself
                Vnew(i,:) = Vtar; % Update actual coordinate
            else % If there IS corresponding vertex inside
                Vpair(i,:) = [Tar,i]; % Pair with the other
                Vnew(i,:) = [0,0]; % Update actual coordinate to 0
            end
            
        elseif V(i,1)>Num && V(i,2)>Num % Top Right
            Vtar(1) = V(i,1) - Num; % Target X
            Vtar(2) = V(i,2) - Num; % Target Y
            Vc(i) = 5;              % Mark this vertex as warped
            
            % Rounding
            Vtar(1) = round(Vtar(1)*100)/100;
            Vtar(2) = round(Vtar(2)*100)/100;
            
            % Check if the target vertex is within the box
            [Mark,Tar] = ismember(Vtar,Vnew,'rows');
            
            if Mark == 0 % If there's NO corresponding vertex inside
                Vpair(i,:) = [i,i]; % Pair with itself
                Vnew(i,:) = Vtar; % Update actual coordinate
            else % If there IS corresponding vertex inside
                Vpair(i,:) = [Tar,i]; % Pair with the other
                Vnew(i,:) = [0,0]; % Update actual coordinate to 0
            end
            
        elseif V(i,1)>Num && V(i,2)<0 % Bottom Right
            Vtar(1) = V(i,1) - Num; % Target X
            Vtar(2) = V(i,2) + Num; % Target Y
            Vc(i) = 6;              % Mark this vertex as warped
            
            % Rounding
            Vtar(1) = round(Vtar(1)*100)/100;
            
            % Check if the target vertex is within the box
            [Mark,Tar] = ismember(Vtar,Vnew,'rows');
            
            if Mark == 0 % If there's NO corresponding vertex inside
                Vpair(i,:) = [i,i]; % Pair with itself
                Vnew(i,:) = Vtar; % Update actual coordinate
            else % If there IS corresponding vertex inside
                Vpair(i,:) = [Tar,i]; % Pair with the other
                Vnew(i,:) = [0,0]; % Update actual coordinate to 0
            end
            
        elseif V(i,1)<0 && V(i,2)>Num % Top Left
            Vtar(1) = V(i,1) + Num; % Target X
            Vtar(2) = V(i,2) - Num; % Target Y
            Vc(i) = 7;              % Mark this vertex as warped
            
            % Rounding
            Vtar(2) = round(Vtar(2)*100)/100;
            
            % Check if the target vertex is within the box
            [Mark,Tar] = ismember(Vtar,Vnew,'rows');
            
            if Mark == 0 % If there's NO corresponding vertex inside
                Vpair(i,:) = [i,i]; % Pair with itself
                Vnew(i,:) = Vtar; % Update actual coordinate
            else % If there IS corresponding vertex inside
                Vpair(i,:) = [Tar,i]; % Pair with the other
                Vnew(i,:) = [0,0]; % Update actual coordinate to 0
            end
            
        elseif V(i,1)<0 && V(i,2)<0 % Bottom Left
            Vtar(1) = V(i,1) + Num; % Target X
            Vtar(2) = V(i,2) + Num; % Target Y
            Vc(i) = 8;              % Mark this vertex as warped
            
            % Check if the target vertex is within the box
            [Mark,Tar] = ismember(Vtar,Vnew,'rows');
            
            if Mark == 0 % If there's NO corresponding vertex inside
                Vpair(i,:) = [i,i]; % Pair with itself
                Vnew(i,:) = Vtar; % Update actual coordinate
            else % If there IS corresponding vertex inside
                Vpair(i,:) = [Tar,i]; % Pair with the other
                Vnew(i,:) = [0,0]; % Update actual coordinate to 0
            end
        end
    end
end
% % * Value Table for Vc (Vertex Warp Info):
% %     Right --- 1; Top --- 2; Left --- 3; Bottom --- 4
% %     Top Right --- 5; Bottom Right --- 6;
% %     Top Left --- 7;  Bottom Left --- 8;

%% Calculate edge lengths
EL = zeros(length(edge),1);
for i = 1:length(edge)
    if edge(i,1)~=0 && edge(i,2)~=0
        EL(i) = norm([(V(edge(i,1),1) - V(edge(i,2),1)),...
                      (V(edge(i,1),2) - V(edge(i,2),2))]);
    end
end
EL = round(EL*10000)/10000;

%% Warp Periodic Edges
% Generate Pair up symmetric edges
Vtar1 = zeros(1,2);
Vtar2 = zeros(1,2);
flag1 = 0; flag2 = 0;
% Warp Direction array "Dir" of each edge:
%    [0 0] --- no warp
%    [1 0] --- +x ;    [-1 0] --- -x
%    [0 1] --- +y ;    [0 -1] --- -y
%    [1 1] --- +x+y;   [1 -1] --- +x-y
%    [-1 1]--- -x+y;   [-1 -1]--- -x-y
Dir = zeros(length(edge),2);
% Symmetry pair of edges
% Epair = zeros(length(edge),2);
Enew = edge; neS = ne; VLocbuff  = zeros(length(edge),1);
for i = 1:length(edge)
    if Enew(i,1)~=0 && Enew(i,2)~=0
        V1 = edge(i,1);
        V2 = edge(i,2);
        flag1 = Vc(V1); flag2 = Vc(V2); % Get warp info on vertex
        
% Part 1 % ======= % Update Edge Info % ======= %
        if flag1 ~= 0 && flag2 == 0
            Enew(i,1) = Vpair(V1,1);
            Enew(i,2) = V2;
        elseif flag1 == 0 && flag2 ~= 0
            Enew(i,1) = V1;
            Enew(i,2) = Vpair(V2,1);
        elseif flag1 == 0 && flag2 == 0
            Enew(i,1) = V1;
            Enew(i,2) = V2;
        else
            Enew(i,1) = Vpair(V1,1);
            Enew(i,2) = Vpair(V2,1);
        end
        
% Part 2 % ======= % Determine Warp Form of edge: % ======= %
        if     (flag1 == 0 && flag2 == 1) || (flag1 == 2 && flag2 == 5) || ...
               (flag1 == 4 && flag2 == 6) || (flag1 == 3 && flag2 == 0) || ...
               (flag1 == 7 && flag2 == 2) || (flag1 == 8 && flag2 == 4)
           % ====== % +X % ====== % (BEGIN)
           Dir(i,:) = [1 0]; % Warp Form of the edge
           % ====== % +X % ====== %
        elseif (flag1 == 0 && flag2 == 3) || (flag1 == 2 && flag2 == 7) || ...
               (flag1 == 4 && flag2 == 8) || (flag1 == 1 && flag2 == 0) || ...
               (flag1 == 5 && flag2 == 2) || (flag1 == 6 && flag2 == 4)
           % ====== % -X % ====== % (BEGIN)
           Dir(i,:) = [-1 0];
           % ====== % -X % ====== %
        elseif (flag1 == 0 && flag2 == 2) || (flag1 == 1 && flag2 == 5) || ...
               (flag1 == 3 && flag2 == 7) || (flag1 == 4 && flag2 == 0) || ...
               (flag1 == 6 && flag2 == 1) || (flag1 == 8 && flag2 == 3)
           % ====== % +Y % ====== % (BEGIN)
           Dir(i,:) = [0 1];
           % ====== % +Y % ====== %
        elseif (flag1 == 0 && flag2 == 4) || (flag1 == 1 && flag2 == 6) || ...
               (flag1 == 3 && flag2 == 8) || (flag1 == 2 && flag2 == 0) || ...
               (flag1 == 5 && flag2 == 1) || (flag1 == 7 && flag2 == 3)
           % ====== % -Y % ====== % (BEGIN)
           Dir(i,:) = [0 -1];
           % ====== % -Y % ====== %
        elseif (flag1 == 0 && flag2 == 5) || (flag1 == 4 && flag2 == 1) || ...
               (flag1 == 8 && flag2 == 0) || (flag1 == 3 && flag2 == 2)
           % ====== % +X +Y % ====== % (BEGIN)
           Dir(i,:) = [1 1];
           % ====== % +X +Y % ====== %
        elseif (flag1 == 0 && flag2 == 6) || (flag1 == 2 && flag2 == 1) || ...
               (flag1 == 3 && flag2 == 4) || (flag1 == 7 && flag2 == 0)
           % ====== % +X -Y % ====== % (BEGIN)
           Dir(i,:) = [1 -1];
           % ====== % +X -Y % ====== %
        elseif (flag1 == 0 && flag2 == 7) || (flag1 == 1 && flag2 == 2) || ...
               (flag1 == 4 && flag2 == 3) || (flag1 == 6 && flag2 == 0)
           % ====== % -X +Y % ====== % (BEGIN)
           Dir(i,:) = [-1 1];
           % ====== % -X +Y % ====== %
        elseif (flag1 == 0 && flag2 == 8) || (flag1 == 1 && flag2 == 4) || ...
               (flag1 == 5 && flag2 == 0) || (flag1 == 2 && flag2 == 3)
           % ====== % -X -Y % ====== % (BEGIN)
           Dir(i,:) = [-1 -1];
           % ====== % -X -Y % ====== %
        end
        
% Part 3 % ======= % Update Edge-Face Info % ======= %
        if valence(i) == 1 && EL(i) ~= 0%&& Enew(i,1) ~= 0 && Enew(i,2) ~= 0
            % Look for mirror vertex by checking edge length
            ELtemp = EL;
            ELtemp(1:i) = 0;
            
            for j = 1:length(EL)
                if valence(j) == 2
                    ELtemp(j) = 0;
                end
            end
            % Find other edges with same length
            [VLoc,~,~] = find(ELtemp == EL(i));
            % If we have more than 1 match, check for target with slope
            if length(VLoc) > 1
%                 XLoc = abs(V(edge(VLoc,1),1) - V(edge(VLoc,2),1));
%                 YLoc = abs(V(edge(VLoc,1),2) - V(edge(VLoc,2),2));
%                 Xi = abs(V(edge(i,1),1) - V(edge(i,2),1));
%                 Yi = abs(V(edge(i,1),2) - V(edge(i,2),2));

                XLoc = (V(edge(VLoc,1),1) - V(edge(VLoc,2),1));
                YLoc = (V(edge(VLoc,1),2) - V(edge(VLoc,2),2));
                Xi = (V(edge(i,1),1) - V(edge(i,2),1));
                Yi = (V(edge(i,1),2) - V(edge(i,2),2));
                
                XLoc = round(XLoc*10000)/10000;
                YLoc = round(YLoc*10000)/10000;
                Xi = round(Xi*10000)/10000;
                Yi = round(Yi*10000)/10000;
                
                for j = 1:length(VLoc)
                    % Only match would be the one with exact the same
                    % slope
                    if (XLoc(j) == Xi && YLoc(j) == Yi) || ...
                            (XLoc(j) == -Xi && YLoc(j) == -Yi)    %XLoc(j) == Xi && YLoc(j) == Yi
                        VLocF = VLoc(j);
                    end
                end
                
                VLoc = []; VLoc = VLocF;
            end
            % Target Edge grabbed!
%             Vmir = edge(VLoc);
            
                XLoc = (V(edge(VLoc,1),1) - V(edge(VLoc,2),1));
                YLoc = (V(edge(VLoc,1),2) - V(edge(VLoc,2),2));
                Xi = (V(edge(i,1),1) - V(edge(i,2),1));
                Yi = (V(edge(i,1),2) - V(edge(i,2),2));
                
                XLoc = round(XLoc*10000)/10000;
                YLoc = round(YLoc*10000)/10000;
                Xi = round(Xi*10000)/10000;
                Yi = round(Yi*10000)/10000;
                
            % Find orientation of the final target edge and update info
%             if ~isempty(VLoc) && VLoc ~= i
%                 Enew(VLoc,:) = [0,0];
%                 VLocbuff(i) = VLoc;
%                 ne(VLoc) = [-1*ne(i)];
%             end
            if ~isempty(VLoc) && VLoc > i%~=
                Enew(VLoc,:) = [0,0];
                VLocbuff(i) = VLoc;
                %             Epair(i,1) = [-1*ne(i)];
                if (XLoc == -Xi && YLoc == -Yi) %|| ... %ne(VLoc) > 0 && 
                    %(ne(VLoc) < 0 && XLoc == Xi && YLoc == Yi)
                    ne(VLoc) = [-1*ne(i)];%Epair(i,1) = [-1*ne(i)];%
                elseif (XLoc == Xi && YLoc == Yi) %|| ... %ne(VLoc) > 0 && 
                    %(ne(VLoc) < 0 && XLoc == -Xi && YLoc == -Yi)
                    ne(VLoc) = [ne(i)];%Epair(i,1) = [ne(i)];%
                end
            end
            
%             if ~isempty(VLoc) && (V(edge(VLoc,1),1) - V(edge(VLoc,2),1)) == -(V(edge(i,1),1) - V(edge(i,2),1)) && ...
%                     (V(edge(VLoc,1),2) - V(edge(VLoc,2),2)) == -(V(edge(i,1),2) - V(edge(i,2),2))
% %                 Enew(VLoc,:) = [0,0];
% %                 EL(VLoc) = 0;
%                 Epair(i,1) = [-1*ne(i)];
%             elseif ~isempty(VLoc) && (V(edge(VLoc,1),1) - V(edge(VLoc,2),1)) == (V(edge(i,1),1) - V(edge(i,2),1)) && ...
%                     (V(edge(VLoc,1),2) - V(edge(VLoc,2),2)) == (V(edge(i,1),2) - V(edge(i,2),2))
% %                 Enew(VLoc,:) = [0,0];
% %                 EL(VLoc) = 0;
%                 Epair(i,1) = [ne(i)];
%             end
        end
        
% % ======= % End of Checking         
    end
end
% Epair(:,2) = ne;
neS = neS';
neS(:,2) = ne';
Eorg = [Enew,ne'];

for i =1:length(Eorg)
    if Eorg(i,1) == 0 && Eorg(i,2) == 0 && Eorg(i,3) == i
        if Vpair(edge(i,1),2) == 0 && Vpair(edge(i,2),2) == 0
            Enew(i,:) = edge(i,:);
        elseif Vpair(edge(i,1),2) ~= 0 && Vpair(edge(i,2),2) == 0% || Vpair(edge(i,2),2) ~= 0
            Vbuff1 = Vpair(edge(i,1),1);
            Vbuff2 = edge(i,2);
            
            [EMark1,ETar1] = ismember([Vbuff1,Vbuff2],Enew,'rows');
            [EMark2,ETar2] = ismember([Vbuff2,Vbuff1],Enew,'rows');
            
            if EMark1 == 0
                Enew(i,:) = Enew(ETar2,:);
            elseif EMark2 == 0
                Enew(i,:) = Enew(ETar1,:);
            end
        elseif Vpair(edge(i,1),2) == 0 && Vpair(edge(i,2),2) ~= 0
            Vbuff1 = edge(i,1);
            Vbuff2 = Vpair(edge(i,2),1);
            
            [EMark1,ETar1] = ismember([Vbuff1,Vbuff2],Enew,'rows');
            [EMark2,ETar2] = ismember([Vbuff2,Vbuff1],Enew,'rows');
            
            if EMark1 == 0
                Enew(i,:) = Enew(ETar2,:);
            elseif EMark2 == 0
                Enew(i,:) = Enew(ETar1,:);
            end
        elseif Vpair(edge(i,1),2) ~= 0 && Vpair(edge(i,2),2) ~= 0
            Vbuff1 = Vpair(edge(i,1),1);
            Vbuff2 = Vpair(edge(i,2),1);
            
            [EMark1,ETar1] = ismember([Vbuff1,Vbuff2],Enew,'rows');
            [EMark2,ETar2] = ismember([Vbuff2,Vbuff1],Enew,'rows');
            
            if EMark1 == 0
                Enew(i,:) = Enew(ETar2,:);
            elseif EMark2 == 0
                Enew(i,:) = Enew(ETar1,:);
            end
        end
    end
end

%% Constructing Faces
% Still need work to make it fit for more distorted case

Face = cell(size(Cell));
Lc = 0; %Face = [];
for CF1 = 1:length(Siz)
    
    % ==== Using array to store Faces ==== %
    %     for CF2 = 1:Siz(CF1)
    %         Face(CF1,CF2) = ne(Lcell + CF2);
    %     end
    % ==== Using array to store Faces ==== %
    
    Face{CF1} = ne(Lc + 1:Lc + Siz(CF1));
    % Using Cell array to store Face information
    Lc = Lc + Siz(CF1);
end

%% Statistical Accessment

figure,
histogram(Vol,20)
grid on
set(gca,'FontSize',14)
title('Statistical Distribution of Cellular Volume')
xlabel('V')
ylabel('N')
saveas(gcf,'VHisto.fig')
saveas(gcf,'VHisto.jpg')

Pref = 2*sqrt(pi*Vol);

for i =1:length(Cell)
    xlen(i) = length(Cell{i});
end
figure,
histogram(xlen)
grid on
set(gca,'FontSize',14)
title('Polygon Type')
xlabel('Number of Edges')
ylabel('N')

% text(4.8,100,'Pentagon')
% text(5.8,370,'Hexagon')
% text(6.8,70,'Heptagon')

saveas(gcf,'Polygon.fig')
saveas(gcf,'Polygon.jpg')

%% Outputing into .FE file
% %% Output Flow --- phase file
% fileID1 = fopen('phaseSim.txt','w');
% 
% fprintf(fileID1,'2\n');
% fprintf(fileID1,'%1d  %1d  %4.3f\n',1,1,0.001);%2*beta - gamma1
% fprintf(fileID1,'%1d  %1d  %4.3f\n',1,2,1);%2*beta - (gamma1 + gamma2)/2
% fprintf(fileID1,'%1d  %1d  %4.3f\n',2,2,2*beta - gamma2);

%% Output Flow --- To Surface Evolver Directly
fileID = fopen('randomGen.fe','w');

% Headers
fprintf(fileID,'STRING\n');
% fprintf(fileID,'quadratic\n');
fprintf(fileID,'space_dimension 2 \n');
fprintf(fileID,'\n');
fprintf(fileID,'torus_filled\n');
fprintf(fileID,'periods\n');
fprintf(fileID,[num2str(Num),' ',num2str(0),'\n']);
fprintf(fileID,[num2str(0),' ',num2str(Num),'\n']);
fprintf(fileID,'\n');
fprintf(fileID,'parameter beta = 0.0 \n');
fprintf(fileID,'parameter gamma = 0.22 \n');
% fprintf(fileID,'\n');
% fprintf(fileID,'\n');

% fprintf(fileID,'quantity SurfArea info_only method edge_length\n');
% fprintf(fileID,'quantity BulkArea info_only method edge_length\n');
fprintf(fileID,'\n');

fprintf(fileID,'define edge attribute Cell real[2]\n');
fprintf(fileID,'define facet attribute Pref real\n');
fprintf(fileID,'\n');
% fprintf(fileID,'define facet attribute warp real\n');

% ======== Energy Functional Apply here! ======== %
% % for i=1:length(Face)
% % %     if Phase(i) == 2
% %         fprintf(fileID,['method_instance Peri',num2str(i),' method edge_general_integral\n']);
% %         fprintf(fileID,['	scalar_integrand: sqrt(x3^2 + x4^2)\n']);
% %         fprintf(fileID,['quantity contract',num2str(i),' energy function 0.5*(Peri',num2str(i),'.value - ',num2str((1+r/2)*Ptar(i)),')^2/',num2str(Ptar(i)),'\n']);
% %         fprintf(fileID,'\n');
% % %     end
% % end
% ======== Energy Functional Apply here! ======== %

% ======== X,Y Moments Calculation ======== %
% % for i=1:length(Face) % % X-Moment
% % %     if Phase(i) == 2
% %         fprintf(fileID,['quantity xmoment',num2str(i),' info_only method facet_scalar_integral\n']);
% %         fprintf(fileID,['	scalar_integrand: y^2\n']);
% %         fprintf(fileID,'\n');
% % %     end
% % end
% % 
% % for i=1:length(Face) % % Y-Moment
% % %     if Phase(i) == 2
% %         fprintf(fileID,['quantity ymoment',num2str(i),' info_only method facet_scalar_integral\n']);
% %         fprintf(fileID,['	scalar_integrand: x^2\n']);
% %         fprintf(fileID,'\n');
% % %     end
% % end

% % for i=1:length(Face) % % X-Moment
% % %     if Phase(i) == 2
% %         fprintf(fileID,['quantity xmoment',num2str(i),' info_only method facet_vector_integral\n']);
% %         fprintf(fileID,['	vector_integrand:\n']);
% %         fprintf(fileID,['	q1: y^2\n']);
% %         fprintf(fileID,['	q2: 0\n']);
% %         fprintf(fileID,'\n');
% % %     end
% % end
% % 
% % for i=1:length(Face) % % Y-Moment
% % %     if Phase(i) == 2
% %         fprintf(fileID,['quantity ymoment',num2str(i),' info_only method facet_vector_integral\n']);
% %         fprintf(fileID,['	vector_integrand:\n']);
% %         fprintf(fileID,['	q1: 0\n']);
% %         fprintf(fileID,['	q2: x^2\n']);
% %         fprintf(fileID,'\n');
% % %     end
% % end
% ======== X,Y Moments Calculation ======== %

% fprintf(fileID,'quantity mci info_only method mean_curvature_integral');
% fprintf(fileID,'\n');

% % Correct Vertices
% for i = 1:length(Vnew)
%     if Vnew(i,1) ~= 0 && Vnew(i,1) ~= 0 && i > 1
%         % Find if this vertex is repeated already
%         for j = 1:i-1
%             if Vnew(j,1) == Vnew(i,1) && Vnew(j,2) == Vnew(i,2)
%                 % If repeated, clear later vertex i and keep j
%                 Vnew(i,:) = [0 0];
%                 % Find and change correspoding edge constitution
%                 for k = 1:length(Enew)
%                     if ~isempty(find(Enew(k,:) == i))
%                         id =  find(Enew(k,:) == i);
%                         Enew(k,id) = j;
%                     end
%                 end
%                 % Check and correct warping array if necessary
%             end
%         end
%     end
% end

% Vertices
Vert = Vnew;
 fprintf(fileID,'vertices \n');
 for i=1:length(Vert)
     if Vert(i,1) ~= 0 || Vert(i,2) ~= 0%Vert(i,:) ~= [0 0]
         fprintf(fileID,'%1d        %4.3f %4.3f \n',i,Vert(i,:));%botcon topcon
     end
 end
 fprintf(fileID,'\n');
 
 % Correct Edges
for i = 1:length(Enew)
    if Enew(i,1) ~= 0 && Enew(i,1) ~= 0 && i > 1
        % Find if this edge is repeated already
        for j = 1:i-1
            if Enew(j,1) == Enew(i,1) && Enew(j,2) == Enew(i,2)
                % If repeated, clear later edge i and keep j
                Enew(i,:) = [0 0];
                % Find and change correspoding face constitution
                for k = 1:length(Face)
                    if ~isempty(find(abs(Face{k}) == i))
                        id =  find(abs(Face{k}) == i);
                        Face{k}(id) = j*(Face{k}(id)/i);
                    end
                end
                % Check and correct warping array if necessary
                
            elseif Enew(j,1) == Enew(i,2) && Enew(j,2) == Enew(i,1)
                % If repeated, clear later edge i and keep j
                Enew(i,:) = [0 0];
                % Find and change correspoding face constitution
                for k = 1:length(Face)
                    if ~isempty(find(abs(Face{k}) == i))
                        id =  find(abs(Face{k}) == i);
                        Face{k}(id) = -j*(Face{k}(id)/i);
                    end
                end
                % Check and correct warping array if necessary
                
            end
        end
    end
end

% Edges
% Warp Direction array "Dir" of each edge:
%    [0 0] --- no warp
%    [1 0] --- +x ;    [-1 0] --- -x
%    [0 1] --- +y ;    [0 -1] --- -y
%    [1 1] --- +x+y;   [1 -1] --- +x-y
%    [-1 1]--- -x+y;   [-1 -1]--- -x-y
fprintf(fileID,'edges \n');
for i=1:length(Enew)
    if Enew(i,1) ~= 0 || Enew(i,2) ~= 0
        
        
        if Dir(i,:) == [1 0]
            array = ' + *';
            fprintf(fileID,['%1d        %1d    %1d',array,' tension beta  \n'],i,Enew(i,1),Enew(i,2));
        elseif Dir(i,:) == [-1  0]
            array = ' - *';
            fprintf(fileID,['%1d        %1d    %1d',array,' tension beta  \n'],i,Enew(i,1),Enew(i,2));
        elseif Dir(i,:) == [ 0  1]
            array = ' * +';
            fprintf(fileID,['%1d        %1d    %1d',array,' tension beta  \n'],i,Enew(i,1),Enew(i,2));
        elseif Dir(i,:) == [ 0 -1]
            array = ' * -';
            fprintf(fileID,['%1d        %1d    %1d',array,' tension beta  \n'],i,Enew(i,1),Enew(i,2));
        elseif Dir(i,:) == [ 1  1]
            array = ' + +';
            fprintf(fileID,['%1d        %1d    %1d',array,' tension beta  \n'],i,Enew(i,1),Enew(i,2));
        elseif Dir(i,:) == [ 1 -1]
            array = ' + -';
            fprintf(fileID,['%1d        %1d    %1d',array,' tension beta  \n'],i,Enew(i,1),Enew(i,2));
        elseif Dir(i,:) == [-1  1]
            array = ' - +';
            fprintf(fileID,['%1d        %1d    %1d',array,' tension beta  \n'],i,Enew(i,1),Enew(i,2));
        elseif Dir(i,:) == [-1 -1]
            array = ' - -';
            fprintf(fileID,['%1d        %1d    %1d',array,' tension beta  \n'],i,Enew(i,1),Enew(i,2));
        else
            array = ' * *';
            fprintf(fileID,['%1d        %1d    %1d',array,' tension beta  \n'],i,Enew(i,1),Enew(i,2));
        end
    end
end
fprintf(fileID,'\n');

% Faces
fprintf(fileID,'faces \n');
for i=1:length(Face)
    array = [];
    for j = 1:length(Face{i})
        array = [array, ' %d'];
    end
    
    fprintf(fileID,['%1d       ', array,'\n'],i, Face{i});
end
fprintf(fileID,'\n');

% Bodies
fprintf(fileID,'bodies \n');
for i=1:length(Face)
    if Sign(i)<0
        fprintf(fileID,'%1d        -%1d       volume  %4.3f\n',i, i,Vol(i));
    elseif Sign(i)>0
        fprintf(fileID,'%1d        %1d       volume  %4.3f\n',i, i,Vol(i));%volume  %4.3f %,Vol(i) %pressure 0
    end
end
fprintf(fileID,'\n');

fprintf(fileID,'read\n');
% fprintf(fileID,'\n');
% fprintf(fileID,'generate := {\n');
% for i = 1:length(Cell)
%     fprintf(fileID,['    body[',num2str(i),'].target := ',num2str(VolTar(i)),';\n']);
% end
% fprintf(fileID,'	}\n');
fprintf(fileID,'\n');
fprintf(fileID,'ini := {\n');
fprintf(fileID,'\n');
fprintf(fileID,['    quadratic;\n']);
fprintf(fileID,'\n');
fprintf(fileID,'    define facet attribute perimeter real function\n');
fprintf(fileID,'        {self.perimeter := sum(self.edge,length);}; //Define facet attribute "perimeter"\n');
fprintf(fileID,'\n');
for i = 1:length(Cell)
    fprintf(fileID,['    facet[',num2str(i),'].Pref := ',num2str(Ptar(i)),';\n']);
end
fprintf(fileID,'\n');
fprintf(fileID,['    foreach edge eee do {eee.Cell[1] := 0; eee.Cell[2] := 0;};\n']);
fprintf(fileID,'\n');
fprintf(fileID,['    foreach edge eee do eee.color := red;\n']);
fprintf(fileID,'\n');
fprintf(fileID,['    foreach facet ff do { exec sprintf "define quantity cell_|03d energy real; ",ff.id;};\n']);
fprintf(fileID,'	}\n');
fprintf(fileID,'\n');

% % fprintf(fileID,'// X-Moments\n');
% % for i = 1:length(Face) % % X-Moment
% % %     if Phase(i) == 2
% %         fprintf(fileID,['    foreach facet[',num2str(i),'].edge eee where valence = 2 do {set eee quantity xmoment',num2str(i),'};\n']);
% % %     end
% % end
% % fprintf(fileID,'\n');
% % 
% % fprintf(fileID,'// Y-Moments\n');
% % for i = 1:length(Face) % % Y-Moment
% % %     if Phase(i) == 2
% %         fprintf(fileID,['    foreach facet[',num2str(i),'].edge eee where valence = 2 do {set eee quantity ymoment',num2str(i),'};\n']);
% % %     end
% % end
% % fprintf(fileID,'\n');

% fprintf(fileID,'// Energy Functionals\n');
% for i = 1:length(Face)
% %     if Phase(i) == 2
%         fprintf(fileID,['    foreach facet[',num2str(i),'].edge eee where valence = 2 do {set eee method_instance Peri',num2str(i),'};\n']);
% %     end
% end

% % ======== Quenched Disorder ======== % %
modifier = normrnd(1,cp,[length(Cell),1]); % generate a normal distributed modifier with given cp
PtarQD = Ptar'.*modifier;

fprintf(fileID,'\n');
fprintf(fileID,'iniQD := {\n');
fprintf(fileID,'\n');
fprintf(fileID,['    quadratic;\n']);
fprintf(fileID,'\n');
fprintf(fileID,'    define facet attribute perimeter real function\n');
fprintf(fileID,'        {self.perimeter := sum(self.edge,length);}; //Define facet attribute "perimeter"\n');
fprintf(fileID,'\n');
for i = 1:length(Cell)
    fprintf(fileID,['    facet[',num2str(i),'].Pref := ',num2str(PtarQD(i)),';\n']);
end
fprintf(fileID,'\n');
fprintf(fileID,['    foreach edge eee do {eee.Cell[1] := 0; eee.Cell[2] := 0;};\n']);
fprintf(fileID,'\n');
fprintf(fileID,['    foreach edge eee do eee.color := red;\n']);
fprintf(fileID,'\n');
fprintf(fileID,['    foreach facet ff do { exec sprintf "define quantity cell_|03d energy real; ",ff.id;};\n']);
fprintf(fileID,'	}\n');
fprintf(fileID,'\n');

% % ======== Unset Energy Functional ======== % %
% fprintf(fileID,'uni := {\n');
% fprintf(fileID,'\n');
% fprintf(fileID,'// Unset\n');
% for i = 1:length(Face)
% %     if Phase(i) == 2
%         fprintf(fileID,['    unset edge method_instance Peri',num2str(i),';\n']);
% %     end
% end
% fprintf(fileID,'\n');
% 
% fprintf(fileID,'	}\n');
% fprintf(fileID,'\n');
% % ======== Unset Energy Functional ======== % %

fprintf(fileID,'gg := {\n');
fprintf(fileID,'\n');
% fprintf(fileID,'    quadratic;\n');
fprintf(fileID,['    // Average Area of Cell',num2str(PtarAvg),'\n']);
fprintf(fileID,['    // Margin ratio Beta',num2str(B),'\n']);
% fprintf(fileID,'	g 100000;	\n');

fprintf(fileID,['    foreach facet do trash := perimeter;\n']);
fprintf(fileID,['    \n']);
fprintf(fileID,['    foreach edge eee do {eee.Cell[1] := 0; eee.Cell[2] := 0;};\n']);
fprintf(fileID,['    \n']);
fprintf(fileID,['    foreach facet fff do \n']);
fprintf(fileID,['        { \n']);
fprintf(fileID,['          foreach fff.edge eee do \n']);
fprintf(fileID,['              { if eee.Cell[1] = 0 then eee.Cell[1] := (fff.perimeter - fff.Pref)/fff.Pref else eee.Cell[2] := (fff.perimeter - fff.Pref)/fff.Pref;};\n']);
fprintf(fileID,['        };\n']);
fprintf(fileID,['    \n']);
fprintf(fileID,['    foreach edge eee do eee.tension := eee.Cell[1] + eee.Cell[2] - gamma;\n']);
fprintf(fileID,['    \n']);
fprintf(fileID,['    foreach facet ff do {\n']);
fprintf(fileID,['        exec sprintf "set cell_|03d volconst (1+gamma/2)^2*|f-(1+gamma/2)*|f; trash := cell_|03d.value; ",ff.id,ff.pref,ff.perimeter,ff.id;\n']);
fprintf(fileID,['    };\n']);
fprintf(fileID,'}\n');
fprintf(fileID,'\n');

fprintf(fileID,'test := {\n');
fprintf(fileID,'     g; gg;\n');
fprintf(fileID,'     \n');
fprintf(fileID,'     var := 0; foreach facet fff do {var := var + (fff.valence - avg(facet,valence))^2;};\n');
fprintf(fileID,['     var := sqrt(var/',num2str(Num^2-1),')/avg(facet,valence);\n']);
fprintf(fileID,'     \n');
fprintf(fileID,'     print "Variance for # of neighbor: "; print var; print "shiftn";\n');
fprintf(fileID,'     print "mean for tension coeff.: "; print avg(edge, log(abs(tension))/log(10)); print "shiftn";\n');
fprintf(fileID,'}\n');
fprintf(fileID,'\n');

fprintf(fileID,'recordQD := {\n');
fprintf(fileID,'     foreach facet fff do\n');
fprintf(fileID,'     {\n');
fprintf(fileID,'     	foreach fff.edge eee do\n');
fprintf(fileID,'        {\n');
fprintf(fileID,'        printf "|d ", eee.oid >> "OUTPUT/Ten022QD/FaceEdge.txt";\n');
fprintf(fileID,'        printf "|d ", eee.oid >> "OUTPUT/QD/FaceEdge_Ten022QD.txt";\n');
fprintf(fileID,'        \n');
fprintf(fileID,'        printf "|d ", eee.id >> "OUTPUT/Ten022QD/EdgeInfo.txt";\n');
fprintf(fileID,'        printf "|d ", eee.id >> "OUTPUT/QD/EdgeInfo_Ten022QD.txt";\n');
fprintf(fileID,'        \n');
fprintf(fileID,'        foreach eee.vertex vvv do \n');
fprintf(fileID,'            {\n');
fprintf(fileID,'            printf "|d ", vvv.id >> "OUTPUT/Ten022QD/EdgeInfo.txt";\n');
fprintf(fileID,'            printf "|d ", vvv.id >> "OUTPUT/QD/EdgeInfo_Ten022QD.txt";\n');
fprintf(fileID,'            };\n');
fprintf(fileID,'        };\n');
fprintf(fileID,'        printf "shiftn" >> "OUTPUT/Ten022QD/FaceEdge.txt";\n');
fprintf(fileID,'        printf "shiftn" >> "OUTPUT/Ten022QD/EdgeInfo.txt";\n');
fprintf(fileID,'        printf "shiftn" >> "OUTPUT/QD/FaceEdge_Ten022QD.txt";\n');
fprintf(fileID,'        printf "shiftn" >> "OUTPUT/QD/EdgeInfo_Ten022QD.txt";\n');
fprintf(fileID,'     };\n');
fprintf(fileID,'     \n');
fprintf(fileID,'     for( inx:=1 ; inx<=400 ; inx:=inx+1 )\n');
fprintf(fileID,'        {\n');
fprintf(fileID,'            sumN := 0;\n');
fprintf(fileID,'            foreach facet[inx].vertex vvv where valence >= 3 do sumN := sumN + 1;\n');
fprintf(fileID,'            printf "|d |d |f shiftn", inx, sumN, body[inx].volume>> "OUTPUT/Ten022QD/Neighbor.txt";\n');
fprintf(fileID,'            printf "|d |d |f shiftn", inx, sumN, body[inx].volume>> "OUTPUT/QD/Neighbor_Ten022QD.txt";\n');
fprintf(fileID,'            \n');
fprintf(fileID,'            printf "|d ", inx >> "OUTPUT/Ten022QD/Vx.txt";\n');
fprintf(fileID,'            printf "|d ", inx >> "OUTPUT/Ten022QD/Vy.txt";\n');
fprintf(fileID,'            printf "|d ", inx >> "OUTPUT/QD/Vx_Ten022QD.txt";\n');
fprintf(fileID,'            printf "|d ", inx >> "OUTPUT/QD/Vy_Ten022QD.txt";\n');
fprintf(fileID,'            \n');
fprintf(fileID,'            foreach facet[inx].vertex vvv do {\n');
fprintf(fileID,'                printf "|f ", vvv.x >> "OUTPUT/Ten022QD/Vx.txt";\n');
fprintf(fileID,'                printf "|f ", vvv.y >> "OUTPUT/Ten022QD/Vy.txt";\n');
fprintf(fileID,'                printf "|f ", vvv.x >> "OUTPUT/QD/Vx_Ten022QD.txt";\n');
fprintf(fileID,'                printf "|f ", vvv.y >> "OUTPUT/QD/Vy_Ten022QD.txt";\n');
fprintf(fileID,'            };\n');
fprintf(fileID,'            printf " shiftn" >> "OUTPUT/Ten022QD/Vx.txt";\n');
fprintf(fileID,'            printf " shiftn" >> "OUTPUT/Ten022QD/Vy.txt";\n');
fprintf(fileID,'            printf " shiftn" >> "OUTPUT/QD/Vx_Ten022QD.txt";\n');
fprintf(fileID,'            printf " shiftn" >> "OUTPUT/QD/Vy_Ten022QD.txt";\n');
fprintf(fileID,'        };\n');
fprintf(fileID,'     \n');
fprintf(fileID,'     foreach vertex vvv do {\n');
fprintf(fileID,'        printf "|d |f |f shiftn", vvv.id, vvv.x, vvv.y >>"OUTPUT/Ten022QD/Vertex.txt";\n');
fprintf(fileID,'        printf "|d |f |f shiftn", vvv.id, vvv.x, vvv.y >>"OUTPUT/QD/Vertex_Ten022QD.txt";};\n');
fprintf(fileID,'     \n');
fprintf(fileID,'     foreach edge eee do printf "|d |f shiftn", eee.id, eee.tension >> "OUTPUT/Ten022QD/tension.txt";\n');
fprintf(fileID,'     foreach edge eee do printf "|d |f shiftn", eee.id, eee.tension >> "OUTPUT/QD/tension_Ten022QD.txt";\n');
fprintf(fileID,'     \n');
fprintf(fileID,'     foreach body bbb do printf "|d |f shiftn", bbb.id, bbb.volume >> "OUTPUT/Ten022QD/volume.txt";\n');
fprintf(fileID,'     foreach body bbb do printf "|d |f shiftn", bbb.id, bbb.volume >> "OUTPUT/QD/volume_Ten022QD.txt";\n');
fprintf(fileID,'     \n');
fprintf(fileID,'     foreach facet fff do printf "|f |f shiftn", fff.perimeter, fff.Pref >> "OUTPUT/Ten022QD/perimeter.txt";\n');
fprintf(fileID,'     foreach facet fff do printf "|f |f shiftn", fff.perimeter, fff.Pref >> "OUTPUT/QD/perimeter_Ten022QD.txt";\n');
fprintf(fileID,'}\n');
fprintf(fileID,'\n');

fprintf(fileID,'oper := {\n');
fprintf(fileID,'     autopop;\n');
fprintf(fileID,'     \n');
fprintf(fileID,'     //for(inx:=1;inx<=20;inx++)\n');
fprintf(fileID,'     //{\n');
fprintf(fileID,'     	test 10;\n');
fprintf(fileID,'     	foreach facet fff do printf "|d |d shiftn", fff.id, fff.valence >> "OUTPUT/Ten022QD/NeiIter.txt";\n');
fprintf(fileID,'        foreach facet fff do {\n');
fprintf(fileID,'            foreach fff.vertex vvv do {\n');
fprintf(fileID,'     			printf "|f ", vvv.x >> "OUTPUT/Ten022QD/VxIter.txt";\n');
fprintf(fileID,'     			printf "|f ", vvv.y >> "OUTPUT/Ten022QD/VyIter.txt";\n');
fprintf(fileID,'     		};\n');
fprintf(fileID,'     		printf " shiftn" >> "OUTPUT/Ten022QD/VxIter.txt";\n');
fprintf(fileID,'     		printf " shiftn" >> "OUTPUT/Ten022QD/VyIter.txt";\n');
fprintf(fileID,'     	};\n');
fprintf(fileID,'     	\n');
fprintf(fileID,'     	printf "|f shiftn", var >> "OUTPUT/Ten022QD/varIter.txt";\n');
fprintf(fileID,'     	printf "|f shiftn", avg(edge, log(abs(tension))/log(10)) >> "OUTPUT/Ten022QD/tenIter.txt";\n');
fprintf(fileID,'     //};\n');
fprintf(fileID,'}\n');
fprintf(fileID,'\n');

fprintf(fileID,'annealQD := {\n');
fprintf(fileID,'     iniqd;gg;\n');
fprintf(fileID,'     mash := 0.076; Nit := 7;\n');
fprintf(fileID,'     for(inx:=1;inx<=Nit;inx++){ mash := mash/2; m mash; autopop; oper 20*2^(inx-1);};\n');
fprintf(fileID,'}\n');

fprintf(fileID,'\n');

fprintf(fileID,'sche := {\n');
fprintf(fileID,'     iniqd;gg;\n');
fprintf(fileID,'     m 0.076; autopop; test 500;//oper 80;\n');
fprintf(fileID,'     // m 0.02; oper 100;\n');
fprintf(fileID,'     m 0.0006; autopop; test 20000;//oper 10000;// m 0.02; oper 100;\n');
fprintf(fileID,'}\n');

% fprintf(fileID,'record := {\n');
% fprintf(fileID,['	for( inx:=1 ; inx<=',num2str(length(Cell)),' ; inx:=inx+1 )\n']);
% fprintf(fileID,'        {\n');
% fprintf(fileID,'        sumN := 0;\n');
% fprintf(fileID,'        foreach facet[inx].vertex vvv where valence == 3 do sumN := sumN + 1;\n');
% fprintf(fileID,'        printf ;\n');
% fprintf(fileID,'        };\n');

% fprintf(fileID,' "%d %d %f\n", inx, sumN, body[inx].volume>> "C:/Users/Ran Li/Documents/Evolver/Periodic/Neighbor.txt" ')


% for i = 1:length(Face) % % X-Moment
% %     if Phase(i) == 2
%         fprintf(fileID,['    printf "%f %f %f\n", xmoment',num2str(i),'.value, ymoment',num2str(i),'.value, body[',num2str(i),'].volume >> "C:/Users/Ran Li/Documents/Evolver/Periodic/Eccentricity.txt";\n ']);  %xmoment',num2str(i),'};\n'
% %     end
% end
% fprintf(fileID,'	\n');
% fprintf(fileID,'}\n');