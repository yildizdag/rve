%MicroStructure Generator:

% Volume Fraction of Carbon Black:
v_cb = 0.02;
% Volume Fraction of Carbon Nanotubes:
v_cnt = 0.02;
% Volume fraction of Matrix:
v_m = 1 - v_cb - v_cnt;
% Particle Radius - Carbon Black:
r_cb = 25E-9;
% Length - Carbon Nanotube:
l_cnt = 1E-6;
% Radius - Carbon Nanotube:
r_cnt = 5E-9;
% Domain sizes per each Sample:
L = 700E-9;
n_cb = ceil((v_cb*L^3)/((4/3)*pi*r_cb^3)); 
% Volume of CNT:
n_cnt = floor((v_cnt*(L^3))/(pi*(r_cnt^2)*l_cnt));
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%GENERATE         %%%%%%%%%%%
%%%%%%%MICROSTRUCTURE       %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%
%Total Surface Area:
totalSurface = n_cb*4*pi*r_cb^2 + n_cnt*(2*pi*r_cnt)*l_cnt;

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% MESH %%%%%%%%%%%
el = 10;
h1 = L/el;
h2 = L/el;
h3 = L/el;
size1 = 0:h1:L;
size2 = 0:h2:L;
size3 = 0:h3:L;
nel = el*el*el;
[X1,Y1,Z1] = meshgrid(size1,size2,size3);
    
nnp1 = length(size1);
nnp2 = length(size2);
nnp3 = length(size3);
%Total #Boundary Nodes:
bNodesNum = 2*nnp1*nnp2+2*nnp1*nnp3+2*nnp2*nnp3;
nodes = zeros(nnp1*nnp2*nnp3,3);
count = 1;
for k = 1:length(size3)
    for i = 1:length(size2)
        for j = 1:length(size1)
            nodes(count,1) = X1(i,j,k);
            nodes(count,2) = Y1(i,j,k);
            nodes(count,3) = Z1(i,j,k);
            count = count + 1;
        end
    end
end
%Connectivity:
conn = zeros(nel,8);
for j = 1:el
    for k = 1:el
        for i = 1:el
            conn(i+el*(k-1)+(el*el)*(j-1),:) = [i+(el+1)*(k-1)+(nnp1*nnp2)*(j-1), i+(el+1)*(k-1)+(nnp1*nnp2)*(j-1)+1, i+(el+1)*(k)+(nnp1*nnp2)*(j-1)+1, i+(el+1)*(k)+(nnp1*nnp2)*(j-1),...
                                                                                i+(el+1)*(k-1)+(nnp1*nnp2)*(j), i+(el+1)*(k-1)+(nnp1*nnp2)*(j)+1, i+(el+1)*(k)+(nnp1*nnp2)*(j)+1, i+(el+1)*(k)+(nnp1*nnp2)*(j)];
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%  MWCNT  %%%%%%%%%%

x_tube = zeros(n_cnt,2);
y_tube = zeros(n_cnt,2);
z_tube = zeros(n_cnt,2);
while i < n_cnt
    for i = 1:n_cnt
        check = 0;
        count = 0;
        disp(i)
        while check == 0
            count = count + 1;
            x_tube(i,:) = L.*rand(1,2);
            y_tube(i,:) = L.*rand(1,2);
            z_tube(i,:) = L.*rand(1,2);
            dist = sqrt((x_tube(2)-x_tube(1))^2 + (y_tube(2)-y_tube(1))^2 + (z_tube(2)-z_tube(1))^2);
            if dist>=0.95*l_cnt && dist <= 1.05*l_cnt
                if i ~= 1
                    dist = zeros(1,i-1);
                    for j = 1:length(dist)
                        distSegment = DistBetween2Segment([x_tube(i,1), y_tube(i,1), z_tube(i,1)], [x_tube(i,2) y_tube(i,2) z_tube(i,2)],...
                                                                                                     [x_tube(i-j,1) y_tube(i-j,1) z_tube(i-j,1)], [x_tube(i-j,2) y_tube(i-j,2) z_tube(i-j,2)]);
                        dist(i-j) = distSegment;
                        rel = find(dist<0.95*2*r_cnt);
                        if isempty(rel) == 0
                            check = 0;
                        elseif isempty(rel) == 1
                            check = 1;
                        end
                    end
                else
                    check = 1;
                end
            else
                check = 0;
            end
            if count == 20000
                break
            end
        end
        if count == 20000
            i=0;
            break
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% CARBON BLACK  %%%%%%%

x_p = zeros(n_cb,1);
y_p = zeros(n_cb,1);
z_p = zeros(n_cb,1);
i = 0;
while i < n_cb
    for i = 1:n_cb
        check = 0;
        count = 0;
        disp(i);
        while check == 0
            count = count + 1;
            x_p(i) = L .* rand(1);
            y_p(i) = L .* rand(1);
            z_p(i) = L .*rand(1);
            %Check Tubes:
            distTubes = zeros(1,n_cnt);
            for kk = 1:n_cnt
                [distT,~,~] = distancePoint2Line([x_tube(kk,1),y_tube(kk,1),z_tube(kk,1)],[x_tube(kk,2),y_tube(kk,2),z_tube(kk,2)],[x_p(i),y_p(i),z_p(i)],'segment');
                distTubes(kk) = distT;
            end
            if n_cnt  == 0
                distT = [];
            end
            relT = find(distT<(0.95*(r_cnt+r_cb)));
            if isempty(relT) == 0
                break;
            elseif isempty(relT) == 1
            end

            if (0 < x_p(i)-r_cb) && (x_p(i)+r_cb < L) && (0 < y_p(i)-r_cb) && (y_p(i)+r_cb < L) && (0 < z_p(i)-r_cb) && (z_p(i)+r_cb < L)
                if i ~= 1
                    dist = zeros(1,i-1);
                    for j = 1:length(dist)
                        dist(i-j) = sqrt((x_p(i)-x_p(i-j))^2 + (y_p(i)-y_p(i-j))^2 + (z_p(i)-z_p(i-j))^2);
                        rel = find(0<dist & dist<0.95*2*r_cb);
                        if isempty(rel) == 0
                            check = 0;
                        elseif isempty(rel) == 1
                            check = 1;
                        end
                    end
                else
                    check = 1;
                end
            else
                check = 0;
            end
            if count == 20000
                break
            end
        end
        if count == 20000
            i=0;
            break
        end 
    end
end
pos_cb = [x_p, y_p, z_p];
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

[x, y, z] = sphere;
figure;
hold on
for i = 1:n_cb
    surf(x*r_cb+x_p(i),y*r_cb+y_p(i),z*r_cb+z_p(i),'EdgeColor','none','FaceColor',[0.6 0 0]);
end

for i = 1:n_cnt
    [X, Y, Z] = cylinder2P(r_cnt,100,[x_tube(i,1),y_tube(i,1),z_tube(i,1)],[x_tube(i,2),y_tube(i,2),z_tube(i,2)]);
    surf(X, Y, Z,'LineStyle','none','EdgeColor','none','FaceLighting','none','FaceColor',[0 0 0.6]);
    %plot3(x_tube(i,:),y_tube(i,:),z_tube(i,:),'LineWidth',4);
end

for i = 1:length(conn)
    vert = nodes(conn(i,:),:);
    fac = [1 2 6 5;2 3 7 6;3 4 8 7;4 1 5 8;1 2 3 4;5 6 7 8];
    patch('Vertices',vert,'Faces',fac,'EdgeColor','k','FaceColor','none','LineWidth',2);
end
hold off
axis equal
view(30,20)
axis off

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%INITIAL CONTACT %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%
%Carbon Black
%%%%%%%%%%%%%%%%%%%%%%%%%
cAreaInit = 0;
for i = 1:n_cb
    for j = 1:n_cb
        if j > i
            dist = sqrt((pos_cb(j,1)-pos_cb(i,1))^2+(pos_cb(j,2)-pos_cb(i,2))^2+(pos_cb(j,3)-pos_cb(i,3))^2);
            if dist <= (1*2*r_cb)
                cAreaInit = cAreaInit + (4*pi*r_cb^2);
            end
        end
    end
    for j = 1:n_cnt
        [distT,~,~] = distancePoint2Line([x_tube(j,1),y_tube(j,1),z_tube(j,1)],[x_tube(j,2),y_tube(j,2),z_tube(j,2)],[x_p(i),y_p(i),z_p(i)],'segment');
        if distT <= (1*(r_cnt+r_cb))
            cAreaInit = cAreaInit + (4*pi*r_cb^2);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%MWCNT
%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:n_cnt
    for j = 1:n_cnt
        if j > i
            distSegment = DistBetween2Segment([x_tube(i,1), y_tube(i,1), z_tube(i,1)], [x_tube(i,2) y_tube(i,2) z_tube(i,2)],...
                                                                                                     [x_tube(j,1) y_tube(j,1) z_tube(j,1)], [x_tube(j,2) y_tube(j,2) z_tube(j,2)]);
            if dist <= 1*2*r_cnt
                %cAreaInit = cAreaInit + 2*pi*r_cnt*l_cnt;
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%DEFORM%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
d_grad = [1 0 0; 0 0 0; 0 0 0];
H = eye(3)+d_grad;
dnodes = zeros(size(nodes,1),3);
dpos_cb = zeros(n_cb,3);
dx_tube = zeros(n_cnt,2);
dy_tube = zeros(n_cnt,2);
dz_tube = zeros(n_cnt,2);
for i = 1:size(nodes,1)
    dnodes(i,:) = transpose(H*transpose(nodes(i,:)));
end
for i = 1:n_cb
    dpos_cb(i,:) = transpose(H*transpose(pos_cb(i,:)));
end
for i = 1:n_cnt
    mid_tube = [0.5*(x_tube(i,1)+x_tube(i,2)); 0.5*(y_tube(i,1)+y_tube(i,2)); 0.5*(z_tube(i,1)+z_tube(i,2))];
    dmid_tube = (H*mid_tube);
    umid = dmid_tube - mid_tube;
    dx_tube(i,1) = x_tube(i,1)+umid(1); dy_tube(i,1) = y_tube(i,1)+umid(2); dz_tube(i,1) = z_tube(i,1)+umid(3);
    dx_tube(i,2) = x_tube(i,2)+umid(1); dy_tube(i,2) = y_tube(i,2)+umid(2); dz_tube(i,2) = z_tube(i,2)+umid(3);
end
cAreaDef = 0;
for i = 1:n_cb
    for j = 1:n_cb
        if j > i
            dist = sqrt((dpos_cb(j,1)-dpos_cb(i,1))^2+(dpos_cb(j,2)-dpos_cb(i,2))^2+(dpos_cb(j,3)-dpos_cb(i,3))^2);
            if dist <= (1*2*r_cb)
                cAreaDef = cAreaDef + (4*pi*r_cb^2);
            end
        end
    end
   	for j = 1:n_cnt
    	[distT,~,~] = distancePoint2Line([dx_tube(j,1),dy_tube(j,1),dz_tube(j,1)],[dx_tube(j,2),dy_tube(j,2),dz_tube(j,2)],[dpos_cb(i,1),dpos_cb(i,2),dpos_cb(i,3)],'segment');
        if distT <= (1*(r_cnt+r_cb))
            cAreaDef = cAreaDef + (4*pi*r_cb^2);
        end
    end
end
for i = 1:n_cnt
    for j = 1:n_cnt
        if j > i
            distSegment = DistBetween2Segment([dx_tube(i,1), dy_tube(i,1), dz_tube(i,1)], [dx_tube(i,2) dy_tube(i,2) dz_tube(i,2)],...
                                                                                                     [dx_tube(j,1) dy_tube(j,1) dz_tube(j,1)], [dx_tube(j,2) dy_tube(j,2) dz_tube(j,2)]);
            if dist <= 1*2*r_cnt
                %cAreaDef = cAreaDef + 2*pi*r_cnt*l_cnt;
            end
        end
    end
end
cLost = 100*(cAreaDef/cAreaInit);

[x, y, z] = sphere;
figure;
hold on
for i = 1:n_cb
    surf(x*r_cb+dpos_cb(i,1),y*r_cb+dpos_cb(i,2),z*r_cb+dpos_cb(i,3),'EdgeColor','none','FaceColor',[0.6 0 0]);
end

for i = 1:n_cnt
    [X, Y, Z] = cylinder2P(r_cnt,100,[dx_tube(i,1),dy_tube(i,1),dz_tube(i,1)],[dx_tube(i,2),dy_tube(i,2),dz_tube(i,2)]);
    surf(X, Y, Z,'LineStyle','none','EdgeColor','none','FaceLighting','none','FaceColor',[0 0 0.6]);
    %plot3(x_tube(i,:),y_tube(i,:),z_tube(i,:),'LineWidth',4);
end

for i = 1:length(conn)
    vert = dnodes(conn(i,:),:);
    fac = [1 2 6 5;2 3 7 6;3 4 8 7;4 1 5 8;1 2 3 4;5 6 7 8];
    patch('Vertices',vert,'Faces',fac,'EdgeColor','k','FaceColor','none','LineWidth',2);
end
hold off
axis equal
view(30,20)
axis off