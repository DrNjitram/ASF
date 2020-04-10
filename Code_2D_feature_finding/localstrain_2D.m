function [out,D2] = localstrain_2D(a0,a1,p_size)

%calculates local strain due to particle jump in an amorphous system.
%a0 and a1 are the lists of all particle positions at times t0 and t1,


% Set the nearest neighbor distance (typically the first minimum in the pair correlation function g(r)) below


% Set the field of view here
% xmin = 2; xmax = 600;
% ymin = 2; ymax = 600;
% zmin = 2; zmax = 200;

N = size(a0,1);

% Set the range for calculating the strain here; typically a nearest neighbor distance, but can be larger to have larger average.
NNacceptance = p_size*2; % corresponds to the first minimum in the g(r)



min = 1;
epsilontensor = zeros(N,4);
D2=zeros(N,1);
for ii = 1:N
    p0 = a0(ii,1:2);     % off-atoms --------------------------
    p1 = a1(ii,1:2);
    %     if p0(1) > xmin & p0(1) < xmax & p0(2) > ymin & p0(2) < ymax & p0(3) > zmin & p0(3) < zmax
    FalkX = zeros(2,2);    % Initialize Tensors ---------------------
    FalkY = zeros(2,2);
    deform = zeros(2,2);
    R = sqrt((a0(:,1)-p0(1)).*(a0(:,1)-p0(1)) + (a0(:,2)-p0(2)) .* (a0(:,2)-p0(2)) );%+ (a0(:,3)-p0(3)) .* (a0(:,3)-p0(3)));
    XX = [a0 a1 R];    % find the neighbors --------------------
    XX(XX(:,5) > NNacceptance,:) = []; %LR original #7 ---find out why XX has 5 colums%
    XX(XX(:,5) < min,:) = [];
    nn = size(XX,1);
    X0 = XX(:,1:2);   
    X1 = XX(:,3:4);%X0 = XX(:,1:3);   X1 = XX(:,4:6);
    if nn >=3
        for i = 1:2     % get Falk's X ------------------------------
            for j = 1:2
                for n = 1:nn
                    FalkX(i,j) = FalkX(i,j) + (X1(n,i) - p1(i)) * (X0(n,j) - p0(j));
                end
            end
        end
        for i = 1:2     %get Falk's Y ----------------------------------
            for j = 1:2
                for n = 1:nn
                    FalkY(i,j) = FalkY(i,j) + (X0(n,i) - p0(i)) * (X0(n,j) - p0(j));
                end
            end
        end
        Yinv = inv(FalkY);
        for i = 1:2       % Calculate strain -------------------------------
            for j = 1:2
                for k = 1:2
                    deform(i,j) = deform(i,j) + FalkX(i,k) * Yinv(j,k);
                end
            end
        end
        deform = deform - eye(2);
        deform1 = 1/2 * (deform + deform');
        epsilontensor(ii,:) = deform1(1:4);
        for n = 1:nn     %get Non-Affine D2 ----------------------------------
            D2_i=0;
            for i = 1:2
                D2_j=0;
                for j = 1:2
                    D2_j=D2_j+((i==j)+deform(i,j))*(X0(n,j) - p0(j));
                end
                D2_i=D2_i+(X1(n,i) - p1(i)- D2_j)^2;
            end
            D2(ii) = D2(ii) + D2_i;
        end
    end
    %     end
end
out = epsilontensor;