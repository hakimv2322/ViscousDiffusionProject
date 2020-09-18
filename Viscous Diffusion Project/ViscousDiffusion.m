% Numerically simulate viscous diffusion described
% by the viscous diffusion equation (equation 3-250)
% using an explicit numerical algorithm.
% All variables are assumed to be non-dimensionalized.

close all
clear all

h = 1/20; % mesh spacing
DeltaT = 1/4000; % time spacing
sigma = DeltaT/h^2; % mesh-size parameter

% We assume the duct width is a = 1.
b = 2; % duct height

[X,Y] = meshgrid(-1:h:1,-b:h:b); % create mesh grid
u = X; % initialize velocity u
[r, c] = size(X); % number of rows and columns

% Implement boundary conditions:
u(1,:) = 0; u(r,:) = 0; % walls at y = -b, b
u(:,1) = 0; u(:,c) = 0; % walls at x = -1, 1
% u(1,:) = 0; u(r,:) = 0; % walls at y = -b, b
% u(:,1) = 0; u(:,c) = 1; % walls at x = -1, 1

% Implement initial conditions:
for it1 = 2:r-1
    for it2 = 2:c-1
        x = X(it1,it2); y = Y(it1,it2); % extract x, y
        u(it1,it2) = u0(x,y,b);
%         u(it1,it2) = 0;
    end
end

uNew = u; % necessary for updates
T = 1; % total (non-dimensionalized) time of evolution

tArray = 0:DeltaT:T; % array of time progression
uArray = zeros(size(tArray)); % array of velocities at point of interest

% Get initial velocity profile along y = 0:
uInitial = u(round(r/2),:);


% Implement explicit numerical algorithm:
it0 = 0;
temp = 1;
for t = tArray
    it0 = it0 + 1;
    uArray(it0) = u(round(r/2), round(c/2)); % extract point of interest
    for it1 = 2:r-1
        for it2 = 2:c-1
            x = X(it1,it2); y = Y(it1,it2); % extract x, y
            A = sigma*(u(it1+1,it2) + u(it1-1,it2) + u(it1,it2+1) + u(it1,it2-1));
            B = (1 - 4*sigma)*u(it1,it2);
            uNew(it1,it2) = A + B;
        end
    end
    u = uNew;
    
    % extract evolution of velocity profile along y = 0:
    if (t > 0.25*T) && (temp == 1)
        u1 = u(round(r/2),:);
        temp = 2;
    end
    if (t > 0.5*T) && (temp == 2)
        u2 = u(round(r/2),:);
        temp = 3;
    end
    if (t > 0.75*T) && (temp == 3)
        u3 = u(round(r/2),:);
        temp = 4;
    end
end
u4 = u(round(r/2),:);



% Plot final results:

% total velocity flow after time T
surf(X,Y,u)
axis equal
xlabel('x-axis')
ylabel('y-axis')

% contour plot
figure()
contour(X,Y,u)
axis equal
xlabel('x-axis')
ylabel('y-axis')

% evolution of center velocity over time
figure()
plot(tArray, uArray, 'LineWidth',1.5)
xlabel('Time t')
ylabel('Center Velocity')

% evolution of velocity profile along y = 0
figure()
plot(X(round(r/2),:),uInitial, 'LineWidth',1.5)
hold on
plot(X(round(r/2),:),u1, 'LineWidth',1.5)
hold on
plot(X(round(r/2),:),u2, 'LineWidth',1.5)
hold on
plot(X(round(r/2),:),u3, 'LineWidth',1.5)
hold on
plot(X(round(r/2),:),u4, 'LineWidth',1.5)
legend('t = 0','t = 0.25','t = 0.5','t = 0.75','t = 1', 'Location', 'northwest')
% legend('t = 0','t = 0.25','t = 0.5','t = 0.75','t = 1', 'Location', 'southeast')
xlabel('x-location')
ylabel('Velocity')



