 % Define constants
dP_dx = -8;
rho = 1;
mu = 1;
h = 1;
delta_t = 0.0001;
delta_y = 0.001;
points_y = ceil(h / delta_y) + 1;

% Calculate lambda
lambda = (mu * delta_t) / (delta_y^2);

% Initialize arrays
u = zeros(points_y, 1);
A = zeros(points_y, points_y);
B_1 = zeros(points_y, 1);
B_2 = zeros(points_y, 1);

% Create Left Hand Side Matrix (A)
for i = 1:(points_y - 1)
  A(i, i - 1) = -lambda / 2;
  A(i, i) = 1 + lambda;
  A(i, i + 1) = -lambda / 2;
end

A(1, 1) = 1;
A(points_y, points_y) = 1;

% Create First Right Hand Side Matrix (B_1)
for i = 1:(points_y - 1)
  B_1(i, 1) = -delta_t / rho * dP_dx;
end

% Create second Right Hand Side Matrix (B_2)
for i = 1:(points_y - 2)
  B_2(i, 1) = lambda / 2 * u(i - 1) + (1 - lambda) * u(i) + lambda / 2 * u(i + 1);
end

% Create automatic controls for when the lines have to be plotted
upper_limit_time = [];
times_list = [0.025, 0.05, 0.1, 0.4];
for i = 1:length(times_list)
  upper_limit_time = [upper_limit_time, int(times_list(i) / delta_t)];
end

% Calculate u_n trough sparse matrix method
list_of_u_n = [];
for x = 0:int(0.4/delta_t)
    % sum B to u_n
    RHS_sum = B_1 + B_2;

    % Matrix multiplication between A and RHS_sum
    u_n = A \ RHS_sum;

    % Recreate second Right Hand Side Matrix in the sparse matrix method
    for i = 1:(points_y - 2)
        B_2(i, 1) = lambda / 2 * u_n(i - 1) + (1 - lambda) * u_n(i) + lambda / 2 * u_n(i + 1);
    end

    for t = upper_limit_time
        if x == t
            plot(u_n, 1:points_y);
            list_of_u_n = [list_of_u_n; u_n.'];
        end
    end
end

% Plot
times_list_str = string(times_list);
for i = 1:length(times_list)
  plot(list_of_u_n(:, i), 1:points_y, '--', 'DisplayName',['Analytical t=', times_list_str(i)]);
end
for i = 1:length(times_list)
  plot(list_of_u_n(:, i + length(times_list)), 1:points_y, 'DisplayName',['Numerical t=', times_list_str(i)]);
end
legend, xlabel('u'), ylabel('y'), title('Poiseuille Flow in the x-direction'), savefig('Matrix attempt.png'), hold off;

% Calculate B_n elements
B_n = @(n) 8 * ((pi * n)*sin(pi*n) + 2 * cos(pi * n) - 2) / (pi ^ 3 * n ^ 3);

% Calculate when B_n < 1e-8
for n = 1:500
    if n % 2 != 0
        if abs(B_n(n))< 1e-8
            n_max = n;
            break;
        end
    end
end

exponential_factor = @(n, t) exp(-n^2*pi^2*mu*t/(h^2));
