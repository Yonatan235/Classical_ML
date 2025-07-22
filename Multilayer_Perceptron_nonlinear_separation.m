% This script implements the perceptron algorithm to classify a 
% semi-randomly generated 2D dataset with two separated arc-shaped classes.

clear all

% === Parameters ===
n = 200;        % Total number of data points
rad = 10;       % Radius of the circular arc for both classes
thk = 5;        % Thickness of each arc (points are scattered within a band of this width)
sep = 1;        % Vertical separation between the two arcs

% === Data Initialization ===
x = zeros(n, 2);    % Stores 2D coordinates of the data points
y = zeros(n, 1);    % Stores class labels (+1 or -1)

% Helper arrays (used for polar coordinate construction of arcs)
r1 = zeros(n,1); theta1 = zeros(n,1);
r2 = zeros(n,1); theta2 = zeros(n,1);

XX = -2.5*rad : 2.5*rad;  % Range for plotting the final decision boundary

% === Data Generation ===
for i = 1:n
    i  % Optional progress display
    
    % Random angle and radial offsets for arc-like distribution
    theta1(i) = pi * rand(1);     % Angle for class -1 (top arc), [0, π]
    theta2(i) = -pi * rand(1);    % Angle for class +1 (bottom arc), [-π, 0]
    r1(i) = thk * rand(1);        % Random thickness for top arc
    r2(i) = thk * rand(1);        % Random thickness for bottom arc
    
    if i <= n/2
        % Class -1: upper arc centered at origin
        x(i,1) = (rad + r1(i)) * cos(theta1(i));  % x-coordinate
        x(i,2) = (rad + r1(i)) * sin(theta1(i));  % y-coordinate
    else
        % Class +1: lower arc shifted slightly right and down
        x(i,1) = thk/2 + rad + (rad + r2(i)) * cos(theta2(i));
        x(i,2) = -sep + (rad + r2(i)) * sin(theta2(i));
    end

    % === Label and Plot the Data Point ===
    if i > n/2
        y(i) = 1;  % Class +1 (blue)
        c = [0 0 255]/256;
        s1 = scatter(x(i,1), x(i,2), [], c); hold on;
    else
        y(i) = -1;  % Class -1 (red)
        c = [255 0 0]/256;
        s2 = scatter(x(i,1), x(i,2), [], c); hold on;
    end

    xlabel('x1'); ylabel('x2');
end

% === Perceptron Algorithm Initialization ===
stop = 0;      % Used to stop early if no misclassified point found
w0 = 0;        % Bias term (intercept)
w1 = 0;        % Weight for x1
w2 = 0;        % Weight for x2
counter = 1;   % Count the number of iterations

% === Perceptron Learning Loop ===
% Loop until all data points are correctly classified
while sum(sign(y(:) .* (w0 + w1 * x(:,1) + w2 * x(:,2)))) ~= n
    jj = randperm(n);  % Shuffle indices to introduce randomness
    j = jj(1);         % Start with a random point
    k = 1;

    % Try to find a misclassified point
    while y(j) * (w0 + w1 * x(j,1) + w2 * x(j,2)) > 0
        if k < n
            k = k + 1;
        else
            stop = 1;  % No misclassified point found
            break
        end
        j = jj(k);  % Try next point
    end

    if stop == 1
        break
    end

    % === Update the perceptron weights ===
    w0 = w0 + y(j);               % Update bias
    w1 = w1 + y(j) * x(j,1);      % Update weight for x1
    w2 = w2 + y(j) * x(j,2);      % Update weight for x2

    % Check again if all points are correctly classified
    if sum(sign(y(:) .* (w0 + w1 * x(:,1) + w2 * x(:,2)))) == n
        i  
        break
    end

    counter = counter + 1
end

XX3 = -2.5*rad : 2.5*rad;
YY3 = (-w1 * XX3 + w0) / w2;   % Solve for x2 in w0 + w1*x1 + w2*x2 = 0

p2 = plot(XX3, YY3, 'black');  
legend([p2 s1 s2], {'Perceptron line', 'y(i)=1', 'y(i)=-1'})
