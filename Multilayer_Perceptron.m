% This script implements the perceptron learning algorithm to classify
% a randomly generated linearly separable dataset in 2D.

% Step 1: Randomly generate the target (ground-truth) linear classifier y = ax + b
k1 = 10^2;                  % Controls range for slope (a) and intercept (b)
k2 = 10^2;                  % Controls the range for the data points
a = -(k1/2) + k1*rand(1);   % Random slope in range (-k1/2, k1/2)
b = -(k1/2) + k1*rand(1);   % Random intercept in range (-k1/2, k1/2)

% Step 2: Generate n random data points
n = 20;
x = zeros(n,2);    
y = zeros(n,1);   

% Compute the target line to plot (y = ax + b)
XX = -k2/2:k2/2;
YY = a*XX + b;

% Step 3: Assign labels to data based on the target function
for i = 1:n
    % Generate a random point in 2D space within the range [-k2/2, k2/2]
    x(i,1) = -(k2/2) + k2*rand(1); 
    x(i,2) = -(k2/2) + k2*rand(1); 

    % Assign label +1 if the point is above the target line, -1 if below
    if x(i,2) > a*x(i,1) + b
        y(i) = 1;  % Above the line
        c = [0 0 255]/256;  % Blue color for y = +1
        s1 = scatter(x(i,1), x(i,2), [], c);  
        xlabel('x1'); ylabel('x2');
        title('Classified Data Set');
        hold on
    elseif x(i,2) < a*x(i,1) + b
        y(i) = -1;  % Below the line
        c = [255 0 0]/256;  % Red color for y = -1
        s2 = scatter(x(i,1), x(i,2), [], c);  
        hold on
    end

    % Plot the target function line
    p1 = plot(XX, YY, 'g'); 
    axis([-k2/2 k2/2 -k2/2 k2/2]);
end

% Step 4: Initialize a random line (the perceptron guess)
stop = 0;  % Control flag to exit loops if needed
aa = -(k1/2) + k1*rand(1);  
bb = -(k1/2) + k1*rand(1);  

% plot initial guess line 
% XX2 = -k2/2:k2/2;
% YY2 = aa*XX2 + bb;
% plot(XX2, YY2, 'c');

% Step 5: Run the perceptron learning algorithm
counter = 1;
% Loop until all points are correctly classified
while sum(sign(y(:).*(-bb - aa.*x(:,1) + x(:,2)))) ~= n
    % Randomly permute indices for stochastic update
    jj = randperm(n);
    j = jj(1);
    k = 1;

    % Find a misclassified point
    while y(j)*(-bb - aa*x(j,1) + x(j,2)) > 0
        if k < n
            k = k + 1;
        else
            stop = 1;  % No misclassified point found
            break
        end
        j = jj(k);
    end

    if stop == 1
        break
    end

    % Update perceptron weights (aa and bb)
    bb = bb - y(j);
    aa = aa - y(j)*x(j,1);

    % Check for convergence
    if sum(sign(y(:).*(-bb - aa.*x(:,1) + x(:,2)))) == n
        i  % Print index of final iteration
        break
    end

    counter = counter + 1;
end

% Step 6: Plot the final perceptron decision boundary
XX3 = -k2/2:k2/2;
YY3 = aa*XX3 + bb;
p2 = plot(XX3, YY3, 'black');

legend([p1 p2 s1 s2], {'Target function', 'Perceptron line', 'y(i)=1', 'y(i)=-1'})
axis([-k2/2 k2/2 -k2/2 k2/2]);

