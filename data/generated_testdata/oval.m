num_regular_points = 20;
num_rand_points = 6;

major_axis = 20;
minor_axis = 10;

theta_regular = linspace(0, 2 * pi, num_regular_points);

theta_rand = 2 * pi * rand(1, num_rand_points);

theta = sort([theta_regular, theta_rand]);

points = zeros(length(theta), 2);

for t = 1:length(theta)
    points(t, 1) = major_axis * cos(theta(t));
    points(t, 2) = minor_axis * sin(theta(t));
end

