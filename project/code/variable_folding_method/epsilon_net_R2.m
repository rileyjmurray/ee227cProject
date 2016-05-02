function net = epsilon_net_R2(epsilon)

    if (nargin == 0)
        epsilon = 0.1; % default.
    end
    N = ceil(pi / acos(1-epsilon^2/2)); % the minimum number of elements in this epsilon net.
    theta = 2*pi / N; % the largest angle permissible between elements of the epsilon-net.
    net = zeros(2, N);
    net(:,1) = [1;0];
    Rot = [cos(theta) -sin(theta); sin(theta) cos(theta)];
    for i = 2:N
        net(:,i) = Rot*net(:,i-1);
    end

end