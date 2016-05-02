function [A_map, phi_map, Z, invZ] = vfm_mappings(uvecs, n, D, epsilon)

    if (nargin == 0)
       n = 500;
       D = 2;
       uvecs = normc(normrnd(0,1,2,n*D)); 
       epsilon = 0.2;
    end
    net = epsilon_net_R2(epsilon);
    
    A_map = containers.Map;
    disc_uvecs = cell(n,1);
    for i = 1:n
        idx = (i-1)*D;
        [~, disc_uvecs{i}] = v_classify(uvecs(:,(idx+1):(idx + D))',net');
        z = num2str(disc_uvecs{i});
        if (isKey(A_map,z))
            A_map(z) = [A_map(z), i];
        else
            A_map(z) = i;
        end
    end

    ks = keys(A_map);
    phi_map = containers.Map('KeyType','double','ValueType','any');
    Z = zeros(1, length(ks));
    for i = 1:length(ks)
        z = ks{i};
        temp = A_map(z);
        Z(i) = min(temp);
        for j = temp
        phi_map(j) = Z(i);
        end
    end
    invZ = containers.Map('KeyType','double','ValueType','double');
    for i = 1:length(Z)
        invZ(Z(i)) = i;
    end

end