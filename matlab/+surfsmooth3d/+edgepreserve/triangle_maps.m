function maps = triangle_maps(level, convention)
%TRIANGLE_MAPS Affine [origin du dv] maps into a coarse reference triangle.
validateattributes(level, {'numeric'}, {'scalar','integer','nonnegative'});
maps = zeros(2, 3, 4^level);
switch convention
    case 'step'
        n = 2^level;
        k = 0;
        for i = 0:n-1
            for j = 0:n-i-1
                k = k + 1;
                vertices = [i i+1 i; j j j+1] / n;
                maps(:,:,k) = [vertices(:,1), diff(vertices(:,1:2),1,2), ...
                    vertices(:,3)-vertices(:,1)];
                if i+j < n-1
                    k = k + 1;
                    vertices = [i+1 i+1 i; j j+1 j+1] / n;
                    maps(:,:,k) = [vertices(:,1), vertices(:,2)-vertices(:,1), ...
                        vertices(:,3)-vertices(:,1)];
                end
            end
        end
    case 'smoother'
        maps(:,:,1) = [0 1 0; 0 0 1];
        children = zeros(2,3,4);
        children(:,:,1) = [0 .5 0; 0 0 .5];
        children(:,:,2) = [.5 -.5 0; .5 0 -.5];
        children(:,:,3) = [.5 .5 0; 0 0 .5];
        children(:,:,4) = [0 .5 0; .5 0 .5];
        for depth = 1:level
            previous = maps(:,:,1:4^(depth-1));
            for parent = 1:size(previous,3)
                origin = previous(:,1,parent);
                linear = previous(:,2:3,parent);
                for child = 1:4
                    maps(:,:,(parent-1)*4+child) = ...
                        [origin + linear*children(:,1,child), ...
                         linear*children(:,2:3,child)];
                end
            end
        end
    otherwise
        error('edgepreserve:mapConvention', 'Unknown subdivision convention.');
end
end
