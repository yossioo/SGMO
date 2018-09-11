function [W_CH,W] = W_CH_from_Contacts(Contacts, origin_point,moment_normalization)
%W_CH_FROM_CONTACTS Generates an array of wrench space vectors
%   Given set of contacts, and origin point, a convex hull of wrench
%   space vectors is generated.
if nargin < 2
    origin_point = [0 0];
    moment_normalization = 1;
end
if nargin <3
    moment_normalization = 1;
end

N_C = numel(Contacts);
W = zeros(3,N_C);
if iscell(Contacts)
    for c_i = 1:N_C
        W(:,c_i) = [Contacts{c_i}.direction_vector(:);
            cross2d(Contacts{c_i}.point_on_the_line-origin_point,...
            Contacts{c_i}.direction_vector(:))/...
            moment_normalization];
    end
elseif strcmp(class(Contacts),'ContactVector')
    for c_i = 1:N_C
        W(:,c_i) = [Contacts(c_i).direction_vector(:);
            cross2d(Contacts(c_i).point_on_the_line-origin_point,...
            Contacts(c_i).direction_vector(:))/...
            moment_normalization];
    end
end
try
    K = convhulln(W');
    
    W_CH = W(:,unique(K));
catch
    W_CH = W;
end
end

