%% Contacts search
%#ok<*SAGROW>
% clc;
clearvars -except DEBUG PolyList finger_d P finger_d Filtered_Contacts_poly_ind Filtered_Contacts
warning off Polygon_mkII:VertexNormalUndefined


d_scale = 1e-3*2.5;
if numel(PolyList) == 1
    return
end
%% Iterate over each polygon

N = numel(PolyList);
ALL = 1:N;
C = {};
C_poly_inds = [];
for p_i = ALL
    % Remove these 2 lines later. Now want to see only for given
    % polygon
    % C = {};
    % C_poly_inds = [];
    %
    others = ALL ~= p_i;
    Unified = Polygon_mkII(get_unified_poly(PolyList(others)));
    Unified.Name = "w/o " + num2str(p_i);
    p = PolyList{p_i};
    
    % Have to search separately for v2e contacts, e2e end contacts,
    % and v2iv contacts
    
    
    % Search for common edge segments
    d = d_scale*sqrt(Unified.Area); % Bufffer measure
    for i = 1:p.N_e
        e = p.Edges([i i+1],:);
        [in,out] = intersect(polybuffer(Unified.Shape,d),e);
        if norm(diff(in)) > d*3
            n = p.find_normal_at_point(mean(in));
            if isempty(n)
                continue
            end
            for in_i = 1:size(in,1)
                C{end+1} = ContactVector(in(in_i,:),n,finger_d*2,p_i);
                C_poly_inds(end+1) = p_i;
            end
        end
    end
    
    d = d_scale*sqrt(p.Area); % Bufffer measure
    for i = 1:Unified.N_e
        e = Unified.Edges([i i+1],:);
        try
            [in,out] = intersect(polybuffer(p.Shape,d),e);
        catch
            in = [];
        end
        if norm(diff(in)) > d*3
            n = p.find_normal_at_point(mean(in));
            if isempty(n)
                continue
            end
            for in_i = 1:size(in,1)
                C{end+1} = ContactVector(in(in_i,:),n,finger_d*2,p_i);
                C_poly_inds(end+1) = p_i;
            end
        end
    end
    % Search vertices of this polygon touch another
    [~,on] = Unified.in_polygon(p.Edges(1:end-1,:));
    points = p.Edges(on,:);
    for pt_i = 1:size(points,1)
        point = points(pt_i,:);
        n = -Unified.find_normal_at_point(point);
        if isempty(n)
            continue
        end
        for n_i = 1:size(n,1)
            if n(n_i,:)
            end
            C{end+1} = ContactVector(point,n(n_i,:),finger_d*2,p_i);
            C_poly_inds(end+1) = p_i;
        end
    end
    % Search vertices of another touch this
    [~,on] = p.in_polygon(Unified.Edges(1:end-1,:));
    points = Unified.Edges(on,:);
    for pt_i = 1:size(points,1)
        point = points(pt_i,:);
        n = p.find_normal_at_point(point);
        for n_i = 1:size(n,1)
            C{end+1} = ContactVector(point,n(n_i,:),finger_d*2,p_i);
            C_poly_inds(end+1) = p_i;
        end
    end
    if DEBUG
        f = figure(15);clf  %p_i
        f.Name = "Inter-object contacts search";
        f.Name = PolyList{p_i}.Name;
        Unified.plot();
        axis equal; grid on; hold on;
        plot(polybuffer(Unified.Shape,d))
        p.plot(p.Polygon_color,false);
        for c = C(C_poly_inds == p_i)
            cont = c{1};
            cont.plot_contact()
        end
%         drawnow
        pause(0.5)
    end
end


%% Remove duplicate contacts ?

% This can be done by iterating over wrenches and `cross`-ing one with
% another. If some cross yields norm smaller than

% norm(cross([1 0 0 ],[1 0.03 0.01]))


% Iterate for each polygon?
Filtered_Contacts = {};
Filtered_Contacts_poly_ind = [];
for p_i = 1:numel(PolyList)
    p = PolyList{p_i};
    Contacts = C(C_poly_inds == p_i);
    W = zeros(numel(Contacts),3);
    for c_i = 1:numel(Contacts)
        W(c_i,:) = [ Contacts{c_i}.direction_vector,...
            cross2d(Contacts{c_i}.point_on_the_line-PolyList{Contacts{c_i}.polygon_num}.Center,...
            Contacts{c_i}.direction_vector)];
    end
    W = round(W,3);
    [~,IA,IC]  = uniquetol(W,0.03,'ByRows',true);
    Filtered_Contacts = [Filtered_Contacts(:)',   Contacts(IA)];
    Filtered_Contacts_poly_ind = [Filtered_Contacts_poly_ind(:);   0*IA+p_i];
    
end

%% Can draw filtered contacts
if DEBUG && 0
    figure(94)
    for p_i = 1:numel(PolyList) %#ok<*UNRCH>
        p = PolyList{p_i};
        Contacts = Filtered_Contacts(Filtered_Contacts_poly_ind == p_i);
        f = figure(p_i);clf
        f.Name = PolyList{p_i}.Name;
        axis equal; grid on; hold on;
        PolyList{p_i}.plot();
        
        for c_i = 1:numel(Contacts)
            Contacts{c_i}.plot_contact()
        end
        %     pause
    end
end
