classdef Polygon_mkII < handle
    %POLYGON_MKII Summary of this class goes here
    %   Detailed explanation goes here
    % The imporved polygon class will consist of a polygonal object along
    % with several additional methods:
    % plot - plotting the polygon
    
    
    properties
        Name; % Can save some ID
        Shape
        Polygon_color = [1 0 0];
        Edges; % [N+1x2] matrix with edges points, clockwise
        Inner_normals;
        Variants2;
        Variants3;
        E_lengths;
        E_norm_lengths
        N_e;
        Area;
        Cum_length;  % Cumulative length of all the perimeter
        Norm_cum_len; % Normalized Cumulative Length
        Center;
        Turn_directions; %1 for left; -1 for right
        Concave_ind; % Indices of concave vertices
        V_angles; % Vertex inner angles
    end
    
    methods
        function obj = Polygon_mkII(edges,name, color)
            %POLYGON_MKII Construct an instance of this class
            %   edges - [Nx2] matrix with [X,Y] coordinates of the polygon
            %   edges.
            if nargin >1 
                obj.Name = name;
            else
                name = randi([1 52],1,5);
                name(name > 26) = name(name > 26)+6;
                obj.Name = char(name+64);
            end
            if nargin > 2
                obj.Polygon_color = color; 
            end
            
            if ~isequal(class(edges),'polyshape')
                obj.Shape = polyshape(edges);
            else
                obj.Shape = edges;
            end
            update_from_shape(obj)
        end
        
        function find_concave_vertices(obj)
            v1 = obj.Edges(1,:)-obj.Edges(obj.N_e,:);
           v2 = obj.Edges(2,:)-obj.Edges(1,:);
           cr = cross2d(v1/norm(v1),v2/norm(v2));
           obj.Turn_directions(1) = sign(cr);
           obj.V_angles(1) = round(180+angle_of_2_vec(v1,v2),1);
           for i = 1:obj.N_e-1
               v1 = obj.Edges(i+1,:)-obj.Edges(i,:);
               v2 = obj.Edges(i+2,:)-obj.Edges(i+1,:);
               cr = cross2d(v1/norm(v1),v2/norm(v2));
               obj.Turn_directions(i+1) = sign(cr);
               obj.V_angles(i+1) = round(180+angle_of_2_vec(v1,v2),1);
           end
           obj.Concave_ind = obj.Turn_directions > 0;
        end
        
        function rotate(obj, theta_d, point)
            if nargin < 3
                point = obj.Center;
            end
            obj.Shape = rotate(obj.Shape, theta_d, point);
            update_from_shape(obj);
        end
                
        function vector = find_vector_from_vertex(obj, i, forward)
            % Given a vertex id, returns a vector along the boundary.
            % If forward is positive returns clockwise vector;
            % otherwise returns counterclockwise vector;
            % backward
            if forward > 0 
                vector = diff(obj.Edges([i,i+1],:));
            else
                if i == 1
                    vector = diff(obj.Edges([end,end-1],:));
                else
                    vector = diff(obj.Edges([i,i-1],:));
                end
            end
            vector = round(vector/norm(vector),5);
        end
        
        
        function update_from_shape(obj)
            obj.Shape.Vertices = round(obj.Shape.Vertices,5);
            [x,y] = boundary(obj.Shape);
            obj.N_e = length(x)-1;
            
            obj.Edges = [x,y];
            obj.find_concave_vertices();
            obj.calculate_normals();
            obj.calculate_lengths();
            obj.Variants2 = nchoosek(1:length(obj.Shape.Vertices),2);
            obj.Variants3 = nchoosek(1:length(obj.Shape.Vertices),3);
            
            obj.Area = obj.Shape.area();
            if obj.N_e > 0
                obj.Cum_length = zeros(obj.N_e, 1);
                for i  = 1:obj.N_e
                    obj.Cum_length(i) = sum(obj.E_lengths(1:i));
                end
                obj.Norm_cum_len = ...
                    obj.Cum_length ./ obj.Cum_length(obj.N_e);
                obj.Norm_cum_len = ...
                    [0; obj.Norm_cum_len(:)];
                obj.Center = mean(obj.Shape.Vertices(~isnan(obj.Shape.Vertices(:,1)),: ));
            end
            obj.E_norm_lengths = obj.E_lengths/obj.Cum_length(end);
        end
        
        function translate(obj, diff_vec)
            new_edges = obj.Edges(1:end-1,:) + diff_vec;
            obj.Shape = polyshape(new_edges);
            update_from_shape(obj)
        end
        
        function inter = find_intersection_normals(obj,indices)
            l = 10; % length of the corridor
            a = indices(1);
            b = indices(2);
            edge_A = obj.Edges([a a+1],:);
            normal_A = obj.Inner_normals(a,:);
            edge_B = obj.Edges([b b+1],:);
            normal_B = obj.Inner_normals(b,:);
            shape_A = polyshape([flipud(edge_A-l*normal_A); edge_A+l*normal_A]);
            shape_B = polyshape([flipud(edge_B-l*normal_B); edge_B+l*normal_B]);
            inter = intersect(shape_A, shape_B);
        end
        
        function calculate_normals(obj)
            E = obj.Edges;
            obj.Inner_normals = zeros(length(E)-1,2);
            for i  = 1:length(E)-1
                p1 = E(i,:);
                p2 = E(i+1,:);
                p = (p2-p1)/norm(p2-p1);
                obj.Inner_normals(i,:) = ([0 1; -1 0]*p(:)).';
            end
        end
        
        function calculate_lengths(obj)
            E = obj.Edges;
            obj.E_lengths = zeros(length(E)-1,1);
            for i  = 1:length(E)-1
                p1 = E(i,:);
                p2 = E(i+1,:);
                p = p2-p1;
                obj.E_lengths(i) = norm(p);
            end
        end
        
        function calculate_azimuths(obj)
            for i  = 1:length(obj.Inner_normals)-1
                n = obj.Inner_normals(i,:);
                obj.Edges_azimuth(i) = atan2(n(2),n(1));
            end
        end
        
        function plot(obj, color,show_text)
            if nargin <2
                color = obj.Polygon_color;
                show_text = 1;
            end
            
            plot(obj.Shape, 'FaceColor',color,'FaceAlpha',.4);
            if show_text
                text(obj.Center(1),obj.Center(2),obj.Name);
            end
        end
        
        function plot_vertex(obj, vertex_id, color)
            if nargin < 2
                error("No vertex index specified.")
            else
                if nargin < 3
                    color = obj.Polygon_color;
                end
                plot(obj.Edges(vertex_id,1),...
                    obj.Edges(vertex_id,2),'-s',...
                    'Color',color,'MarkerSize',15);
            end
        end
        
        function p = find_fartherst_point_from_x(obj,x)
            dist = zeros(size(obj.Shape.Vertices));
            for i = 1:length(obj.Shape.Vertices)
                vertex = obj.Shape.Vertices(i,:);
                dist(i) = norm(vertex-x);
            end
            ind = find(dist == max(dist));
            p = obj.Shape.Vertices(ind(1),:);
        end
        
        function [inside,on] = in_polygon(obj,points)
            % points have to be Nx2 array
            [inside,on] = inpolygon(points(:,1),points(:,2),...
                obj.Edges(:,1),obj.Edges(:,2));
        end
        
        function p_out = project_p_on_edge(obj, p, edge_num)
            edge = obj.Edges([edge_num edge_num+1],:);
            p_out = project_point_on_line(p, edge);
            
        end
        
        function point = point_from_edgePosition(obj,edge_num, position_on_edge)
            %             position = obj.Norm_cum_len(edge_num)+obj.E_norm_lengths(edge_num)*position_on_edge;
            if isempty(position_on_edge)
                point = [];
            else
                if edge_num > obj.N_e
                    error("Edge number exceeds number of edges")
                end
                v = obj.Edges([edge_num, edge_num+1],:);
                point = v(1,:)+diff(v)*position_on_edge;
            end
        end
        
        function n = find_normal_at_point(obj, point)
            % What about case of e2e contact, looking for normal of
            % the end-vertex?
            n = [];
            dist_from_vert = vecnorm((obj.Edges - point)');
            dist_from_vert( isnan(dist_from_vert)) = Inf;
            ind = find(dist_from_vert<sqrt(obj.Area)*1e-3);
            if ind
                warning("Polygon_mkII:VertexNormalUndefined","Vertex detected, skipping the point")
                return
            end
            for i = 1:obj.N_e
                is_on_edge = dot(point-obj.Edges(i,:),point-obj.Edges(i+1,:)) <1e-5;
                % Check whether points are collinear
                if is_on_edge
                    D  = det([point 1; obj.Edges([i i+1],:) [1;1]]);
                    if abs(D) < obj.Area *1e-4
                        ni = obj.Inner_normals(i,:);
                        ni = ni./norm(ni);
                        n = [n;ni];
                    end
                end
            end
        end
        
        function n_vecs = find_contacts_for_positions(obj, positions)
            % positions specify locations along the polygon boundary [0..1]
            n_vecs = zeros(numel(positions),4);
            tol = 1e-5;
            
            for p_i = 1:numel(positions)
                for i = 1:obj.N_e
                    if obj.Norm_cum_len(i+1) > positions(p_i)
                        break
                    end
                end
                p1 = obj.Edges(i,:);
                p2 = obj.Edges(i+1,:);
                v = p2-p1;
                v = v./norm(v);
                %                 points(p_i,:) = p1 + ...
                %                     v * (positions(p_i) - ...
                %                     obj.Norm_cum_len(i))*obj.Cum_length(N);
                %                 directions(p_i,:) = obj.Inner_normals(i,:);
                
                % Rewrite it:
                n_vecs(p_i,:) = [p1 + ...
                    v * (positions(p_i) - ...
                    obj.Norm_cum_len(i))*obj.Cum_length(obj.N_e), obj.Inner_normals(i,:)];
            end
%             n_vecs = [points, directions];
            n_vecs(abs(n_vecs) < tol) = 0;
        end
        
        function pos = point_position(obj,points)
            for point_ind = 1:size(points,1)
                point = points(point_ind,:);
                for i = 1:obj.N_e
                    proj_p = obj.project_p_on_edge(point,i);
                    if abs(point - proj_p) < obj.Area/1e3
                        break
                    end
                end
                % Now i is the edge number for given point
                p1 = obj.Edges(i,:);
                p2 = obj.Edges(i+1,:);
                v = p2-p1;
                v = v./norm(v);
                pos(point_ind) = dot((proj_p-p1), v)/obj.Cum_length(end) + obj.Norm_cum_len(i);
            end
        end
    end
end

