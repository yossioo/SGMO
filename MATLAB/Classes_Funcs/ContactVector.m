classdef ContactVector
    %CONTACTVECTOR plucker vector
    %   Detailed explanation goes here
    
    properties
        point_on_the_line = [0 0];
        direction_vector = [1 0];
        length = 1;
        default_color = [1 0 0];
        polygon_num
        acting_polygon_num
    end
    
    methods
        function obj = ContactVector(point, direction, length, polygon_num, acting_polygon_num)
            %VECTOR Construct an instance of this class
            %   Detailed explanation goes here
            switch nargin
                case 1
                    obj.point_on_the_line = reshape(point(1:2),1,2);
                    obj.direction_vector = reshape(point(3:4),1,2);
                    obj.length = 1;
                case 2
                    obj.point_on_the_line = point(:).';
                    obj.direction_vector = direction(:).' ./ norm(direction);
                    obj.length = 1;
                case 3
                    obj.length = length;
                    obj.point_on_the_line = point(:).';
                    obj.direction_vector = direction(:).' ./ norm(direction);
                case 4
                    obj.length = length;
                    obj.point_on_the_line = point(:).';
                    obj.direction_vector = direction(:).' ./ norm(direction);
                    obj.polygon_num = polygon_num;
                case 5
                    obj.length = length;
                    obj.point_on_the_line = point(:).';
                    obj.direction_vector = direction(:).' ./ norm(direction);
                    obj.polygon_num = polygon_num;
                    obj.acting_polygon_num = acting_polygon_num;
            end
        end
        
        function CCW = cross_around_point(obj, point)
            % Returns the moment of the vector aroung a point
            r = obj.point_on_the_line - point;
            CCW = r(1)*obj.direction_vector(2) - r(2)*obj.direction_vector(1);
%             CCW = CCW./norm(CCW);
        end
        
        function line_segment = inf_line(obj)
            % Returns line at distance of 2e4
            line_segment = obj.point_on_the_line + 1e4*[-obj.direction_vector; obj.direction_vector];
        end
        
        function line_segment = semi_inf_line(obj)
            % Returns line at distance of 2e4
            line_segment = obj.point_on_the_line + 1e4*[0 0; obj.direction_vector];
        end
        
        function plot_contact(obj, color)
            if nargin == 1
                color = [1 0 0];
            end
            pA = obj.point_on_the_line-obj.direction_vector*obj.length;
            q = quiver(pA(1),pA(2),...
                obj.direction_vector(1)*obj.length, ...
                obj.direction_vector(2)*obj.length,...
                0);
            q.MaxHeadSize = .5;
            q.Color = color;
            q.LineWidth = 1;
        end
        
        function p = get_finger_center(obj, finger_diameter)
            p = obj.point_on_the_line - finger_diameter/2*obj.direction_vector;
        end
        
        function draw_inf_line(obj, color)
            if nargin == 1
                color = [.3 .3 .3];
            end
            points = obj.inf_line();
            l = line(points(:,1),points(:,2));
            l.Color = color;
            l.LineStyle = '--';
        end
        
        function draw(obj, color)
            if nargin == 1
                color = obj.default_color;
            end
            warning("draw is deprecated, use plot instead.")
            obj.plot(color)
        end
        
        function plot(obj, color)
            if nargin == 1
                color = obj.default_color;
            end
            points = obj.point_on_the_line + obj.length*[[0 0]; obj.direction_vector];
            thickness = 2;
            head_size = 3;
            pA = points(1,:);
            pB = points(2,:);
            v = pB-pA;
            quiver(pA(1),pA(2),v(1),v(2),...
                'LineWidth',thickness,...
                'Color',color,...
                'AutoScale','off',...
                'MaxHeadSize',head_size);
        end
    end
end

