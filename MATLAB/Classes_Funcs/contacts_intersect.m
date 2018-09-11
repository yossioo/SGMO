function [outputArg1,outputArg2] = contacts_intersect(contact_set)
%CONTACTS_INTERSECT Summary of this function goes here
%   Detailed explanation goes here
num_contacts = numel(contact_set);
points = [contact_set.point_on_the_line]
directions = [contact_set.point_on_the_line]
Mat = [1 0 -1 0; 0 1 0 -1];
end

