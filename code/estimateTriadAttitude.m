function Cbn = estimateTriadAttitude(a,b,p,q)
% Estimate the attitude matrix from two observed and two reference vectors.

a = a(:);
b = b(:);
p = p(:);
q = q(:);
if any([numel(a),numel(b),numel(p),numel(q)] ~= 3)
    error('estimateTriadAttitude:InvalidVectorSize', ...
        'All input vectors must have three elements.');
end

body_cross = cross(a,b);
reference_cross = cross(p,q);
if norm(a) == 0 || norm(p) == 0 || norm(body_cross) == 0 ...
        || norm(reference_cross) == 0
    error('estimateTriadAttitude:DegenerateVectors', ...
        'Each vector pair must be nonzero and non-parallel.');
end

body_axis_1 = a / norm(a);
body_axis_2 = body_cross / norm(body_cross);
body_axis_3 = cross(body_axis_1,body_axis_2);
reference_axis_1 = p / norm(p);
reference_axis_2 = reference_cross / norm(reference_cross);
reference_axis_3 = cross(reference_axis_1,reference_axis_2);
Cbn = [body_axis_1,body_axis_2,body_axis_3] ...
    * [reference_axis_1,reference_axis_2,reference_axis_3]';
end
