function index = cyclicIndex(axis_index,offset)
% Return a wrapped component index in the three-dimensional vector.
index = mod(axis_index + offset - 1,3) + 1;