function width_centerline = width_along_centerline(x_center,y_center,x_left,y_left,width,ellipsoid)

   
    width_centerline = zeros(numel(x_center),1);

    for i = 1:numel(x_center)
        edge_dis = distance(y_center(i)*ones(numel(x_left),1),x_center(i)*ones(numel(x_left),1),y_left, x_left, ellipsoid);

        [~, edge_ind] = min(edge_dis);
        width_centerline(i) = width(edge_ind);
    end

end