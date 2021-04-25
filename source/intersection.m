function sol = intersection(xx, yy)

    M = [(xx(2) - xx(1))/(yy(2) - yy(1)), 1; ...
         (xx(4) - xx(3))/(yy(4) - yy(3)), 1];
    b = 0.5 * [yy(1) + yy(2) + (xx(1) + xx(2)) * M(1, 1);
               yy(3) + yy(4) + (xx(3) + xx(4)) * M(2, 1)];
        
    sol = M\b;