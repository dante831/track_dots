function image = masking(x0, y0, X, Y, r, image, inner, outer)

    ind = (X - x0).^2 + (Y - y0).^2 > (r*outer)^2 | (X - x0).^2 + (Y - y0).^2 < (r*inner)^2;
    image(ind) = NaN;
end