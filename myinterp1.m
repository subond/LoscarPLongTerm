%% Interpolate with no extrapolation
%% The values outside of the range are clipped at their minimum and maximum
%% values

function vq = myinterp1(x,v,xq)

    vq = interp1(x,v,xq);

    [XMax, idxVMax] = max(x);
    [XMin, idxVMin] = min(x);

    idxMax = xq > XMax;
    idxMin = xq < XMin;

    vq(idxMax) = v(idxVMax);
    vq(idxMin) = v(idxVMin);

end