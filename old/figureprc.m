function hcell=figureprc(fn)
    hcell=openfig([f(1:end-8) '.fig']);
    p=hcell.Position;
    set(hcell, 'Position', [p(1) p(2) p(3)*2 p(4)*2])
end