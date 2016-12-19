sizeage = [261 135  87 115 158];
pbin = [0.1379 0.1185 0.0919 0.2521 0.2468];
plblist = [];
publist = [];
for i=1:length(agesize)
    [ plb pub ] = calError( sizeage(i), pbin(i) )
    plblist(i) = plb;
    publist(i) = pub; 
end
