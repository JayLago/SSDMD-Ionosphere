function smata = varscale(smat)

    smata = mean(smat, 1);
    smata = smata - mean(smata);
    var = sqrt(mean(smata.^2));
    smata = smata/var;
        