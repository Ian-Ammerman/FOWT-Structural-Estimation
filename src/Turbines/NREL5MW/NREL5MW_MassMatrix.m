function M = NREL5MW_MassMatrix(azimuth, tower, rotor)

persistent  mpt mtb mpb mpp mtt mbb IDT init_flag

azimuth = azimuth + [0, 2*pi/3, 4*pi/3];


if isempty(init_flag)

    % Some system constants
    Ht = 87.6;
    mNH = 110000 + 240000;
    Mnh = mNH;
    N = 97;
    
    Ip = 2.56193E+09;
    IpAM = 7204372480.00000;
    
    Irotor = 35447886 + 115926;
    Igen = 534.116;
    IDT = Irotor + Igen * N^2;

    %%% ===== Ian A Derivation ===== %%%
    mpp = trapz(tower.Ht, tower.Ht.^2 .* tower.TMassDen) + Mnh*Ht^2 + 3*Ht^2*trapz(rotor.BladeStations, rotor.BMassDen) + (Ip + IpAM);
    mtt = trapz(tower.Ht, tower.mode1.^2 .* tower.TMassDen) + Mnh + 3*trapz(rotor.BladeStations, rotor.BMassDen);
    mbb = trapz(rotor.BladeStations, rotor.BMassDen .* rotor.BladeFmode1.^2);

    mpt = trapz(tower.Ht, tower.mode1 .* tower.Ht .* tower.TMassDen) + Mnh*Ht + 3*Ht*trapz(rotor.BladeStations, rotor.BMassDen);
    mpb = Ht * trapz(rotor.BladeStations, rotor.BladeFmode1 .* rotor.BMassDen);
    mtb = trapz(rotor.BladeStations, rotor.BladeFmode1 .* rotor.BMassDen);    

    % Set init flag
    init_flag = 1;
end

% Form mass matrix
M = [   mpp, mpt,    mpb,    mpb,    mpb,   0;
        mpt, mtt,    mtb,    mtb,    mtb,   0;
        mpb, mtb,    mbb,      0,      0,   0;
        mpb, mtb,      0,    mbb,      0,   0;
        mpb, mtb,      0,      0,    mbb,   0;
          0,   0,      0,      0,      0, IDT];

end