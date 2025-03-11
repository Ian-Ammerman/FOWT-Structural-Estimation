function M = NREL5MW_MassMatrix(azimuth, tower, rotor, floater)

persistent  mpt mtb mpb mpp mtt mbb init_flag

% azimuth = azimuth + [0, 2*pi/3, 4*pi/3];


if isempty(init_flag)

    % Some system constants
    Ht = max(tower.Ht);
    Mnh = tower.NacHubMass;
    
    Ip = floater.Ip;
    IpAM = floater.IpAM;

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
M = [   mpp, mpt,    mpb,    mpb,    mpb;
        mpt, mtt,    mtb,    mtb,    mtb;
        mpb, mtb,    mbb,      0,      0;
        mpb, mtb,      0,    mbb,      0;
        mpb, mtb,      0,      0,    mbb];

end