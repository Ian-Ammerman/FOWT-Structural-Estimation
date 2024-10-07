function M = FloatingMassMatrix(azimuth, omega)

persistent  mpt mtb mpb mpp mtt mbb IDT init_flag

azimuth = azimuth + [0, 2*pi/3, 4*pi/3];


if isempty(init_flag)

    % Some system constants
    [btable,ttable] = load5MWOC4SemiDistributed;
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
    mpp = trapz(ttable.Ht, ttable.Ht.^2 .* ttable.TMassDen) + Mnh*Ht^2 + 3*Ht^2*trapz(btable.Rb, btable.BMassDen) + (Ip + IpAM);
    mtt = trapz(ttable.Ht, ttable.mode1.^2 .* ttable.TMassDen) + Mnh + 3*trapz(btable.Rb, btable.BMassDen);
    mbb = trapz(btable.Rb, btable.BMassDen .* btable.Fmode1.^2);

    mpt = trapz(ttable.Ht, ttable.mode1 .* ttable.Ht .* ttable.TMassDen) + Mnh*Ht + 3*Ht*trapz(btable.Rb, btable.BMassDen);
    mpb = Ht * trapz(btable.Rb, btable.Fmode1 .* btable.BMassDen);
    mtb = trapz(btable.Rb, btable.Fmode1 .* btable.BMassDen);
    

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