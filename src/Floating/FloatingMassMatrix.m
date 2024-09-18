function M = FloatingMassMatrix(btable, ttable, azimuth, omega)

persistent mpp_con mss mst msp msb mhh mpp_mul1 mpp_mul2 mpt mpb_con mpb_mul mtb mpb mpp mtt mbb IDT init_flag

azimuth = azimuth + [0, 2*pi/3, 4*pi/3];

Ht = 87.6;
mNH = 110000 + 240000;
Mnh = mNH;
Mp = 3.85218E+06;
N = 97;

Ip = 2.56193E+09;
IpAM = 7204372480.00000;

Irotor = 35447886 + 115926;
Igen = 534.116;
IDT = Irotor + Igen * N^2;

if isempty(init_flag)
    %%% ===== Trinity Paper As-Written ===== %%%
    % % Platform Pitch
    % mpp_con = 12 * Ht^2 * trapz(btable.Rb, btable.BMassDen) + Ht^2*mNH + Ip + IpAM + trapz(ttable.Ht, ttable.Ht .* ttable.TMassDen);
    % mpp_mul1 = 4 * Ht * trapz(btable.Rb, btable.Rb .* btable.BMassDen);
    % mpp_mul2 = trapz(btable.Rb, btable.Rb.^2 .* btable.BMassDen);
    % 
    % % Pitch-tower coupling term
    % mpt = 6 * Ht * trapz(btable.Rb, btable.BMassDen) + trapz(ttable.Ht, ttable.Ht .* ttable.TMassDen .* ttable.mode1);
    % 
    % % Tower-blade coupling term
    % mtb = trapz(btable.Rb, btable.Fmode1 .* btable.BMassDen);
    % 
    % % Tower mass term
    % mtt = 3 * trapz(btable.Rb, btable.BMassDen) + trapz(ttable.Ht, ttable.mode1.^2 .* ttable.TMassDen) + mNH;
    % 
    % % Blade mass term (THIS IS A GUESS FOR NOW)
    % mbb = trapz(btable.Rb, btable.Fmode1.^2 .* btable.BMassDen);
    % 
    % % Pitch-blade coupling terms
    % mpb_con = 2 * Ht * trapz(btable.Rb, btable.Fmode1 .* btable.BMassDen);
    % mpb_mul = trapz(btable.Rb, btable.Rb .* btable.Fmode1 .* btable.BMassDen);

    % %%% ===== Ian A Derivation ===== %%%
    mpp = trapz(ttable.Ht, ttable.Ht.^2 .* ttable.TMassDen) + Mnh*Ht^2 + 3*Ht^2*trapz(btable.Rb, btable.BMassDen) + (Ip + IpAM);
    mtt = trapz(ttable.Ht, ttable.mode1.^2 .* ttable.TMassDen) + Mnh + 3*trapz(btable.Rb, btable.BMassDen);
    mbb = trapz(btable.Rb, btable.Fmode1.^2 .* btable.BMassDen);

    mpt = trapz(ttable.Ht, ttable.mode1 .* ttable.Ht .* ttable.TMassDen) + Mnh*Ht + 3*Ht*trapz(btable.Rb, btable.BMassDen);
    mpb = Ht * trapz(btable.Rb, btable.Fmode1 .* btable.BMassDen);
    mtb = trapz(btable.Rb, btable.Fmode1 .* btable.BMassDen);

    % Set init flag
    init_flag = 1;
end

%%% ===== Trinity Paper As-Written ===== %%%
% % Platform Pitch
% mpp = mpp_con + sum(mpp_mul1 * cos(azimuth) + mpp_mul2 * cos(azimuth).^2);
% 
% % Pitch-blade coupling terms
% mpb = cos(azimuth) * mpb_mul + mpb_con;
% 
% % Form mass matrix
% M = [   mpp, mpt, mpb(1), mpb(2), mpb(3),   0;
%         mpt, mtt,    mtb,    mtb,    mtb,   0;
%      mpb(1), mtb,    mbb,      0,      0,   0;
%      mpb(2), mtb,      0,    mbb,      0,   0;
%      mpb(3), mtb,      0,      0,    mbb,   0;
%           0,   0,      0,      0,      0, IDT];

%%% ===== Ian A - First Derivation ===== %%%
% Form mass matrix
M = [   mpp, mpt,    mpb,    mpb,    mpb,   0;
        mpt, mtt,    mtb,    mtb,    mtb,   0;
        mpb, mtb,    mbb,      0,      0,   0;
        mpb, mtb,      0,    mbb,      0,   0;
        mpb, mtb,      0,      0,    mbb,   0;
          0,   0,      0,      0,      0, IDT];