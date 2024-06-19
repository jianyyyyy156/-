
% Arakawa 1966 Eq(46), but note that i, j in the equation are to the East and North respectively,
%  but here jy is to the south therefore j+1��in the original eq 46 is jy-1 in our martrix, j - 1 is jy+1 to the
% South,   % J1=Jacobian(zeta,psi)

J0(jyy,ixx)=.......
    (psi(jyyp,ixx)+psi(jyyp,ixxp)-psi(jyym,ixx)-psi(jyym,ixxp)).*(zeta(jyy,ixxp)+zeta(jyy,ixx))-....
    (psi(jyyp,ixxm)+psi(jyyp,ixx)-psi(jyym,ixxm)-psi(jyym,ixx)).*(zeta(jyy,ixx)+zeta(jyy,ixxm))+.....
    (psi(jyy,ixxp)+psi(jyym,ixxp)-psi(jyy,ixxm)-psi(jyym,ixxm)).*(zeta(jyym,ixx)+zeta(jyy,ixx))-.....
    (psi(jyyp,ixxp)+psi(jyy,ixxp)-psi(jyyp,ixxm)-psi(jyy,ixxm)).*(zeta(jyy,ixx)+zeta(jyyp,ixx))+.....
    (psi(jyy,ixxp)-psi(jyym,ixx)).*(zeta(jyym,ixxp)+zeta(jyy,ixx))-.....
    (psi(jyyp,ixx)-psi(jyy,ixxm)).*(zeta(jyy,ixx)+zeta(jyyp,ixxm))+.....
    (psi(jyym,ixx)-psi(jyy,ixxm)).*(zeta(jyym,ixxm)+zeta(jyy,ixx))-.....
    (psi(jyy,ixxp)-psi(jyyp,ixx)).*(zeta(jyy,ixx)+zeta(jyyp,ixxp));

