
tMax = 60;
tDur = 1e5;
% tMax = lt;
% tDur = Dt;

% total - initial = excess moles of carbon during perturbation period

FinTotal = sum(Findt(1:tMax)*Aoc)-Fint(1)*Aoc*tDur; 
FSiTotal = sum(FSidt(1:tMax)*Aoc)-FSit(1)*Aoc*tDur;


for kn = 1 : tMax-1
    myFindt   (kn) = ((Fint   (kn+1)+Fint   (kn))*dtv(kn)/2)*Aoc; % mol C
    myFSidt   (kn) = ((FSit   (kn+1)+FSit   (kn))*dtv(kn)/2)*Aoc; % mol C
    
    myFprdtA  (kn) = ((FprtA  (kn+1)+FprtA(  kn))*dtv(kn)/2)*A(1); % mol C
    myFprdtI  (kn) = (FprtI  (kn+1)+FprtI  (kn))*dtv(kn)/2*A(2); % mol C
    myFprdtP  (kn) = (FprtP  (kn+1)+FprtP  (kn))*dtv(kn)/2*A(3); % mol C
    myFprdtT  (kn) = (FprtT  (kn+1)+FprtT  (kn))*dtv(kn)/2*A(11); % mol C
end

myFinT   =  sum(myFindt);                    % mol
myFSiT   =  sum(myFSidt);                    % mol
myFprT   =  sum(  sum(myFprdtA) ...        % mol
    +sum(myFprdtI) ...
    +sum(myFprdtP) ...
    +sum(myFprdtT));
myFprExc = myFprT - sum(FprtA(1)*A(1)*tDur+FprtI(1)*A(2)*tDur + FprtP(1)*A(3)*tDur + FprtT(1)*A(11)*tDur);

for i = 1 : tMax-1
    for k=1:Ns
        myDisstA(i,k)   = ...
            (dissvtA(i+1,k)+dissvtA(i,k))*dtv(i)/2; % kg
        myDisstI(i,k)   = ...
            (dissvtI(i+1,k)+dissvtI(i,k))*dtv(i)/2; % kg
        myDisstP(i,k)   = ...
            (dissvtP(i+1,k)+dissvtP(i,k))*dtv(i)/2; % kg
    end
end

myDissA    = sum(sum(DisstA))  /VsA; % kg /m3 sed
myDissI    = sum(sum(DisstI))  /VsI; % kg /m3 sed
myDissP    = sum(sum(DisstP))  /VsP; % kg /m3 sed
myDissT    = sum(sum(DisstT))  /VsT; % kg /m3 sed

myDissCA   = myDissA*VsA/Voc/m2kg;     % mol/m3 ocean
myDissCI   = myDissI*VsI/Voc/m2kg;     % mol/m3 ocean
myDissCP   = myDissP*VsP/Voc/m2kg;     % mol/m3 ocean
myDissCT   = myDissT*VsT/Voc/m2kg;     % mol/m3 ocean

myDDissC   = sum(myDissCA+myDissCI+myDissCP+myDissCT);

riverineFluxes = myFinT*tDur + myFSiT*tDur