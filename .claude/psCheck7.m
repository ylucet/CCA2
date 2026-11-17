function psCheck7()
% Verify MATLAB algRowFails on the 7 trace rows that the Python port now certifies (psTangency
% True, upper/side>0, P+ a negative square). psTangency is a pure function of the row and fires
% before the box recursion, so a trivial box/inBox isolates the certificate.
    addpath MATLAB;

    % ---- row 1 (idx 1 and idx 6 share num/den; idx 1's FP) ----
    FP = zeros(5,5,3);
    FP(1,1,2)=4;  FP(1,1,3)=1;  FP(1,2,1)=8;  FP(1,3,1)=-8;  FP(1,3,2)=-2;  FP(1,5,1)=1;
    FP(2,1,2)=-4; FP(2,2,1)=-16; FP(2,3,1)=12; FP(3,1,1)=-4;  FP(3,1,2)=-2;
    FP(3,2,1)=8;  FP(3,3,1)=-2;  FP(4,1,1)=4;  FP(5,1,1)=1;
    a = struct('FP',FP,'FB',-1,'num',[0 0 0 0 0 0 0 2 0 -1],'den',1,'side',1);
    fprintf('row idx1  algRowFails = %d\n', algRowFails(a,[0 1 0 1],@(~)true,6));

    % ---- row 6 ----
    FP = zeros(5,5,3);
    FP(1,1,3)=1;  FP(1,2,2)=4;  FP(1,3,1)=4;  FP(1,3,2)=-2;  FP(1,4,1)=-4;  FP(1,5,1)=1;
    FP(3,1,2)=-2; FP(3,2,1)=-4; FP(3,3,1)=-2; FP(5,1,1)=1;
    a = struct('FP',FP,'FB',-1,'num',[0 0 0 0 0 0 0 2 0 -1],'den',1,'side',1);
    fprintf('row idx6  algRowFails = %d\n', algRowFails(a,[0 1 0 1],@(~)true,6));

    % ---- row 22 / 36 (same FP, num=[..2 1 -1], den=1) ----
    FP = zeros(5,4,3);
    FP(1,1,3)=1;  FP(1,2,2)=-2; FP(1,3,1)=1;  FP(1,3,2)=-1;  FP(1,4,1)=1;
    FP(2,2,2)=2;  FP(2,3,1)=-2; FP(3,1,2)=-2; FP(3,2,1)=2;   FP(3,3,1)=1;
    FP(4,2,1)=-2; FP(5,1,1)=1;
    a = struct('FP',FP,'FB',-1,'num',[0 0 0 0 0 0 0 2 1 -1],'den',1,'side',1);
    fprintf('row idx22 algRowFails = %d\n', algRowFails(a,[0 1 0 1],@(~)true,6));

    % ---- row 23 / 38 / 41 (same FP, num=[..4 4 -1], den=4) ----
    FP = zeros(5,4,3);
    FP(1,1,3)=1;  FP(1,2,2)=-2; FP(1,3,1)=1;  FP(1,3,2)=-1;  FP(1,4,1)=1;
    FP(2,2,2)=2;  FP(2,3,1)=-2; FP(3,1,2)=-2; FP(3,2,1)=2;   FP(3,3,1)=1;
    FP(4,2,1)=-2; FP(5,1,1)=1;
    a = struct('FP',FP,'FB',-1,'num',[0 0 0 0 0 0 0 4 4 -1],'den',4,'side',1);
    fprintf('row idx23 algRowFails = %d\n', algRowFails(a,[0 1 0 1],@(~)true,6));
end
