function rowTangProbe()
% Call the REAL algRowFails on the 2x-1 face row of diagonalKink(-2) (the zero-area faces 2 & 4)
% and report. Also the bottom-root variant (FB = -2 style: hi = false) as a sanity check that
% psTangency only fires on the tangential side.
addpath MATLAB;
FP = zeros(5,5,3);
FP(5,1,1)= 1;  FP(4,1,1)= 4;  FP(3,3,1)=-2;  FP(3,2,1)= 8;  FP(3,1,1)=-4;
FP(2,3,1)= 12; FP(2,2,1)=-16; FP(1,5,1)= 1;  FP(1,3,1)=-8;  FP(1,2,1)= 8;
FP(3,1,2)=-2;  FP(2,1,2)=-4;  FP(1,3,2)=-2;  FP(1,1,2)= 4;  FP(1,1,3)= 1;
aTop = struct('FP', FP, 'FB', -1, 'num', [0 0 0 0 0 0 0 2 0 -1], 'den', 1, 'side', 1);
box = [0 1 0 1];
inBox = @(b) true;
fprintf('algRowFails(top row)  = %d\n', algRowFails(aTop, box, inBox, 6));
aBot = struct('FP', FP, 'FB', 1, 'num', [0 0 0 0 0 0 0 2 0 -1], 'den', 1, 'side', -1);
fprintf('algRowFails(bot row)  = %d\n', algRowFails(aBot, box, inBox, 6));
end
