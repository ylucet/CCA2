function ovf_where()
% Where does the intermediate blow up, and by how much? The resultant of two conics is a 4x4
% determinant of quadratics, so its entries are degree-8 products of the inputs: coefficients of
% size c give entries of size ~c^8. Measure that directly.
    fprintf('%-28s %12s %12s %14s\n', 'two conics with coeffs ~', 'input max', 'resultant max', 'ratio to 2^53');
    for c = [3 10 30 100 300 1000 3000 10000]
        rng(1);
        worst = 0;
        for t = 1:20
            A = randi([-c c], 1, 6);  B = randi([-c c], 1, 6);
            if all(A(1:3)==0) || all(B(1:3)==0), continue, end
            try
                [~, info] = conicMeet(A, B);
                if ~isempty(info.quartic), worst = max(worst, max(abs(info.quartic))); end
            catch ME
                if strcmp(ME.identifier, 'ratQ:overflow')
                    worst = inf;  break
                end
            end
        end
        if isinf(worst)
            fprintf('%-28d %12d %12s %14s\n', c, c, 'OVERFLOW', '>1');
        else
            fprintf('%-28d %12d %12.3g %14.2g\n', c, c, worst, worst/2^53);
        end
    end
end
