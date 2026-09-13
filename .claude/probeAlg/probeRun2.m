function probeRun2()
% Follow-up to probeRun: can the ROOT of the lattice be CONCRETE?
% RatCon is abstract today and QuaCon is the concrete cell at (Qua, Conic). At (Alg, Conic) the
% root and the concrete cell are the SAME cell, so either AlgCon is concrete or a second name is
% needed. This asks whether MATLAB lets an abstract class inherit a concrete one, and whether the
% leaf-only constructor protocol still holds when the ROOT has a real constructor.
  fprintf('--- root CONCRETE, with a real constructor ---\n');
  b = CAlgCon([1 2 3], [7;8]);
  fprintf('  CAlgCon built directly   : f=%s FP=%s  isa Alg %d  isa Rat %d\n', ...
      mat2str(b.f), mat2str(b.FP), isa(b,'PAlg'), isa(b,'PRat'));
  q = CQuaCon([4 5 6]);
  fprintf('  CQuaCon (abstract mid)   : f=%s  isa Alg %d  isa Rat %d  isa Qua %d\n', ...
      mat2str(q.f), isa(q,'PAlg'), isa(q,'PRat'), isa(q,'PQua'));
  fprintf('  the ROOT constructor did NOT clobber the leaf: f = %s (want [4 5 6])\n', mat2str(q.f));
  try
     q.FP = [1;2];
     fprintf('  FAIL: a Rat accepted FP\n');
  catch ME
     fprintf('  pinning still refuses by name: %s\n', ME.identifier);
  end
end
