function probeRun()
  fprintf('--- can the lattice be expressed with Alg inserted ABOVE Rat? ---\n');
  q = PQuaCon([1 2 3]);
  fprintf('PQuaCon built            : f=%s den=%s\n', mat2str(q.f), mat2str(q.den));
  fprintf('  isa Alg   %d   (a rational function IS algebraic -- must be TRUE)\n', isa(q,'PAlg'));
  fprintf('  isa Rat   %d   isa Qua %d   isa Conic %d\n', isa(q,'PRat'), isa(q,'PQua'), isa(q,'PConic'));
  fprintf('  kind()    %s\n', q.kind());

  b = PBiconj([4 5 6], [1;2]);
  fprintf('PBiconj built            : f=%s FK=%s\n', mat2str(b.f), mat2str(b.FK));
  fprintf('  isa Alg   %d   isa Rat %d  (must be TRUE then FALSE)\n', isa(b,'PAlg'), isa(b,'PRat'));

  fprintf('--- the PINNING: a Rat must not be able to carry FK ---\n');
  try
     q.FK = [1;2];
     fprintf('  FAIL: a PQuaCon accepted FK\n');
  catch ME
     fprintf('  refused by name: %s\n', ME.identifier);
  end

  fprintf('--- does the leaf-only constructor protocol still hold through the deeper chain? ---\n');
  e = PQuaCon.empty();
  fprintf('  PQuaCon.empty() ok, size %s\n', mat2str(size(e)));
  z = PQuaCon();
  fprintf('  no-arg PQuaCon ok, f empty = %d\n', isempty(z.f));
end
