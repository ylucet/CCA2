# Port scope -- measured, not guessed

Line coverage of the production `.m` files under the END-TO-END exact-engine suites, MATLAB R2024a,
2026-11-10 (`python/tools/coverage_scope.m`). "conj" = conjQTest + conjQExamplesSweepTest + maxQTest (57/0);
"biconj" = biconjQTest + biconjAtomsTest + erratumExamplesTest + erratumBiconjTest (59/0).
A file with no executed line under either is OUT of the port. Test files and oracle-only files
(legacy engine called by the tests as a cross-check) are listed but marked.

| file | conj hit/total | biconj hit/total |
|---|---|---|
| ratQ.m | 306/2699 | 1526/2699 |
| QuaPol.m | 898/1340 | 856/1340 |
| conjQTest.m | 1167/1169 | 0/1169 |
| conjQ.m | 543/598 | 308/598 |
| biconjQTest.m | 0/788 | 788/788 |
| algRowFails.m | 0/815 | 624/815 |
| envelopeSubdivisionAlg.m | 0/706 | 515/706 |
| biconjAtoms.m | 0/481 | 466/481 |
| biconjQ.m | 0/510 | 437/510 |
| conicMeet.m | 167/178 | 167/178 |
| assembleQuaConCells.m | 179/184 | 133/184 |
| QuaPar.m | 288/1028 | 18/1028 |
| erratumExamplesTest.m | 0/308 | 262/308 |
| RatPol.m | 2/796 | 252/796 |
| QuaCon.m | 124/248 | 120/248 |
| cellInterior.m | 42/216 | 199/216 |
| bridgeFaceAlg.m | 0/451 | 226/451 |
| convEnvCPLQ.m | 0/234 | 155/234 |
| envelopeSubdivision.m | 0/197 | 139/197 |
| bridgeFaceRat.m | 0/194 | 126/194 |
| algBranch.m | 0/157 | 121/157 |
| erratumBiconjTest.m | 0/126 | 120/126 |
| cornerFaceAlg.m | 0/130 | 115/130 |
| maxQTest.m | 106/106 | 0/106 |
| cellHasInterior.m | 52/61 | 53/61 |
| cornerWedgeRows.m | 0/125 | 103/125 |
| AlgCrv.m | 0/126 | 102/126 |
| envelopeCellRegions.m | 0/95 | 92/95 |
| bridgeWitness.m | 0/176 | 90/176 |
| envelopeEval.m | 0/108 | 85/108 |
| compactCells.m | 41/47 | 41/47 |
| envelopeAssemble.m | 0/115 | 78/115 |
| isConvexQ.m | 0/96 | 72/96 |
| biconjCPLQ.m | 0/178 | 71/178 |
| bigZ.m | 0/202 | 68/202 |
| conjPieceCPLQ.m | 66/337 | 0/337 |
| bridgeMember.m | 0/72 | 65/72 |
| checkQuaConConsistent.m | 31/44 | 31/44 |
| mergeAdjacentCells.m | 31/31 | 31/31 |
| AlgCon.m | 25/109 | 35/109 |
| algFaceBelow.m | 0/66 | 59/66 |
| maxQ.m | 29/29 | 27/29 |
| assembleAlgCrv.m | 0/76 | 56/76 |
| EnvelopeMin.m | 0/56 | 56/56 |
| domainIsConvexD.m | 26/30 | 26/30 |
| envelopePrune.m | 0/52 | 49/52 |
| conjCPLQ.m | 45/346 | 0/346 |
| cellSamplePoint.m | 0/47 | 43/47 |
| conjQExamplesSweepTest.m | 42/44 | 0/44 |
| envelopeCellPairs.m | 0/50 | 40/50 |
| biconjAtomsTest.m | 0/40 | 40/40 |
| isolatedRootSign.m | 0/43 | 37/43 |
| envelopeRationalCorners.m | 0/33 | 33/33 |
| polySignOn.m | 0/37 | 31/37 |
| conjEdgeLowerBound.m | 25/32 | 0/32 |
| bridgePhi.m | 0/33 | 22/33 |
| mergeSameQuadFaces.m | 6/47 | 6/47 |
| bridgeArcCut.m | 0/69 | 10/69 |
