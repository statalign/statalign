package com.ppfold.algo;

/**
 * Calculation of outside variables for a job.
 * 
 * @author Z.Sukosd
 * @see CYKJob
 * @see Inside
 * @see ExpectationMatrixCalc
 */

public class Outside {

	public static JobResults buildOutside(CYKJob cYKJob) {
		
		//final long starttime = System.currentTimeMillis();
		
		PointRes tmp = new PointRes(0, 0); // stores temporary PointRes objects
		PointRes tmp2 = new PointRes(0, 0); // fetcher PointRes

		PointRes Fvalue = new PointRes(0, 0);
		PointRes Lvalue = new PointRes(0, 0);

		PointRes Lval = new PointRes(0, 0);
		PointRes Fval = new PointRes(0, 0);
		PointRes Sval = new PointRes(0, 0);

		int dim = cYKJob.getN(); // dim is dimension of results matrices;
		double[][] prob = cYKJob.getProb();
		int minj = cYKJob.getMinj();
		int mini = cYKJob.getMini();
		int nts = cYKJob.getSeqlength(); // number of nucleotides available to
		// the job

		int totallength = cYKJob.getFullseqlength(); // length of full input
		// sequence.
		// X is container for data we want to find -- all results for the job
		JobResults X = new JobResults(cYKJob.getN());
		ResMatrix L = X.L;
		ResMatrix S = X.S;
		ResMatrix F = X.F;

		ResMatrix partialS = new ResMatrix(dim);
		ResMatrix partialL = new ResMatrix(dim);

		// job.aboveS/F: outside matrices right above job
		// job.diagaboveS/F: outside matrices diagonally above job
		// job.diagbelow: inside L matrix for combination with diagaboveS/F
		// job.below: inside S matrix for combination with aboveS/F
		// coordinate directions are the same for all above/diagabove matrices.

		// first create matrix multiplication results
		// First half of sums: vertical
		for (int c = 0; c <= cYKJob.getAboveF().size() - 1; c++) {
			// S->LS
			partialL.incrementLWithDirectProduct(cYKJob.getAboveS().get(c),
					cYKJob.getBelow().get(c + 1), prob[0][0], tmp, tmp2);
			// F->LS
			partialL.incrementLWithDirectProduct(cYKJob.getAboveF().get(c),
					cYKJob.getBelow().get(c + 1), prob[2][1], tmp, tmp2);
		}

		// Second half of sums: diagonal
		for (int c = 0; c <= cYKJob.getDiagaboveF().size() - 1; c++) {
			// S->LS
			partialS.incrementSWithDirectProduct(cYKJob.getDiagaboveS().get(c),
					cYKJob.getDiagbelow().get(c + 1), prob[0][0], tmp, tmp2);
			// F->LS
			partialS.incrementSWithDirectProduct(cYKJob.getDiagaboveF().get(c),
					cYKJob.getDiagbelow().get(c + 1), prob[2][1], tmp, tmp2);
			;
		}

		for (int t = dim - 1; t >= 0; t--) {
			for (int s = dim - 1; s >= 0; s--) {
				int i = dim - s - 1;
				int j = minj - dim + 1 + (t + s);
				if (j < 0 || j >= -i + nts) {
					// don't calculate points outside range
					continue;
				}
				if (j == totallength - 1) {
					tmp.add(1);
					S.setProb(s, t, tmp);

					tmp.setToFloat(0);
					tmp.add(prob[0][1]);
					L.setProb(s, t, tmp);
					continue;
				}
				if (mini + i + j <= totallength - 2 && i + mini >= 1) {
					if (minj + t < nts) {
						// L -> dFd
						// F -> dFd
						if (s < dim - 1 && t < dim - 1) {
							Fvalue.copyFrom(F.fetchProb(s + 1, t + 1, tmp));
							Lvalue.copyFrom(L.fetchProb(s + 1, t + 1, tmp));
						} else if (t == dim - 1 && s != dim - 1) {
							// get F from above
							Fvalue = cYKJob.getVertF(s + 1);
							Lvalue = cYKJob.getVertL(s + 1);
						} else if (s == dim - 1 & t != dim - 1) {
							// get F from diagabove
							Fvalue = cYKJob.getDiagF(t + 1);
							Lvalue = cYKJob.getDiagL(t + 1);
						} else if (s == dim - 1 && t == dim - 1) {
							// special F
							Fvalue = cYKJob.specialF;
							Lvalue = cYKJob.specialL;
						} else {
							System.out.println("Ooops, I didn't expect this!");
						}

						tmp.copyFrom(Lvalue);
						tmp.multiply(cYKJob.getouterBPProb(s, t).multiply(
								prob[1][0]));
						F.addToProb(s, t, tmp, tmp2);
						
						tmp.copyFrom(Fvalue);
						tmp.multiply(cYKJob.getinnerBPProb(s, t).multiply(
								prob[2][0]));
						F.addToProb(s, t, tmp, tmp2);
					}
				}

				// BIFURICATION
				// INSIDE OWN JOB
				if (minj == -1) {
					// vertical parts: adding to value of S at point
					for (int tracs = s + 1; tracs <= dim - 1; tracs++) {
						// F->LS
						Fval.copyFrom(F.fetchProb(tracs, t, tmp2));
						Lval.copyFrom(cYKJob.insideresultsL.fetchProb(tracs,
								dim - s - 1, tmp2));
						tmp.copyFrom(Fval);
						tmp.multiply(Lval, prob[2][1]);
						S.addToProb(s, t, tmp, tmp2);

						// S->LS
						Fval.copyFrom(S.fetchProb(tracs, t, tmp));
						tmp.copyFrom(Fval);
						tmp.multiply(Lval, prob[0][0]);
						S.addToProb(s, t, tmp, tmp2);
					}
					// diagonal parts: adding to value of L at point
					for (int tract = t + 1; tract <= dim - 1; tract++) {
						// F->LS
						Fval.copyFrom(F.fetchProb(s, tract, tmp2));
						Sval.copyFrom(cYKJob.insideresultsS.fetchProb(dim - t
								- 1, tract, tmp2));
						tmp.copyFrom(Fval);
						tmp.multiply(Sval, prob[2][1]);
						L.addToProb(s, t, tmp, tmp2);
						// S->LS
						Fval.copyFrom(S.fetchProb(s, tract, tmp2));
						tmp.copyFrom(Fval);
						tmp.multiply(Sval, prob[0][0]);
						L.addToProb(s, t, tmp, tmp2);
					}
					S.addToProb(s, t, partialS.fetchProb(s, t, tmp), tmp2);
					L.addToProb(s, t, partialL.fetchProb(s, t, tmp), tmp2);
				} else {
					// diagonal parts, adding to value of L
					for (int tracs = s + 1; tracs <= dim - 1; tracs++) {
						// F->LS
						Fval.copyFrom(F.fetchProb(tracs, t, tmp2));
						Lval.copyFrom(cYKJob.getDiagbelow().get(0).fetchProb(
								tracs, dim - s - 1, tmp2));
						tmp.copyFrom(Fval);
						tmp.multiply(Lval, prob[2][1]);
						S.addToProb(s, t, tmp, tmp2);
						Fval.copyFrom(S.fetchProb(tracs, t, tmp2));
						tmp.copyFrom(Fval);
						tmp.multiply(Lval, prob[0][0]);
						S.addToProb(s, t, tmp, tmp2);
					}

					// vertical parts: adding to value of S at point
					for (int tract = t + 1; tract <= dim - 1; tract++) {
						// F->LS
						Fval.copyFrom(F.fetchProb(s, tract, tmp2));
						Sval.copyFrom(cYKJob.getBelow().get(0).fetchProb(
								dim - t - 1, tract, tmp2));
						tmp.copyFrom(Fval);
						tmp.multiply(Sval, prob[2][1]);
						L.addToProb(s, t, tmp, tmp2);

						// S->LS
						Fval.copyFrom(S.fetchProb(s, tract, tmp2)); // this is
						// actually
						// outside S
						tmp.copyFrom(Fval);
						tmp.multiply(Sval, prob[0][0]);
						L.addToProb(s, t, tmp, tmp2);
					}
					// big matrix products
					tmp.copyFrom(partialS.fetchProb(s, t, tmp2));
					S.addToProb(s, t, tmp, tmp2);
					tmp.copyFrom(partialL.fetchProb(s, t, tmp2));
					L.addToProb(s, t, tmp, tmp2);
				}

				// S->L
				tmp.copyFrom(S.fetchProb(s, t, tmp2));
				tmp.multiply(prob[0][1]);
				L.addToProb(s, t, tmp, tmp2);
			}

		}
	//	System.out.println(cYKJob.sectorid + " needs " + cYKJob.getDataTransportRequirement()
	//			+  " time-intensiveness: " + cYKJob.getExecutionTimeEstimate()
	//			+ " actual exec time " + (System.currentTimeMillis()-starttime));		
		return X;
	}
}