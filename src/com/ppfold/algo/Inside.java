package com.ppfold.algo;

/**
 * Calculation of inside variables for a job.
 * 
 * @author Z.Sukosd
 * @see CYKJob
 * @see Outside
 * @see ExpectationMatrixCalc
 */

public class Inside {

	public static JobResults buildInside(CYKJob cYKJob) {
		
		//final long starttime = System.currentTimeMillis();

		PointRes tmp = new PointRes(0, 0); // stores temporary PointRes objects
		PointRes tmp2 = new PointRes(0, 0); // fetcher PointRes

		PointRes Fvalue = new PointRes(0, 0);
		PointRes Lval = new PointRes(0, 0);
		PointRes Sval = new PointRes(0, 0);

		int dim = cYKJob.getN(); // dim is dimension of results matrices;
		double[][] prob = cYKJob.getProb();
		int minj = cYKJob.getMinj();
		int mini = cYKJob.getMini();
		int nts = cYKJob.getSeqlength(); // number of nucleotides available to
		// the job

		// X is container for data we want to find -- all results for the job
		JobResults X = new JobResults(cYKJob.getN());
		ResMatrix L = X.L;
		ResMatrix S = X.S;
		ResMatrix F = X.F;

		ResMatrix partialS = new ResMatrix(dim);
		ResMatrix partialF = new ResMatrix(dim);

		// first create matrix multiplication results
		if (cYKJob.getBelow().size() != 0) {
			for (int c = 1; c <= cYKJob.getBelow().size() - 1; c++) {
				partialS.incrementWithLS(cYKJob.getBelow().get(
						cYKJob.getBelow().size() - c - 1), cYKJob
						.getDiagbelow().get(c - 1), prob[0][0], tmp, tmp2);
				partialF.incrementWithLS(cYKJob.getBelow().get(
						cYKJob.getBelow().size() - c - 1), cYKJob
						.getDiagbelow().get(c - 1), prob[2][1], tmp, tmp2);
			}
		}

		for (int t = 0; t < dim; t++) {
			for (int s = 0; s <= dim - 1; s++) {
				int i = dim - s - 1;
				int j = minj - dim + 1 + (t + s);
				// don't calculate points outside range
				if (j < 0 || j >= -i + nts) {
					continue;
				}
				// initialize 0th row
				if (j == 0) {
					// L->s
					tmp.copyFrom(L.fetchProb(s, t, tmp2));
					tmp.add(cYKJob.getBPProb(s, t).multiply(
							prob[1][1]));
					
					//tmp.add(cYKJob.basepairs.getProb(s, t).multiply(
					//				prob[1][1]));
					L.setProb(s, t, tmp);
				}
				// can do basepairing here
				if (j >= 2) {
					if (minj + t < nts) {
						// L -> dFd
						// F -> dFd
						if (s >= 1 && t >= 1) {
							// get F from own matrix
							Fvalue.copyFrom(F.fetchProb(s - 1, t - 1, tmp2));
						} else if (t == 0 && s != 0) {
							// get F from below
							Fvalue = cYKJob.getVertF(s - 1);
						} else if (s == 0 & t != 0) {
							// get F from diagbelow
							Fvalue = cYKJob.getDiagF(t - 1);
						} else if (s == 0 && t == 0) {
							// get the 'special' F
							Fvalue = cYKJob.specialF;
						}

						//L->dFd (outer pairs)
						tmp.copyFrom(Fvalue);
						tmp.multiply(prob[1][0]).multiply(
								cYKJob.getouterBPProb(s, t));
						L.addToProb(s, t, tmp, tmp2);
						
						//F->dFd (inner pairs)
						tmp.copyFrom(Fvalue);
						tmp.multiply(prob[2][0]).multiply(
								cYKJob.getinnerBPProb(s, t));
						F.addToProb(s, t, tmp, tmp2);
					}
				}
				// Bifurcation
				// first row, need to bifurcate inside own matrix
				if (cYKJob.getBelow().size() == 0) {
					for (int k = dim - s; k <= t - 1; k++) {
						Lval.copyFrom(L.fetchProb(s, k, tmp2));
						Sval.copyFrom(S.fetchProb(dim - k - 1, t, tmp2));
						tmp.copyFrom(Lval);
						tmp.multiply(Sval, prob[0][0]);
						S.addToProb(s, t, tmp, tmp2);
						tmp.copyFrom(Lval);
						tmp.multiply(Sval, prob[2][1]);
						F.addToProb(s, t, tmp, tmp2);
					}
				}
				if (cYKJob.getBelow().size() > 0) {
					// extra1
					for (int k = s - 1; k >= 0; k--) {
						Lval.copyFrom(cYKJob.getBelow().get(
								cYKJob.getBelow().size() - 1).fetchProb(s,
								dim - k - 1, tmp2));
						Sval.copyFrom(S.fetchProb(k, t, tmp2));
						tmp.copyFrom(Lval);
						tmp.multiply(Sval, prob[0][0]);
						S.addToProb(s, t, tmp, tmp2);
						tmp.copyFrom(Lval);
						tmp.multiply(Sval, prob[2][1]);
						F.addToProb(s, t, tmp, tmp2);
					}
					// big matrix products
					S.addToProb(s, t, partialS.fetchProb(s, t, tmp), tmp2);
					F.addToProb(s, t, partialF.fetchProb(s, t, tmp), tmp2);
					// extra2
					for (int k = 0; k <= t - 1; k++) {
						Sval.copyFrom(cYKJob.getDiagbelow().get(
								cYKJob.getDiagbelow().size() - 1).fetchProb(
								dim - k - 1, t, tmp2));
						Lval.copyFrom(L.fetchProb(s, k, tmp2));
						tmp.copyFrom(Lval);
						tmp.multiply(Sval, prob[0][0]);
						S.addToProb(s, t, tmp, tmp2);
						tmp.copyFrom(Lval);
						tmp.multiply(Sval, prob[2][1]);
						F.addToProb(s, t, tmp, tmp2);
					}
				}
				// S->L
				tmp.copyFrom(L.fetchProb(s, t, tmp2));
				tmp.multiply(prob[0][1]);
				S.addToProb(s, t, tmp, tmp2);
			}
		}
		
//		System.out.println(cYKJob.sectorid + " needs memory: " + cYKJob.getDataTransportRequirement() +
//				" time-intensiveness: " + cYKJob.getExecutionTimeEstimate()
//				+ " actual exec time " + (System.currentTimeMillis()-starttime));
		return X;
	}
}