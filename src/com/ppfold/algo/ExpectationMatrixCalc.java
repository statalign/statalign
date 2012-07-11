package com.ppfold.algo;

/**
 * Static class for expectation matrix calculations.
 * 
 * @author Z.Sukosd
 * @see CYKJob
 * @see Inside
 * @see Outside
 */

public class ExpectationMatrixCalc {

	public static JobResults buildExpectation(CYKJob job) {
		//final long starttime = System.currentTimeMillis();
		
		PointRes tmp = new PointRes(0, 0);
		PointRes tmp2 = new PointRes(0, 0);

		PointRes Evalue = new PointRes(0, 0);
		PointRes Eval1 = new PointRes(0, 0);
		PointRes Eval2 = new PointRes(0, 0);

		int dim = job.getN(); // dim is dimension of results matrices;
		int minj = job.getMinj();
		int nts = job.getSeqlength(); // number of nucleotides available to the
		// job

		// X is container for data we want to find -- all results for the job
		JobResults X = new JobResults(job.getN());
		ResMatrix E = X.S;
		ResMatrix Bp = job.insideresultsS; // own matrix basepairs

		ResMatrix partialE = new ResMatrix(dim);

		// first create matrix multiplication results
		if (job.getBelow().size() != 0) {
			for (int c = 1; c <= job.getBelow().size() - 1; c++) {
				// below: E matrix of that sector
				// diagbelow: E matrix of that sector
				partialE.replaceWithMaxSum(job.getBelow().get(
						job.getBelow().size() - c - 1), job.getDiagbelow().get(
						c - 1), tmp, tmp2);
			}
		}

		for (int t = 0; t < dim; t++) {
			for (int s = 0; s <= dim - 1; s++) {
				int i = dim - s - 1;
				int j = minj - dim + 1 + (t + s);
				// don't calculate points that are outside the range
				if (j < 0 || j >= -i + nts) {
					continue;
				}
				// initialize 0th row
				if (j == 0) {
					tmp.copyFrom(Bp.fetchProb(s, t, tmp2));
					E.setProb(s, t, tmp);
				}
				// in this case it is possible for the nucleotide to basepair.
				if (j >= 2) {
					// basepairing happens: E(i,j) = E(i+1, j-1) + Bp(i,j)
					if (minj + t < nts) {
						if (s >= 1 && t >= 1) {
							// get E from own matrix
							Evalue.copyFrom(E.fetchProb(s - 1, t - 1, tmp));
						} else if (t == 0 && s != 0) {
							// get E from below
							Evalue = job.getVertF(s - 1);
						} else if (s == 0 & t != 0) {
							// get E from diagbelow
							Evalue = job.getDiagF(t - 1);
						} else if (s == 0 && t == 0) {
							// get the 'special' E
							Evalue = job.specialF;
						}

						tmp.copyFrom(Bp.fetchProb(s, t, tmp2));
						tmp.multiply(2);
						tmp.add(Evalue);

						if (E.fetchProb(s, t, tmp2).isLessThan(tmp)) {
							E.setProb(s, t, tmp);
						}
					}
				}

				// Bifurcation
				// first row, need to bifurcate inside own matrix
				if (job.getBelow().size() == 0) {
					for (int k = dim - s; k <= t - 1; k++) {
						Eval1.copyFrom(E.fetchProb(s, k, tmp));
						Eval2.copyFrom(E.fetchProb(dim - k - 1, t, tmp));
						tmp.copyFrom(Eval1);
						tmp.add(Eval2);

						if (E.fetchProb(s, t, tmp2).isLessThan(tmp)) {
							E.setProb(s, t, tmp);
						}
					}
				}

				// extra1
				if (job.getBelow().size() > 0) {
					for (int k = s - 1; k >= 0; k--) {
						Eval1.copyFrom(job.getBelow().get(
								job.getBelow().size() - 1).fetchProb(s,
								dim - k - 1, tmp));
						Eval2.copyFrom(E.fetchProb(k, t, tmp));
						tmp.copyFrom(Eval1);
						tmp.add(Eval2);

						if (E.fetchProb(s, t, tmp2).isLessThan(tmp)) {
							E.setProb(s, t, tmp);
						}
					}

					// big matrix products

					tmp.copyFrom(partialE.fetchProb(s, t, tmp2));

					if (E.fetchProb(s, t, tmp2).isLessThan(tmp)) {
						E.setProb(s, t, tmp);
					}

					// extra2
					for (int k = 0; k <= t - 1; k++) {
						Eval1.copyFrom(job.getDiagbelow().get(
								job.getDiagbelow().size() - 1).fetchProb(
								dim - k - 1, t, tmp));
						Eval2.copyFrom(E.fetchProb(s, k, tmp));
						tmp.copyFrom(Eval1);
						tmp.add(Eval2);

						if (E.fetchProb(s, t, tmp2).isLessThan(tmp)) {
							E.setProb(s, t, tmp);
						}
					}
				}

			}
		}
	//	System.out.println(job.sectorid + " needs " + job.getDataTransportRequirement()
	//			+ " time-intensiveness: " + job.getExecutionTimeEstimate()
	//			+ " actual exec time " + (System.currentTimeMillis()-starttime));
		return X;
	}
}
