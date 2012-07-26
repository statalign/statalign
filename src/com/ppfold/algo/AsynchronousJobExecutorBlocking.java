package com.ppfold.algo;

import java.util.List;

/**
 * Asynchronous job execution, blocking class. Jobs are executed sequentially.
 * 
 * @author M.Vaerum
 */

public class AsynchronousJobExecutorBlocking extends AsynchronousJobExecutor {
	@Override
	public String getDescription() {
		return "Serial (blocking) execution";
	}

	@Override
	public void startExecution(CYKJob cYKJob, JobListener listener) {
		if (cYKJob.isType() == 0) {
			JobResults res = Inside.buildInside(cYKJob);
			listener.jobFinished(res);
		} else if (cYKJob.isType() == 1) {
			JobResults res = Outside.buildOutside(cYKJob);
			listener.jobFinished(res);
		} else {
			JobResults res = ExpectationMatrixCalc.buildExpectation(cYKJob);
			listener.jobFinished(res);
		}
	}

	@Override
	public void startExecution(PhyloJob job, JobListener listener) {
		if (job.isType() == false) {
			double[][] res = PhyloCalc.calcSingleColumn(job);
			listener.jobFinished(res);
		} else {
			double[][] res = PhyloCalc.calcDoubleColumn(job);
			listener.jobFinished(res);
		}
	}
	
	@Override
	public void startExecution(PhyloJobFuzzy job, JobListener listener) {
		if (job.isType() == false) {
			double[][] res = PhyloCalcFuzzy.calcSingleColumn(job);
			listener.jobFinished(res);
		} else {
			double[][] res = PhyloCalcFuzzy.calcDoubleColumn(job);
			listener.jobFinished(res);
		}
	}

	@Override
	public String getId() {
		return this.getClass().getName();
	}

	public void shutDown(){
		//Nothing happens here as we only have one thread. 
	}
	
	public boolean isTerminated() {
		return false;
	}
	
	// @Override
	// public AsynchronousJobExecutor createJobExecutor(AlgoParameters param) {
	// //We have no state
	// return this;
	// }

}
