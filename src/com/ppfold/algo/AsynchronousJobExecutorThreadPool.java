package com.ppfold.algo;

import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

/**
 * Asynchronous job execution, thread pool class. Jobs are executed in a fixed
 * thread pool.
 * 
 * @author M.Vaerum
 */

public class AsynchronousJobExecutorThreadPool extends AsynchronousJobExecutor {
	private ExecutorService threadPool;
	private int threadCount;

	public AsynchronousJobExecutorThreadPool(int threadCount) {
		this.threadCount = threadCount;
		threadPool = Executors.newFixedThreadPool(threadCount);
	}

	@Override
	public void startExecution(final CYKJob cYKJob, final JobListener listener) {
		threadPool.execute(new Runnable() {
			public void run() {
				JobResults res = null;
				if (cYKJob.isType() == 0) {
					// System.out.println(Thread.currentThread() +
					// " executes job no "+ cYKJob.getSectorid());
					//long before = System.currentTimeMillis();
					// System.out.println(Thread.currentThread() +
					// " starts   job no "+ cYKJob.getSectorid());
					// System.out.println(Thread.currentThread() +
					// " starts   job no "+ cYKJob.getSectorid()+ " at " +
					// System.nanoTime());
					res = Inside.buildInside(cYKJob);
					// System.out.println(Thread.currentThread() +
					// " finishes job no "+ cYKJob.getSectorid()+ ": "+
					// (System.currentTimeMillis() - before) +" ms");
					// System.out.println(Thread.currentThread() +
					// " finishes job no "+ cYKJob.getSectorid()+ " at " +
					// System.nanoTime());
					listener.jobFinished(res);
				} else if (cYKJob.isType() == 1) {
					res = Outside.buildOutside(cYKJob);
					listener.jobFinished(res);
				} else {
					res = ExpectationMatrixCalc.buildExpectation(cYKJob);
					listener.jobFinished(res);
				}
			}
		});
	}

	@Override
	public void startExecution(final PhyloJob job, final JobListener listener) {
		threadPool.execute(new Runnable() {
			public void run() {
				if (job.isType() == false) {
					double[][] res = PhyloCalc.calcSingleColumn(job);
					listener.jobFinished(res);
				} else {
					double[][] res = PhyloCalc.calcDoubleColumn(job);
					listener.jobFinished(res);
				}
			}
		});
	}

	@Override
	public String getDescription() {
		return "Local thread execution. Bounded by " + threadCount + " threads";
	}

	@Override
	public String getId() {
		return this.getClass().getName() + threadCount;
	}

	public void shutDown(){
		this.threadPool.shutdownNow(); 
	}

	public boolean isTerminated() {
		return this.threadPool.isTerminated();
	}
	
	// @Override
	// public AsynchronousJobExecutor createJobExecutor(AlgoParameters param) {
	// //we don't have state
	// return this;
	// }

}
