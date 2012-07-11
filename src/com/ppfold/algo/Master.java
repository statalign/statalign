package com.ppfold.algo;

import java.util.concurrent.BlockingQueue;
import java.util.concurrent.LinkedBlockingQueue;

/**
 * Created by FoldingProject, steers asynchronous execution with queues
 * 
 * @author Z.Sukosd
 * @author M.Vaerum
 * @see FoldingProject
 */

public class Master {
	Sector top;
	Sector bottom;
	double[][] prob;

	BlockingQueue<CYKJob> jobChannelInside = new LinkedBlockingQueue<CYKJob>();
	BlockingQueue<CYKJob> jobChannelOutside = new LinkedBlockingQueue<CYKJob>();
	BlockingQueue<CYKJob> jobChannelExpectation = new LinkedBlockingQueue<CYKJob>();

	public Master() {
	}

	public Master(Sector topsec, double[][] prob) {
		top = topsec;
		bottom = findBottom();
		this.prob = prob;
	}

	public boolean outsideFinished() {
		// returns true if outside is finished, otherwise false
		Sector current = bottom;
		while (current.pos[1] == -1 && current.next != null) {
			if (current.getOutsideResults() == null) {
				return false;
			}
			current = current.next;
		}
		return true;
	}

	public void CreateInsideJobChannel() {

		Sector sec = this.bottom;
		// loop to create inside jobchannel
		while (sec != null && sec.below == null) {
			// first job to have below != null is the first one outside the
			// zeroth row
			// ie. enough to check sec.below, loop exists at the right time
			try {
				jobChannelInside.put(createInsideJobFromSector(sec));
				sec.setIsBeingInsideProcessed(true);
			} catch (InterruptedException e) {
				//System.out.println("Process interrupted!");
			}
			sec = sec.next;
		}
	}

	public CYKJob takeNextInsideJob() throws InterruptedException {
		return jobChannelInside.take();
	}

	public void CreateOutsideJobChannel() {
		// create outside jobchannel - consists of top job to start with
		try {
			jobChannelOutside.put(createOutsideJobFromSector(top));
		} catch (InterruptedException e) {
			//System.out.println("Process interrupted!");
		}
	}

	public CYKJob takeNextOutsideJob() throws InterruptedException {
		return jobChannelOutside.take();
	}

	public void CreateExpectationJobChannel() {
		Sector sec = this.bottom;
		// loop to create expectation jobchannel
		while (sec != null && sec.below == null) {
			// first job to have below != null is the first one outside the
			// zeroth row
			// ie. enough to check sec.below, loop exists at the right time
			try {
				jobChannelExpectation.put(createExpectationJobFromSector(sec));
				sec.setIsBeingExpectationProcessed(true);
			} catch (InterruptedException e) {
				//System.out.println("Process interrupted!");
			}
			sec = sec.next;
		}
	}

	public CYKJob takeNextExpectationJob() throws InterruptedException {
		return jobChannelExpectation.take();
	}

	public synchronized boolean unProcessedInsideSectors() {
		return !jobChannelInside.isEmpty()
				|| (top.getInsideResults() == null && !top
						.isBeingInsideProcessed());
	}

	public synchronized boolean unProcessedOutsideSectors() {
		// returns true if:
		// there are jobs in outside queue
		// there is at least one job with getOutsideResults()==null
		// there is at least one job with isBeingOutsideProcessed() = true
		boolean value = false;
		Sector current = bottom;
		while (current.pos[1] == -1 && current.next != null) {
			if (current.getOutsideResults() == null
					&& !current.isBeingOutsideProcessed()) {
				value = true;
			}
			current = current.next;
		}

		return !jobChannelOutside.isEmpty() || value;
	}

	public synchronized boolean unProcessedExpectationSectors() {
		return !jobChannelExpectation.isEmpty()
				|| (top.getExpectationResults() == null && !top
						.isBeingExpectationProcessed());
	}

	public synchronized void setInsideResult(int sectorNumber,
			JobResults results) {
		// Traverse sector structure and find the sector
		Sector sector = findSector(sectorNumber);
		// Update sector with result
		sector.setInsideResults(results);
		sector.setIsBeingInsideProcessed(false);
		// Check to see if a new job is ready and insert into queue
		if (sector.above != null && sector.above.readyForInsideProcessing()) {
			processInside(sector.above);

		}
		if (sector.diagabove != null
				&& sector.diagabove.readyForInsideProcessing()) {
			processInside(sector.diagabove);
		}

	}

	public synchronized void setOutsideResult(int sectorNumber,
			JobResults results) {
		// Traverse sector structure and find the sector
		Sector sector = findSector(sectorNumber);
		// Update sector with result
		sector.setOutsideResults(results);
		sector.setIsBeingOutsideProcessed(false);
		// Check to see if a new job is ready and insert into queue
		if (sector.below != null && sector.below.readyForOutsideProcessing()) {
			processOutside(sector.below);

		}
		if (sector.diagbelow != null
				&& sector.diagbelow.readyForOutsideProcessing()) {
			processOutside(sector.diagbelow);
		}
	}

	public synchronized void setExpectationResult(int sectorNumber,
			JobResults results) {
		// Traverse sector structure and find the sector
		Sector sector = findSector(sectorNumber);
		// Update sector with result
		sector.setExpectationResults(results);
		sector.setIsBeingExpectationProcessed(false);
		// Check to see if a new job is ready and insert into queue
		if (sector.above != null
				&& sector.above.readyForExpectationProcessing()) {
			processExpectation(sector.above);

		}
		if (sector.diagabove != null
				&& sector.diagabove.readyForExpectationProcessing()) {
			processExpectation(sector.diagabove);
		}
	}

	private void processInside(Sector sector) {
		sector.setIsBeingInsideProcessed(true);
		CYKJob nextJob = createInsideJobFromSector(sector);
		// System.out.println("Job created from sector " + sector.sectorid);
		try {
			jobChannelInside.put(nextJob);
		} catch (InterruptedException e) {
			//System.out.println("Process interrupted!");
		}
	}

	private void processOutside(Sector sector) {
		sector.setIsBeingOutsideProcessed(true);
		CYKJob nextJob = createOutsideJobFromSector(sector);
		// System.out.println("Job created from sector " + sector.sectorid);
		try {
			jobChannelOutside.put(nextJob);
		} catch (InterruptedException e) {
			//System.out.println("Process interrupted!");
		}
	}

	private void processExpectation(Sector sector) {
		sector.setIsBeingExpectationProcessed(true);
		CYKJob nextJob = createExpectationJobFromSector(sector);
		// System.out.println("Job created from sector " + sector.sectorid);
		try {
			jobChannelExpectation.put(nextJob);
		} catch (InterruptedException e) {
			//System.out.println("Process interrupted!");
		}
	}

	private CYKJob createInsideJobFromSector(Sector sector) {
		PointRes tmp = new PointRes(0, 0);

		CYKJob cYKJob = new CYKJob(sector.sectorid, sector.dim, 0);
		Sector current = sector;
		while (current.below != null) {
			current = current.below;
			cYKJob.addBelowL(current, false);

		}
		current = sector;
		while (current.diagbelow != null) {
			current = current.diagbelow;
			cYKJob.addDiagBelowS(current, false);
		}
		current = sector.below;
		if (current != null) {
			for (int i = 0; i < sector.dim; i++) {
				cYKJob.setInsideVertF(current, i, tmp);
			}
		}
		current = sector.diagbelow;
		if (current != null) {
			for (int i = 0; i < sector.dim; i++) {
				cYKJob.setInsideDiagF(current, i, tmp);
			}
			Sector other = current.below;
			if (other != null) {
				// cYKJob.specialF =
				// sector.diagbelow.below.getMatrixInsideFTpVal(sector.dim-1,
				// sector.dim-1).clone();
				cYKJob.specialF.copyFrom(sector.diagbelow.below
						.getMatrixInsideFTpVal(sector.dim - 1, sector.dim - 1,
								tmp));

			}
		}

		cYKJob.setData(sector.dim, sector.getMini(), sector.getMinj(),
				sector.seqlength, prob);
		cYKJob.basepairs = sector.getBasePairs();
		if(sector.diffbp){
			cYKJob.basepairs2 = sector.getBasePairs2();
			cYKJob.setDiffBasepairs(true);
		}
		
		return cYKJob;
	}

	private CYKJob createOutsideJobFromSector(Sector sector) {
		PointRes tmp = new PointRes(0, 0);

		CYKJob cYKJob = new CYKJob(sector.sectorid, sector.dim, 1);
		Sector current = sector;

		cYKJob.setFullseqlength(sector.fullseqlength);
		// System.out.println("Pushing INSIDE of job nr. " + sector.sectorid +
		// " into insideresultsL, insideresultsS");
		cYKJob.insideresultsL = sector.getInsideMatrixL();
		cYKJob.insideresultsS = sector.getInsideMatrixS();

		// adding outside jobresults
		while (current.above != null) {
			current = current.above;
			cYKJob.addAboveS(current, true);
			cYKJob.addAboveF(current, true);
		}

		current = sector;
		while (current.diagabove != null) {
			current = current.diagabove;
			cYKJob.addDiagAboveS(current, true);
			cYKJob.addDiagAboveF(current, true);
		}

		// adding inside jobresults
		// 1. adding "far vertical" = "below"

		current = sector;
		while (current.diagbelow != null) {
			// seek to bottom job in sector's diag column
			current = current.diagbelow;

		}
		// add this one: S, inside
		cYKJob.addBelowS(current, false);

		// add all above's of this one
		while (current.above != null) {
			current = current.above;
			cYKJob.addBelowS(current, false);
		}

		// 2. adding "far diagonal" = "diagbelow"
		current = sector;
		while (current.below != null) {
			// seek to bottom job in sector's column
			current = current.below;

		}
		// add this one: L, inside
		cYKJob.addDiagBelowL(current, false);

		// add all diagabove's of this one
		while (current.diagabove != null) {
			current = current.diagabove;
			cYKJob.addDiagBelowL(current, false);
		}

		// adding F values
		current = sector.above;
		if (current != null) {
			for (int i = 0; i < sector.dim; i++) {
				cYKJob.setOutsideVertF(current, i, tmp);
				cYKJob.setOutsideVertL(current, i, tmp);
			}
		}

		current = sector.diagabove;
		if (current != null) {
			for (int i = 0; i < sector.dim; i++) {
				cYKJob.setOutsideDiagF(current, i, tmp);
				cYKJob.setOutsideDiagL(current, i, tmp);
			}
			Sector other = current.above;
			if (other != null) {
				cYKJob.specialF.copyFrom(sector.diagabove.above
						.getMatrixOutsideFTpVal(0, 0, tmp)); // bottom corner
				// this time
				cYKJob.specialL.copyFrom(sector.diagabove.above
						.getMatrixOutsideLTpVal(0, 0, tmp));
			}
		}
		current = sector.above;
		if (current != null) {
		}
		cYKJob.setData(sector.dim, sector.getMini(), sector.getMinj(),
				sector.seqlength, prob);
		cYKJob.basepairs = sector.getBasePairs();
		if(sector.diffbp){
			cYKJob.basepairs2 = sector.getBasePairs2();
			cYKJob.setDiffBasepairs(true);
		}
		return cYKJob;
	}

	private CYKJob createExpectationJobFromSector(Sector sector) {
		PointRes tmp = new PointRes(0, 0);

		CYKJob cYKJob = new CYKJob(sector.sectorid, sector.dim, 2);
		Sector current = sector;

		cYKJob.insideresultsS = sector.getBasePairs(); // own matrix basepairs

		while (current.below != null) {
			current = current.below;
			cYKJob.addBelowE(current);
		}
		current = sector;
		while (current.diagbelow != null) {
			current = current.diagbelow;
			cYKJob.addDiagBelowE(current);
		}
		current = sector.below;
		if (current != null) {
			for (int i = 0; i < sector.dim; i++) {
				cYKJob.setExpectationVertE(current, i, tmp);
			}
		}
		current = sector.diagbelow;
		if (current != null) {
			for (int i = 0; i < sector.dim; i++) {
				cYKJob.setExpectationDiagE(current, i, tmp);
			}
			Sector other = current.below;
			if (other != null) {
				cYKJob.specialF.copyFrom(sector.diagbelow.below
						.getMatrixExpectationTpVal(sector.dim - 1,
								sector.dim - 1, tmp));
			}
		}

		cYKJob.setData(sector.dim, sector.getMini(), sector.getMinj(),
				sector.seqlength, prob);

		return cYKJob;
	}

	public Sector findSector(int sector) {
		Sector current = this.bottom;
		while (current.sectorid != sector) {
			current = current.next;
		}
		return current;
	}

	private Sector findBottom() {
		Sector current = this.top;
		while (current.below != null) {
			current = current.below;
		}
		return current;
	}

	/**
	 * @return the top
	 */
	public Sector getTop() {
		return top;
	}

}