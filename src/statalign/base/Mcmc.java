package statalign.base;

import java.awt.Cursor;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;

import com.ppfold.algo.AlignmentData;
import com.ppfold.algo.FuzzyAlignment;

import mpi.MPI;
import statalign.MPIUtils;
import statalign.base.thread.Stoppable;
import statalign.base.thread.StoppedException;
import statalign.distance.AlignmentSample;
import statalign.distance.Distance;
import statalign.postprocess.PostprocessManager;
import statalign.postprocess.plugins.contree.CNetwork;
import statalign.postprocess.utils.RNAFoldingTools;
import statalign.ui.ErrorMessage;
import statalign.ui.MainFrame;
import statalign.utils.SimpleStats;

/**
 * 
 * This class handles an MCMC run.
 * 
 * The class extends <tt>Stoppable</tt>, it may be terminated/suspended in
 * graphical mode.
 * 
 * @author miklos, novak
 * 
 */
public class Mcmc extends Stoppable {

	// Constants

	// int samplingMethod = 1; //0: random sampling, 1: total sampling
	double[] weights; // for selecting internal tree node
	final static double LEAFCOUNT_POWER = 1.0;
	final static double SELTRLEVPROB[] = { 0.9, 0.6, 0.4, 0.2, 0 };
	final static int FIVECHOOSE[] = { 35, 5, 15, 35, 10 }; // edge, topology,
	// indel parameter,
	// alignment,
	// substitutionparameter
	final static int FOURCHOOSE[] = { 35, 5, 25, 35 }; // edge, topology, indel
	// parameter, alignment
	private static final int SAMPLE_RATE_WHEN_DETERMINING_THE_SPACE = 100;
	private static final int BURNIN_TO_CALCULATE_THE_SPACE = 25000;


	// Parallelization

	/** Is this a parallel chain? By-default false. */
	private boolean isParallel = false;

	/** The number of processes. */
	private int noOfProcesses;

	/** The rank of the process. */
	private int rank;

	/** When this variable reaches 0 we do a swap. */
	private int swapCounter;

	/** The random number generator used for swapping. */
	private Random swapGenerator;

	// Non parallelization

	public CNetwork network; 

	/** Current tree in the MCMC chain. */
	public Tree tree;

	/**
	 * MCMC parameters including the number of burn-in steps, the total number
	 * of steps in the MCMC and the sampling rate.
	 */
	public MCMCPars mcmcpars;

	public McmcStep mcmcStep = new McmcStep();

	/** PostprocessManager that handles the postprocessing modules. */
	public PostprocessManager postprocMan;

	/** True while the MCMC is in the burn-in phase. */
	public boolean burnin;

	public Mcmc(Tree tree, MCMCPars pars, PostprocessManager ppm) {
		postprocMan = ppm;
		ppm.mcmc = this;
		this.tree = tree;
		weights = new double[tree.vertex.length];
		mcmcpars = pars;
		this.tree.heat = 1.0d;
	}

	public Mcmc(Tree tree, MCMCPars pars, PostprocessManager ppm,
			int noOfProcesses, int rank, double heat) {
		this(tree, pars, ppm);
		this.noOfProcesses = noOfProcesses;
		this.rank = rank;
		this.tree.heat = heat;

		// Is parallel!
		isParallel = true;
	}

	private int alignmentSampled = 0;
	private int alignmentAccepted = 0;
	private int edgeSampled = 0;
	private int edgeAccepted = 0;
	private int topologySampled = 0;
	private int topologyAccepted = 0;
	private int indelSampled = 0;
	private int indelAccepted = 0;
	private int substSampled = 0;
	private int substAccepted = 0;

	/**
	 * In effect starts an MCMC run. It first performs a prescribed number of
	 * burn-in steps, then makes the prescribed number of steps after burn-in,
	 * drawing samples with the prescribes frequency. It also calls the
	 * appropriate functions of the PostpocessManager <tt>postprocMan</tt> to
	 * trigger data transfer to postprocessing modules when necessary
	 */
	public void doMCMC() {
		if (isParallel) {
			String str = String.format(
					"Starting MCMC chain no. %d/%d (heat: %.2f)\n\n", 
					rank + 1, noOfProcesses, tree.heat);
			MPIUtils.println(rank, str);
			swapGenerator = new Random(mcmcpars.swapSeed);
		} else {
			System.out.println("Starting MCMC...\n");
		}

		MainFrame frame = postprocMan.mainManager.frame;
		Utils.generator = new Random(mcmcpars.seed + rank);
		long currentTime, start = System.currentTimeMillis();

		// Triggers a /before first sample/ of the plugins.
		if ((isParallel && MPIUtils.isMaster(rank)) || !isParallel) {
			postprocMan.beforeFirstSample();
		}

		ArrayList<Double> logLikeList = new ArrayList<Double>();

		try {
			//only to use if AutomateParameters.shouldAutomate() == true
			ArrayList<String[]> alignmentsFromSamples = new ArrayList<String[]>(); 
			int burnIn = mcmcpars.burnIn;
			boolean stopBurnIn = false;


			if(AutomateParameters.shouldAutomateBurnIn()){
				burnIn = 10000000;
			} 

			if(AutomateParameters.shouldAutomateStepRate()){
				burnIn += BURNIN_TO_CALCULATE_THE_SPACE;
			}
			//final int SAMPLE_RATE_BURNIN = ((int)((burnIn-burnIn * BURNIN_RATIO_OF_ALIGNM_SPACE_DETERM)/NUMBER_OF_ALIGNMENTS_APPROX));

			burnin = true;
			for (int i = 0; i < burnIn; i++) {

				sample(0);

				// Triggers a /new step/ and a /new peek/ (if appropriate) of
				// the plugins.
				if ((isParallel && MPIUtils.isMaster(rank)) || !isParallel) {
					// TODO do above inside sample() and add more info
					mcmcStep.newLogLike = tree.getLogLike();
					mcmcStep.burnIn = burnin;
					postprocMan.newStep(mcmcStep);
					if (i % mcmcpars.sampRate == 0) {
						postprocMan.newPeek();
					}
				}

				if(AutomateParameters.shouldAutomateBurnIn() && i % 50 == 0){
					logLikeList.add(getState().logLike);
					if(!stopBurnIn){
						stopBurnIn = AutomateParameters.shouldStopBurnIn(logLikeList);
						if(AutomateParameters.shouldAutomateStepRate() && stopBurnIn){
							burnIn = i + BURNIN_TO_CALCULATE_THE_SPACE;
						}else if (stopBurnIn){
							burnIn = i;
						}
					}
				}
				currentTime = System.currentTimeMillis();
				if (frame != null) {
					String text = "";
					if((i > burnIn - BURNIN_TO_CALCULATE_THE_SPACE ) && AutomateParameters.shouldAutomateStepRate()){
						text = "Burn In to get the space: " + (i-(burnIn - BURNIN_TO_CALCULATE_THE_SPACE) + 1) ;
					}else{
						text = "Burn In: " + (i + 1);
					}
					frame.statusText.setText(text);
				} else if (i % 1000 == 999) {
					System.out.println("Burn in: " + (i + 1));
				}


				if( AutomateParameters.shouldAutomateStepRate() && (i >= burnIn - BURNIN_TO_CALCULATE_THE_SPACE) && i % SAMPLE_RATE_WHEN_DETERMINING_THE_SPACE == 0)   {
					String[] align = getState().getLeafAlign();
					alignmentsFromSamples.add(align);
				}	
			}
			
			burnin = false;

			int period;
			if(AutomateParameters.shouldAutomateNumberOfSamples()){
				period = 1000000;
			}else{
				period = mcmcpars.cycles / mcmcpars.sampRate;
			}

			int sampRate;
			if(AutomateParameters.shouldAutomateStepRate()){
				frame.statusText.setText("Calculating the sample rate");
				ArrayList<Double> theSpace = Distance.spaceAMA(alignmentsFromSamples);
				sampRate = AutomateParameters.getSampleRateOfTheSpace(theSpace,SAMPLE_RATE_WHEN_DETERMINING_THE_SPACE);

			}else{
				sampRate = mcmcpars.sampRate;
			}


			int swapNo = 0; // TODO: delete?
			swapCounter = mcmcpars.swapRate;
			AlignmentData alignment = new AlignmentData(getState().getLeafAlign());
			ArrayList<AlignmentData> allAlignments = new ArrayList<AlignmentData>();
			ArrayList<Double> distances = new ArrayList<Double>();

			boolean shouldStop = false;
			double currScore = 0;
			for (int i = 0; i < period && !shouldStop; i++) {
				for (int j = 0; j < sampRate; j++) {
					// Samples.
					sample(0);

					//FuzzyAlignment fuzzyAlignment2 = FuzzyAlignment.getFuzzyAlignmentAndProject(alignments, "");

					// Proposes a swap.
					if (isParallel) {
						swapCounter--;
						if (swapCounter == 0) {
							swapNo++;
							swapCounter = mcmcpars.swapRate;

							doSwap(swapNo);
						}
					}

					// Triggers a /new step/ and a /new peek/ (if appropriate)
					// of the plugins.
					if ((isParallel && MPIUtils.isMaster(rank)) || !isParallel) {
						mcmcStep.newLogLike = tree.getLogLike();
						mcmcStep.burnIn = burnin;
						postprocMan.newStep(mcmcStep);
						if (burnIn + i * period + j % mcmcpars.sampRate == 0) {
							postprocMan.newPeek();
						}
					}

					currentTime = System.currentTimeMillis();
					if (frame != null) {
						String text = "Samples taken: " + Integer.toString(i);
						//remainingTime((currentTime - start)
						//		* ((period - i - 1) * sampRate
						//				+ sampRate - j - 1)
						//				/ (burnIn + i * sampRate + j + 1))

						text += "   The sampling rate: " + sampRate;
						if(AutomateParameters.shouldAutomateNumberOfSamples()){
							text +=  "  Similarity between the last two alignments: " + currScore;
						}
						frame.statusText.setText(text );
					}
				}
				if (frame == null && !isParallel) {
					System.out.println("Sample: " + (i + 1));
				}
				if(AutomateParameters.shouldAutomateNumberOfSamples()){
					alignment = new AlignmentData(getState().getLeafAlign());
					allAlignments.add(alignment);
					if (allAlignments.size() >1){
						FuzzyAlignment Fa = FuzzyAlignment.getFuzzyAlignmentAndProject(allAlignments.subList(0, allAlignments.size()-1), 0);
						FuzzyAlignment Fb = FuzzyAlignment.getFuzzyAlignmentAndProject(allAlignments, 0);
						//System.out.println(Fa);
						//System.out.println("xxxx");
						//System.out.println(Fb);
						currScore = FuzzyAlignment.AMA(Fa, Fb);
						System.out.println(currScore);
						distances.add(currScore);
						if (allAlignments.size() >5){
							shouldStop = AutomateParameters.shouldStopSampling(distances);
						}

					}
				}
				// Report the results of the sample.
				report(i, period);
			}
		} catch (StoppedException ex) {
			// stopped: report and save state
			// should we still call afterLastSample?
		}

		if(Utils.DEBUG) {
			System.out.println("Times spent in each MCMC step type (ms):");
			System.out.println(ali);
			System.out.println(top);
			System.out.println(edge);
			System.out.println(ind);
			System.out.println(sub);
		}

		// Triggers a /after first sample/ of the plugins.
		if ((isParallel && MPIUtils.isMaster(rank)) || !isParallel) {
			postprocMan.afterLastSample();
		}
		if (frame != null) {
			frame.statusText.setText(MainFrame.IDLE_STATUS_MESSAGE);
		}

	}

	private void doSwap(int swapNo) {
		int swapA, swapB;
		swapA = swapGenerator.nextInt(noOfProcesses);
		do {
			swapB = swapGenerator.nextInt(noOfProcesses);
		} while (swapA == swapB);

		System.out.printf("SwapNo: %d - SwapA: %d - SwapB: %d\n", swapNo,
				swapA, swapB);

		double swapAccept = swapGenerator.nextDouble();

		if (rank == swapA || rank == swapB) {
			double[] myStateInfo = new double[3];
			myStateInfo[0] = tree.getLogLike();
			myStateInfo[1] = tree.getLogPrior();
			myStateInfo[2] = tree.heat;

			double[] partnerStateInfo = new double[3];

			mpi.Request send, recieve;

			if (rank == swapA) {
				send = MPI.COMM_WORLD.Isend(myStateInfo, 0, 3, MPI.DOUBLE,
						swapB, 0);
				recieve = MPI.COMM_WORLD.Irecv(partnerStateInfo, 0, 3,
						MPI.DOUBLE, swapB, 1);
			} else {
				send = MPI.COMM_WORLD.Isend(myStateInfo, 0, 3, MPI.DOUBLE,
						swapA, 1);
				recieve = MPI.COMM_WORLD.Irecv(partnerStateInfo, 0, 3,
						MPI.DOUBLE, swapA, 0);
			}

			mpi.Request.Waitall(new mpi.Request[] { send, recieve });

			System.out
			.printf("[Worker %d] Heat: [%f] - Sent: [%f,%f,%f] - Recv: [%f,%f,%f]\n",
					rank, tree.heat, myStateInfo[0], myStateInfo[1],
					myStateInfo[2], partnerStateInfo[0],
					partnerStateInfo[1], partnerStateInfo[2]);

			double myLogLike = myStateInfo[0];
			double myLogPrior = myStateInfo[1];
			double myTemp = myStateInfo[2];
			double hisLogLike = partnerStateInfo[0];
			double hisLogPrior = partnerStateInfo[1];
			double hisTemp = partnerStateInfo[2];

			double acceptance = myTemp * (hisLogLike + hisLogPrior) + hisTemp
					* (myLogLike + myLogPrior);
			acceptance -= hisTemp * (hisLogLike + hisLogPrior) + myTemp
					* (myLogLike + myLogPrior);

			MPIUtils.println(rank,
					"Math.log(swapAccept): " + Math.log(swapAccept));
			MPIUtils.println(rank, "acceptance:           "
					+ acceptance);

			if (acceptance > Math.log(swapAccept)) {
				MPIUtils.println(rank,
						"Just swapped heat with my partner. New heat: "
								+ hisTemp);
				tree.heat = hisTemp;
			}

			// MPI.COMM_WORLD.Send(myStateInfo, 0, 3, MPI.DOUBLE,
			// swapB, 0);
			// statalign.Utils.printLine(swapA, "Just sent " + swapB
			// + " my state.");
		}

	}

	private static String remainingTime(long x) {
		x /= 1000;
		return String.format("Estimated time left: %d:%02d:%02d", x / 3600,
				(x / 60) % 60, x % 60);
	}

	SimpleStats ali;
	SimpleStats top;
	SimpleStats edge;
	SimpleStats ind;
	SimpleStats sub;

	{
		if(Utils.DEBUG) {
			ali = new SimpleStats("Alignment");
			top = new SimpleStats("Topology");
			edge = new SimpleStats("Edge len");
			ind = new SimpleStats("Indel param");
			sub = new SimpleStats("Subst param");
		}
	}

	private void sample(int samplingMethod) throws StoppedException {
		if (samplingMethod == 0) {
			long timer;
			stoppable();
			switch (tree.substitutionModel.params != null
					&& tree.substitutionModel.params.length > 0 ? Utils
							.weightedChoose(FIVECHOOSE) : Utils
							.weightedChoose(FOURCHOOSE)) {
							case 0:
								if(Utils.DEBUG)
									timer = -System.currentTimeMillis();
								sampleEdge();
								if(Utils.DEBUG) {
									timer += System.currentTimeMillis();
									edge.addData(timer);
								}
								break;
							case 1:
								if(Utils.DEBUG)
									timer = -System.currentTimeMillis();
								sampleTopology();
								if(Utils.DEBUG) {
									timer += System.currentTimeMillis();
									top.addData(timer);
								}
								break;
							case 2:
								if(Utils.DEBUG)
									timer = -System.currentTimeMillis();
								sampleIndelParameter();
								if(Utils.DEBUG) {
									timer += System.currentTimeMillis();
									ind.addData(timer);
								}
								break;
							case 3:
								if(Utils.DEBUG)
									timer = -System.currentTimeMillis();
								sampleAlignment();
								if(Utils.DEBUG) {
									timer += System.currentTimeMillis();
									ali.addData(timer);
								}
								break;
							case 4:
								if(Utils.DEBUG)
									timer = -System.currentTimeMillis();
								sampleSubstParameter();
								if(Utils.DEBUG) {
									timer += System.currentTimeMillis();
									sub.addData(timer);
								}
								break;
			}
		} else {
			stoppable();
			sampleEdge();
			sampleTopology();
			sampleIndelParameter();
			sampleSubstParameter();
			sampleAlignment();
		}
	}

	private void sampleEdge() {
		edgeSampled++;
		// System.out.print("Edge: ");
		int i = Utils.generator.nextInt(tree.vertex.length - 1);
		double oldEdge = tree.vertex[i].edgeLength;
		double oldLogLikelihood = tree.getLogLike();
		while ((tree.vertex[i].edgeLength = oldEdge
				+ Utils.generator.nextDouble() * Utils.EDGE_SPAN
				- (Utils.EDGE_SPAN / 2.0)) < 0.01)
			;
		tree.vertex[i].edgeChangeUpdate();
		// Vertex actual = tree.vertex[i];
		// while(actual != null){
		// actual.calcFelsen();
		// actual.calcOrphan();
		// actual.calcIndelLogLike();
		// actual = actual.parent;
		// }
		tree.vertex[i].calcAllUp();
		double newLogLikelihood = tree.getLogLike();
		if (Utils.generator.nextDouble() < (Math.exp((newLogLikelihood
				- oldLogLikelihood - tree.vertex[i].edgeLength + oldEdge)
				* tree.heat) * (Math.min(oldEdge - 0.01, Utils.EDGE_SPAN / 2.0) + Utils.EDGE_SPAN / 2.0))
				/ (Math.min(tree.vertex[i].edgeLength - 0.01,
						Utils.EDGE_SPAN / 2.0) + Utils.EDGE_SPAN / 2.0)) {
			// acceptance, do nothing
			// System.out.println("accepted (old: "+oldLogLikelihood+" new: "+newLogLikelihood+")");
			edgeAccepted++;
		} else {
			// reject, restore
			// System.out.print("Rejected! i: "+i+"\tOld likelihood: "+oldLogLikelihood+"\tNew likelihood: "+newLogLikelihood);
			tree.vertex[i].edgeLength = oldEdge;
			tree.vertex[i].edgeChangeUpdate();
			// actual = tree.vertex[i];
			// while(actual != null){
			// actual.calcFelsen();
			// actual.calcOrphan();
			// actual.calcIndelLogLike();
			// actual = actual.parent;
			// }
			tree.vertex[i].calcAllUp();
			// System.out.println("rejected (old: "+oldLogLikelihood+" new: "+newLogLikelihood+")");

		}
	}

	// this is the old
	/*
	 * private void sampleTopology(){ int vnum = tree.vertex.length;
	 * 
	 * if(vnum <= 3) return;
	 * 
	 * System.out.print("Topology: "); double oldLogLi = tree.getLogLike();
	 * 
	 * int vertId, rnd = Utils.generator.nextInt(vnum-3); vertId =
	 * tree.getTopVertexId(rnd); if(vertId != -1) { int lastId[] = new int[3],
	 * num = 0, newId = vertId;
	 * 
	 * for(int i = vnum-3; i < vnum; i++) { int id = tree.getTopVertexId(i);
	 * if(id == -1) lastId[num++] = i; else if(id < vertId) newId--; } rnd =
	 * lastId[newId]; } Vertex nephew = tree.vertex[rnd]; Vertex uncle =
	 * nephew.parent.brother();
	 * 
	 * // for(vertId = 0; vertId < vnum; vertId++) { //
	 * if(tree.getTopVertexId(vertId) == -1) { // vertex eligible // if(rnd-- ==
	 * 0) // break; // } // } // Vertex nephew = tree.vertex[vertId];
	 * 
	 * double bpp = nephew.swapWithUncle();
	 * 
	 * double newLogLi = tree.getLogLike();
	 * 
	 * // tree.root.calcFelsRecursivelyWithCheck();
	 * //tree.root.calcIndelRecursivelyWithCheck();
	 * 
	 * if(Math.log(Utils.generator.nextDouble()) < bpp+newLogLi-oldLogLi) { //
	 * accepted
	 * System.out.println("accepted (old: "+oldLogLi+" new: "+newLogLi+")"); }
	 * else { // refused uncle.swapBackUncle();
	 * System.out.println("rejected (old: "+oldLogLi+" new: "+newLogLi+")"); }
	 * 
	 * //tree.root.calcFelsRecursivelyWithCheck();
	 * //tree.root.calcIndelRecursivelyWithCheck(); }
	 */
	private void sampleTopology() {
		int vnum = tree.vertex.length;

		if (vnum <= 3)
			return;

		topologySampled++;
		// System.out.println("\n\n\t***\t***\t***\n\n\n");
		// System.out.print("Topology: ");
		// tree.printAllPointers();
		double oldLogLi = tree.getLogLike();

		int vertId, rnd = Utils.generator.nextInt(vnum - 3);
		vertId = tree.getTopVertexId(rnd);
		if (vertId != -1) {
			int lastId[] = new int[3], num = 0, newId = vertId;

			for (int i = vnum - 3; i < vnum; i++) {
				int id = tree.getTopVertexId(i);
				if (id == -1)
					lastId[num++] = i;
				else if (id < vertId)
					newId--;
			}
			rnd = lastId[newId];
		}
		Vertex nephew = tree.vertex[rnd];
		Vertex uncle = nephew.parent.brother();

		// for(vertId = 0; vertId < vnum; vertId++) {
		// if(tree.getTopVertexId(vertId) == -1) { // vertex eligible
		// if(rnd-- == 0)
		// break;
		// }
		// }
		// Vertex nephew = tree.vertex[vertId];

		// String[] s = tree.root.printedMultipleAlignment();
		// System.out.println("Alignment before topology changing: ");
		// for(int i = 0; i < s.length; i++){
		// System.out.println(s[i]);
		// }
		double bpp = nephew.fastSwapWithUncle();
		// double bpp = nephew.swapWithUncle();
		// s = tree.root.printedMultipleAlignment();
		// System.out.println("Alignment after topology changing: ");
		// for(int i = 0; i < s.length; i++){
		// System.out.println(s[i]);
		// }

		double newLogLi = tree.getLogLike();

		// tree.root.calcFelsRecursivelyWithCheck();
		// tree.root.calcIndelRecursivelyWithCheck();

		if (Math.log(Utils.generator.nextDouble()) < bpp
				+ (newLogLi - oldLogLi) * tree.heat) {
			// accepted
			// System.out.println("accepted (old: "+oldLogLi+" new: "+newLogLi+")");
			topologyAccepted++;
		} else {
			// rejected
			// System.out.println("Checking pointer integrity before changing back topology: ");
			for (int i = 0; i < tree.vertex.length; i++) {
				if (tree.vertex[i].left != null && tree.vertex[i].right != null) {
					tree.vertex[i].checkPointers();
					AlignColumn p;
					// checking pointer integrity
					for (AlignColumn c = tree.vertex[i].left.first; c != null; c = c.next) {
						p = tree.vertex[i].first;
						while (c.parent != p && p != null)
							p = p.next;
						if (p == null)
							throw new Error(
									"children does not have a parent!!!"
											+ tree.vertex[i] + " "
											+ tree.vertex[i].print());
					}
					for (AlignColumn c = tree.vertex[i].right.first; c != null; c = c.next) {
						p = tree.vertex[i].first;
						while (c.parent != p && p != null)
							p = p.next;
						if (p == null)
							throw new Error(
									"children does not have a parent!!!"
											+ tree.vertex[i] + " "
											+ tree.vertex[i].print());
					}

				}
			}

			uncle.fastSwapBackUncle();
			// System.out.println("Checking pointer integrity after changing back topology: ");
			for (int i = 0; i < tree.vertex.length; i++) {
				if (tree.vertex[i].left != null && tree.vertex[i].right != null) {
					tree.vertex[i].checkPointers();
					AlignColumn p;
					// checking pointer integrity
					for (AlignColumn c = tree.vertex[i].left.first; c != null; c = c.next) {
						p = tree.vertex[i].first;
						while (c.parent != p && p != null)
							p = p.next;
						if (p == null)
							throw new Error(
									"children does not have a parent!!!"
											+ tree.vertex[i] + " "
											+ tree.vertex[i].print());
					}
					for (AlignColumn c = tree.vertex[i].right.first; c != null; c = c.next) {
						p = tree.vertex[i].first;
						while (c.parent != p && p != null)
							p = p.next;
						if (p == null)
							throw new Error(
									"children does not have a parent!!!"
											+ tree.vertex[i] + " "
											+ tree.vertex[i].print());
					}
				}
			}
			// uncle.swapBackUncle();
			// s = tree.root.printedMultipleAlignment();
			// System.out.println("Alignment after changing back the topology: ");
			// for(int i = 0; i < s.length; i++){
			// System.out.println(s[i]);
			// }
			// System.out.println("rejected (old: "+oldLogLi+" new: "+newLogLi+")");
		}

		// tree.printAllPointers();
		// System.out.println("\n\n\t***\t***\t***\n\n\n");
		tree.root.calcFelsRecursivelyWithCheck();
		tree.root.calcIndelRecursivelyWithCheck();
	}

	private void sampleIndelParameter() {
		indelSampled++;
		switch (Utils.generator.nextInt(3)) {
		case 0:
			// System.out.print("Indel param R: ");
			double oldR = tree.hmm2.params[0];
			double oldLogLikelihood = tree.root.orphanLogLike
					+ tree.root.indelLogLike;
			while ((tree.hmm2.params[0] = oldR + Utils.generator.nextDouble()
					* Utils.R_SPAN - Utils.R_SPAN / 2.0) <= 0.0
					|| tree.hmm2.params[0] >= 1.0)
				;
			for (int i = 0; i < tree.vertex.length; i++) {
				tree.vertex[i].updateHmmMatrices();
			}
			tree.root.calcIndelLikeRecursively();
			double newLogLikelihood = tree.root.orphanLogLike
					+ tree.root.indelLogLike;
			if (Utils.generator.nextDouble() < Math
					.exp((newLogLikelihood - oldLogLikelihood) * tree.heat)
					* (Math.min(1.0 - oldR, Utils.R_SPAN / 2.0) + Math.min(
							oldR, Utils.R_SPAN / 2.0))
							/ (Math.min(1.0 - tree.hmm2.params[0], Utils.R_SPAN / 2.0) + Math
									.min(tree.hmm2.params[0], Utils.R_SPAN / 2.0))) {
				// accept, do nothing
				// System.out.println("accepted (old: "+oldLogLikelihood+" new: "+newLogLikelihood+")");
				indelAccepted++;
			} else {
				// restore
				tree.hmm2.params[0] = oldR;
				for (int i = 0; i < tree.vertex.length; i++) {
					tree.vertex[i].updateHmmMatrices();
				}
				tree.root.calcIndelLikeRecursively();
				// System.out.println("rejected (old: "+oldLogLikelihood+" new: "+newLogLikelihood+")");
			}

			break;
		case 1:
			// ///////////////////////////////////////////////
			// System.out.print("Indel param Lambda: ");
			double oldLambda = tree.hmm2.params[1];
			oldLogLikelihood = tree.root.orphanLogLike + tree.root.indelLogLike;
			while ((tree.hmm2.params[1] = oldLambda
					+ Utils.generator.nextDouble() * Utils.LAMBDA_SPAN
					- Utils.LAMBDA_SPAN / 2.0) <= 0.0
					|| tree.hmm2.params[1] >= tree.hmm2.params[2])
				;
			for (int i = 0; i < tree.vertex.length; i++) {
				tree.vertex[i].updateHmmMatrices();
			}
			tree.root.calcIndelLikeRecursively();
			newLogLikelihood = tree.root.orphanLogLike + tree.root.indelLogLike;
			if (Utils.generator.nextDouble() < Math.exp((newLogLikelihood
					- oldLogLikelihood - tree.hmm2.params[1] + oldLambda)
					* tree.heat)
					* (Math.min(Utils.LAMBDA_SPAN / 2.0, tree.hmm2.params[2]
							- oldLambda) + Math.min(oldLambda,
									Utils.LAMBDA_SPAN / 2.0))
									/ (Math.min(Utils.LAMBDA_SPAN / 2.0, tree.hmm2.params[2]
											- tree.hmm2.params[1]) + Math.min(
													tree.hmm2.params[1], Utils.LAMBDA_SPAN / 2.0))) {
				// accept, do nothing
				// System.out.println("accepted (old: "+oldLogLikelihood+" new: "+newLogLikelihood+" oldLambda: "+oldLambda+" newLambda: "+tree.hmm2.params[1]+")");
				indelAccepted++;
			} else {
				// restore
				tree.hmm2.params[1] = oldLambda;
				for (int i = 0; i < tree.vertex.length; i++) {
					tree.vertex[i].updateHmmMatrices();
				}
				tree.root.calcIndelLikeRecursively();
				// System.out.println("rejected (old: "+oldLogLikelihood+" new: "+newLogLikelihood+" oldLambda: "+oldLambda+" newLambda: "+tree.hmm2.params[1]+")");
			}
			break;
		case 2:
			// ///////////////////////////////////////////////////////
			// System.out.print("Indel param Mu: ");
			double oldMu = tree.hmm2.params[2];
			oldLogLikelihood = tree.getLogLike();
			while ((tree.hmm2.params[2] = oldMu + Utils.generator.nextDouble()
					* Utils.MU_SPAN - Utils.MU_SPAN / 2.0) <= tree.hmm2.params[1])
				;
			for (int i = 0; i < tree.vertex.length; i++) {
				tree.vertex[i].updateHmmMatrices();
			}
			tree.root.calcIndelLikeRecursively();
			newLogLikelihood = tree.getLogLike();
			if (Utils.generator.nextDouble() < Math.exp((newLogLikelihood
					- oldLogLikelihood - tree.hmm2.params[2] + oldMu)
					* tree.heat)
					* (Utils.MU_SPAN / 2.0 + Math.min(oldMu
							- tree.hmm2.params[1], Utils.MU_SPAN / 2.0))
							/ (Utils.MU_SPAN / 2.0 + Math.min(tree.hmm2.params[2]
									- tree.hmm2.params[1], Utils.MU_SPAN / 2.0))) {
				// accept, do nothing
				// System.out.println("accepted (old: "+oldLogLikelihood+" new: "+newLogLikelihood+")");
				indelAccepted++;
			} else {
				// restore
				tree.hmm2.params[2] = oldMu;
				for (int i = 0; i < tree.vertex.length; i++) {
					tree.vertex[i].updateHmmMatrices();
				}
				tree.root.calcIndelLikeRecursively();
				// System.out.println("rejected (old: "+oldLogLikelihood+" new: "+newLogLikelihood+")");
			}
			break;
		}
	}

	private void sampleSubstParameter() {
		substSampled++;
		if (tree.substitutionModel.params.length == 0)
			return;
		else {
			double mh = tree.substitutionModel.sampleParameter();
			double oldlikelihood = tree.root.orphanLogLike;
			for (int i = 0; i < tree.vertex.length; i++) {
				tree.vertex[i].updateTransitionMatrix();
			}
			tree.root.calcFelsRecursively();
			double newlikelihood = tree.root.orphanLogLike;
			if (Utils.generator.nextDouble() < Math.exp(mh
					+ (Math.log(tree.substitutionModel.getPrior())
							+ newlikelihood - oldlikelihood))
							* tree.heat) {
				// System.out.println("Substitution parameter: accepted (old: "+oldlikelihood+" new: "+newlikelihood+")");
				substAccepted++;
			} else {
				tree.substitutionModel.restoreParameter();
				for (int i = 0; i < tree.vertex.length; i++) {
					tree.vertex[i].updateTransitionMatrix();
				}
				tree.root.calcFelsRecursively();
				// System.out.println("Substitution parameter: rejected (old: "+oldlikelihood+" new: "+newlikelihood+")");
			}
		}
	}

	private void sampleAlignment() {
		alignmentSampled++;
		for (int i = 0; i < tree.vertex.length; i++) {
			tree.vertex[i].selected = false;
		}
		// System.out.print("Alignment: ");
		double oldLogLi = tree.getLogLike();
		// System.out.println("fast indel before: "+tree.root.indelLogLike);
		tree.countLeaves(); // calculates recursively how many leaves we have
		// below this node
		for (int i = 0; i < weights.length; i++) {
			weights[i] = Math.pow(tree.vertex[i].leafCount, LEAFCOUNT_POWER);
		}
		int k = Utils.weightedChoose(weights, null);
		// System.out.println("Sampling from the subtree: "+tree.vertex[k].print());
		tree.vertex[k].selectSubtree(SELTRLEVPROB, 0);
		double bpp = tree.vertex[k].selectAndResampleAlignment();
		double newLogLi = tree.getLogLike();

		// String[] printedAlignment = tree.printedAlignment("StatAlign");
		// for(String i: printedAlignment)
		// System.out.println(i);
		//
		// System.out.println("-----------------------------------------------------------------------------");
		// double fastFels = tree.root.orphanLogLike;
		// double fastIns = tree.root.indelLogLike;
		// report();
		// tree.root.first.seq[0] = 0.0;
		// System.out.println("Old before: "+tree.root.old.indelLogLike);
		// tree.root.calcFelsRecursivelyWithCheck();
		// tree.root.calcIndelRecursivelyWithCheck();
		// tree.root.calcIndelLikeRecursively();
		// System.out.println("Old after: "+tree.root.old.indelLogLike);
		// System.out.println("Check logli: "+tree.getLogLike()+" fastFels: "+fastFels+" slowFels: "+tree.root.orphanLogLike+
		// " fastIns: "+fastIns+" slowIns: "+tree.root.indelLogLike);
		// System.out.println("selected subtree: "+tree.vertex[k].print());
		// System.out.println("bpp: "+bpp+"old: "+oldLogLi+"new: "+newLogLi +
		// "heated diff: " + ((newLogLi - oldLogLi) * tree.heat));
		if (Math.log(Utils.generator.nextDouble()) < bpp
				+ (newLogLi - oldLogLi) * tree.heat) {
			// accepted
			// System.out.println("accepted (old: "+oldLogLi+" new: "+newLogLi+")");
			alignmentAccepted++;
		} else {
			// refused
			// String[] s = tree.printedAlignment();
			tree.vertex[k].alignRestore();
			// s = tree.printedAlignment();
			// System.out.println("rejected (old: "+oldLogLi+" new: "+newLogLi+")");
			// System.out.println("after reject fast: "+tree.root.indelLogLike);
			// tree.root.calcIndelRecursivelyWithCheck();
			// System.out.println(" slow: "+tree.root.indelLogLike);
		}
		// tree.root.calcFelsRecursivelyWithCheck();
		// tree.root.calcIndelRecursivelyWithCheck();
	}

	/**
	 * Returns a string representation describing the acceptance ratios of the current MCMC run.
	 * @return a string describing the acceptance ratios.
	 */
	public String getInfoString() {
		return String.format("Acceptances: [Alignment: %f, Edge: %f, Topology: %f, Indel: %f, Substitution: %f]",
				(alignmentSampled == 0 ? 0 : (double) alignmentAccepted
						/ (double) alignmentSampled),
						(edgeSampled == 0 ? 0 : (double) edgeAccepted
								/ (double) edgeSampled),
								(topologySampled == 0 ? 0 : (double) topologyAccepted
										/ (double) topologySampled),
										(indelSampled == 0 ? 0 : (double) indelAccepted
												/ (double) indelSampled),
												(substSampled == 0 ? 0 : (double) substAccepted
														/ (double) substSampled));
	}

	/**
	 * Returns a {@link State} object that describes the current state of the
	 * MCMC. This can then be passed on to other classes such as postprocessing
	 * plugins.
	 */
	public State getState() {
		return tree.getState();
	}

	private boolean isColdChain() {
		return tree.heat == 1.0d;
	}

	private State MPIStateReceieve(int peer) {
		// Creates a new, uninitialized state and initializes the variables.
		State state = new State(tree.vertex.length);

		// We already know the names
		for (int i = 0; i < state.nl; i++) {
			state.name[i] = tree.vertex[i].name;
		}

		int nn = state.nn;
		int tag = 0;

		// left
		MPI.COMM_WORLD.Recv(state.left, 0, nn, MPI.INT, peer, tag++);
		// right
		MPI.COMM_WORLD.Recv(state.right, 0, nn, MPI.INT, peer, tag++);
		// parent
		MPI.COMM_WORLD.Recv(state.parent, 0, nn, MPI.INT, peer, tag++);
		// edgeLen
		MPI.COMM_WORLD.Recv(state.edgeLen, 0, nn, MPI.DOUBLE, peer, tag++);

		// sequences
		int[] seqLengths = new int[nn];
		MPI.COMM_WORLD.Recv(seqLengths, 0, nn, MPI.INT, peer, tag++);

		for (int i = 0; i < nn; i++) {
			char[] c = new char[seqLengths[i]];
			MPI.COMM_WORLD.Recv(c, 0, seqLengths[i], MPI.CHAR, peer, tag++);
			state.seq[i] = new String(c);
		}

		// align
		Object[] recvObj = new Object[1];
		MPI.COMM_WORLD.Recv(recvObj, 0, 1, MPI.OBJECT, peer, tag++);
		state.align = (int[][]) recvObj[0];

		// felsen
		MPI.COMM_WORLD.Recv(recvObj, 0, 1, MPI.OBJECT, peer, tag++);
		state.felsen = (double[][][]) recvObj[0];

		// indelParams
		final int noOfIndelParameter = 3;
		state.indelParams = new double[noOfIndelParameter];
		MPI.COMM_WORLD.Recv(state.indelParams, 0, noOfIndelParameter,
				MPI.DOUBLE, peer, tag++);

		// substParams
		int l = tree.substitutionModel.params.length;
		state.substParams = new double[l];
		MPI.COMM_WORLD.Recv(state.substParams, 0, l, MPI.DOUBLE, peer, tag++);

		// log-likelihood
		double[] d = new double[1];
		MPI.COMM_WORLD.Recv(d, 0, 1, MPI.DOUBLE, peer, tag++);
		state.logLike = d[0];

		// root
		int[] root = new int[1];
		MPI.COMM_WORLD.Recv(root, 0, 1, MPI.INT, peer, tag++);
		state.root = root[0];

		return state;
	}

	private void MPIStateSend(State state) {

		String[] seq = state.seq;
		int[][] align = state.align;
		double[][][] felsen = state.felsen;
		int nn = state.nn;
		int tag = 0;

		// left
		MPI.COMM_WORLD.Send(state.left, 0, nn, MPI.INT, 0, tag++);
		// right
		MPI.COMM_WORLD.Send(state.right, 0, nn, MPI.INT, 0, tag++);
		// parent
		MPI.COMM_WORLD.Send(state.parent, 0, nn, MPI.INT, 0, tag++);
		// edgeLen
		MPI.COMM_WORLD.Send(state.edgeLen, 0, nn, MPI.DOUBLE, 0, tag++);

		// TODO: START OF OPTIMIZATION.

		// sequences
		int[] seqLength = new int[nn];
		char[][] seqChars = new char[nn][];
		for (int i = 0; i < nn; i++) {
			seqLength[i] = seq[i].length();
			seqChars[i] = seq[i].toCharArray();
		}
		MPI.COMM_WORLD.Send(seqLength, 0, nn, MPI.INT, 0, tag++);
		for (int i = 0; i < nn; i++) {
			MPI.COMM_WORLD.Send(seqChars[i], 0, seqLength[i], MPI.CHAR, 0, tag++);
		}

		// align
		Object[] alignObj = new Object[1];
		alignObj[0] = align;
		MPI.COMM_WORLD.Send(alignObj, 0, 1, MPI.OBJECT, 0, tag++);
		/*
		 * int[] alignLength = new int[align.length]; for (int i = 0; i <
		 * seq.length; i++) { alignLength[i] = align[i].length; }
		 * MPI.COMM_WORLD.Send(alignLength, 0, nn, MPI.INT, 0, tag++); for (int
		 * i = 0; i < align.length; i++) { MPI.COMM_WORLD.Send(align[i], 0,
		 * alignLength[i], MPI.INT, 0, tag++); }
		 */

		// felsen
		Object[] felsenObj = new Object[] { felsen };
		MPI.COMM_WORLD.Send(felsenObj, 0, 1, MPI.OBJECT, 0, tag++);

		// indelParams
		MPI.COMM_WORLD.Send(state.indelParams, 0, 3, MPI.DOUBLE, 0, tag++);

		// substParams
		MPI.COMM_WORLD.Send(state.substParams, 0, state.substParams.length,
				MPI.DOUBLE, 0, tag++);

		// loglikelihood
		MPI.COMM_WORLD.Send(new double[] { state.logLike }, 0, 1, MPI.DOUBLE,
				0, tag++);

		// root
		MPI.COMM_WORLD.Send(new int[] { state.root }, 0, 1, MPI.INT, 0, tag++);

		// TODO: END OF OPTIMIZATION.

	}

	private void report(int no, int total) {

		int coldChainLocation = -1;

		if (isParallel) {
			// Get rank of cold chain.
			int[] ranks = new int[] { (isColdChain() ? rank : 0) };
			int[] coldChainLoc = new int[1];
			MPI.COMM_WORLD.Reduce(ranks, 0, coldChainLoc, 0, 1, MPI.INT, MPI.SUM, 0);
			coldChainLocation = coldChainLoc[0];

			// TODO: Remove - for debugging purposes
			if (MPIUtils.isMaster(rank)) {
				MPIUtils.println(rank, "Cold chain is at: " + coldChainLocation);
			}

			if (isColdChain() && MPIUtils.isMaster(rank)) {
				// Sample normally.
				postprocMan.newSample(getState(), no, total);
			} else if (isColdChain() && !MPIUtils.isMaster(rank)) {
				// Send state.
				State state = getState();
				MPIStateSend(state);
			} else if (!isColdChain() && MPIUtils.isMaster(rank)) {
				// Receive state.
				State state = MPIStateReceieve(coldChainLocation);
				postprocMan.newSample(state, no, total);
			}

		} else {
			postprocMan.newSample(getState(), no, total);
		}

		// Log the accept ratios/params to the (.log) file. TODO: move to a plugin.
		try {
			if ((isParallel && MPIUtils.isMaster(rank)) || !isParallel) {
				postprocMan.logFile.write(getInfoString() + "\n");
				postprocMan.logFile.write("Report\tLogLikelihood\t"
						+ (tree.root.orphanLogLike + tree.root.indelLogLike)
						+ "\tR\t" + tree.hmm2.params[0] + "\tLamda\t"
						+ tree.hmm2.params[1] + "\tMu\t" + tree.hmm2.params[2]
								+ "\t" + tree.substitutionModel.print() + "\n");
				if (isParallel) {
					postprocMan.logFile.write("Cold chain location: " + coldChainLocation + "\n");
				}

			}
		} catch (IOException e) {
			if (postprocMan.mainManager.frame != null) {
				new ErrorMessage(null, e.getLocalizedMessage(), true);
			} else {
				e.printStackTrace(System.out);
			}
		}

		// alignmentSampled = 0;
		// alignmentAccepted = 0;
		// edgeSampled = 0;
		// edgeAccepted = 0;
		// topologySampled = 0;
		// topologyAccepted = 0;
		// indelSampled = 0;
		// indelAccepted = 0;
		// substSampled = 0;
		// substAccepted = 0;

	}



	/**
	 * This function is only for testing and debugging purposes.
	 * 
	 * @param args
	 *            Actually, we do not use these parameters, as this function is
	 *            for testing and debugging purposes. All necessary input data
	 *            is directly written into the function.
	 * 
	 */
	// public static void main(String[] args) {
	// try {
	// Tree tree = new Tree(new String[] { "kkkkkkwwwwwwwwlidwwwwwkkk",
	// "kkkwwwwwwwlidwwwwwkkk", "kkkwwwwwwwlidwwwwwkkk",
	// "kkkwwwwwwwlidwwwwwkkk", "kkkwwwwwlidwwwwwkkkddkldkl",
	// "kkkwwwwwlidwwwwwkkkeqiqii", "kkkwwwwwlidwwwwwkkkddkidkil",
	// "kkkwwwwwlidwwwwwkkkeqiq", "kkkwwwwwlidwwwwwkkkddkldkll",
	// "kkkwwwwwlidwwwwwkkkddkldkil" }, new String[] { "A", "B",
	// "C", "D", "E", "F", "G", "H", "I", "J" }, new Dayhoff(), new Blosum62(),
	// "");
	// for (int i = 0; i < tree.vertex.length; i++) {
	// // tree.vertex[i].edgeLength = 0.1;
	// tree.vertex[i].updateTransitionMatrix();
	// }
	// tree.root.calcFelsRecursively();
	// System.out.println(tree.printedTree());
	// // Mcmc mcmc = new Mcmc(tree, new MCMCPars(0, 10000, 10, 1L), new
	// // PostprocessManager(null));
	// // mcmc.doMCMC();
	// } catch (StoppedException e) {
	// // stopped during tree construction
	// } catch (IOException e) {
	// }
	// }

}