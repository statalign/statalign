package statalign.base;

import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Locale;
import java.util.Random;

import mpi.MPI;

import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;

import statalign.MPIUtils;
import statalign.base.mcmc.AlignmentMove;
import statalign.base.mcmc.AllEdgeMove;
import statalign.base.mcmc.CoreMcmcModule;
import statalign.base.mcmc.EdgeMove;
import statalign.base.mcmc.IndelMove;
import statalign.base.mcmc.LOCALTopologyMove;
import statalign.base.mcmc.LambdaMove;
import statalign.base.mcmc.MuMove;
import statalign.base.mcmc.PhiMove;
import statalign.base.mcmc.RMove;
import statalign.base.mcmc.RhoMove;
import statalign.base.mcmc.SilentIndelMove;
import statalign.base.mcmc.SubstMove;
import statalign.base.mcmc.ThetaMove;
import statalign.base.mcmc.TopologyMove;
import statalign.base.thread.Stoppable;
import statalign.base.thread.StoppedException;
import statalign.mcmc.BetaPrior;
import statalign.mcmc.GammaPrior;
import statalign.mcmc.GaussianProposal;
import statalign.mcmc.LogisticProposal;
import statalign.mcmc.McmcCombinationMove;
import statalign.mcmc.McmcModule;
import statalign.mcmc.McmcMove;
import statalign.mcmc.MultiplicativeProposal;
import statalign.mcmc.UniformProposal;
import statalign.model.ext.ModelExtManager;
import statalign.postprocess.PostprocessManager;
import statalign.postprocess.plugins.contree.CNetwork;
import statalign.ui.ErrorMessage;
import statalign.ui.MainFrame;
import statalign.utils.BetaDistribution;

/**
 * 
 * This class handles an MCMC run.
 * 
 * The class extends <tt>Stoppable</tt>, it may be terminated/suspended in
 * graphical mode.
 * 
 * @author miklos, novak, herman
 * 
 */
public abstract class Mcmc extends Stoppable {

	/** Is this a parallel chain? By-default false. */
	protected boolean isParallel = false;
	
	private boolean simulatedAnnealing = false;

	public CNetwork network; 

	/** Current tree in the MCMC chain. */
	public Tree tree;
	
	protected double heat, heatWhenHot;

	/** Total log-likelihood of the current state, cached for speed */
	protected double totalLogLike;
	
	/** 
	 * If this variable is true, all proposed moves are accepted. 
	 * The purpose of this is to allow:
	 * a) an initial run of moves to randomise the starting configuration
	 * b) a series of moves to be proposed before choosing whether to
	 *    accept or reject the final configuration, with Hastings ratio
	 *    equal to the product of the Hastings ratios along the way.
	 */
	boolean acceptAllMoves = false;
	
	/**
	 * Number of steps in which the chain will be allowed to move randomly
	 * with all moves accepted, in order to create a random starting
	 * configuration.
	 */
	int randomisationPeriod = 0;
	/** 
	 * To be used in order to allow a series of moves to be proposed before
	 * deciding whether to accept or reject the final state, with Hastings
	 * ratio given by this cumulative value.
	 */
	double cumulativeLogProposalRatio = 0.0;

	/**
	 * MCMC parameters including the number of burn-in steps, the total number
	 * of steps in the MCMC and the sampling rate.
	 */
	public MCMCPars mcmcpars;

	public McmcStep mcmcStep = new McmcStep();

	/** PostprocessManager that handles the postprocessing modules. */
	public PostprocessManager postprocMan;

	/** Manager that handles model extension plugins */
	public ModelExtManager modelExtMan;
	
	/** McmcModule containing the moves for the core components of
	 * the model, i.e. the indel parameters, substitution model parameters,
	 * alignment, topology and edge lengths.
	 * The coreModel also decides whether to execute MCMC moves from the 
	 * ModelExtension modules.
	 */
	public McmcModule coreModel;
	
	/** 
	 * Interval (in terms of number of samples) 
	 * at which current postprocessing information is flushed to file
	 * and MCMC info is printed to stdout. 
	 */
	int LOG_INTERVAL = 100; 
	
	/** True while the MCMC is in the burn-in phase. */
	public boolean burnin;
	/** True while the MCMC is in the first half of burn-in phase. */
	public boolean firstHalfBurnin;

	public Mcmc(Tree tree, MCMCPars pars, PostprocessManager ppm, ModelExtManager modelExtMan) {
		postprocMan = ppm;
		this.modelExtMan = modelExtMan;
		ppm.mcmc = this;
		this.modelExtMan.setMcmc(this);
		this.tree = tree;
		
		tree.owner = this;
		mcmcpars = pars;
		heat = 1.0d;
		randomisationPeriod = mcmcpars.randomisationPeriod;
	}	
	
	private static final DecimalFormat df = new DecimalFormat("0.0000");

	/** 
	 * Initialises the core McmcModule. This method is to be implemented by
	 * specific instances of the Mcmc class.
	 * @param tree
	 */
	protected abstract void initCoreModel(Tree tree);	
	
	/**
	 * Triggers a call to <tt>coreModel</tt> to propose a move from one
	 * of the McmcMove objects, selected according to its weight.
	 * If there are active ModelExtension plugins, then <tt>coreModel</tt>
	 * will delegate the sampling to the ModelExtensionManager with
	 * probability proportional to the sum of the weights of the
	 * active plugins.
	 * 
	 * @param samplingMethod Currently unused.
	 * @throws StoppedException
	 */
	private void sample(int samplingMethod) throws StoppedException {
		stoppable();
		if(Utils.DEBUG) {
			System.out.println("tree.getLogLike() (BEFORE) = "+tree.getLogLike());
			tree.recomputeCheckLogLike(); 
			if(Math.abs(modelExtMan.totalLogLike(tree)-totalLogLike) > 1e-5) {
				System.out.println("\nBefore: "+modelExtMan.totalLogLike(tree)+" "+totalLogLike);
				throw new Error("Log-likelihood inconsistency at start of sample()");
			}
		}
		boolean accepted = coreModel.proposeParamChange(tree);
		if (accepted) {
			if (Utils.DEBUG) System.out.println("\t\tMove accepted.");
			totalLogLike = coreModel.curLogLike;
		}
		else {

			if (Utils.DEBUG) System.out.println("Move rejected.");
			coreModel.setLogLike(totalLogLike);
		}
		if(Utils.DEBUG) {
			tree.recomputeCheckLogLike();
        	tree.checkPointers();
			if(Math.abs(modelExtMan.totalLogLike(tree)-totalLogLike) > 1e-5) {
				System.out.println("After: "+modelExtMan.totalLogLike(tree)+" "+totalLogLike);
				throw new Error("Log-likelihood inconsistency at end of sample()");
			}
		}
	}
	
	/**
	 * 
	 * @param aligned Vector indicating which characters are aligned to the current
	 * column in the subtrees below.
	 * @return Logarithm of emission probability for subtrees
	 */
	double calcEm(int[] aligned) {
		return modelExtMan.calcLogEm(aligned);
	}
	/**
	 * This function is called by the McmcMove objects in order to determine whether 
	 * the proposed moves are to be accepted.
	 *   
	 * @param logProposalRatio This also includes the contribution from the prior densities.
	 * It is assumed that any dependencies between the priors and other parameters will be
	 * handled inside the McmcMove objects.
	 * @return true if the move is accepted
	 */
	public boolean isParamChangeAccepted(double logProposalRatio,McmcMove m) {
		double newLogLike = coreModel.curLogLike;
		if (Utils.SHAKE_IF_STUCK && firstHalfBurnin && (m.lowCounts > Utils.LOW_COUNT_THRESHOLD)) {			
			newLogLike *= Math.pow(Utils.LOW_COUNT_MULTIPLIER,m.lowCounts-Utils.LOW_COUNT_THRESHOLD);			
		}
		acceptAllMoves = firstHalfBurnin && m.acceptAllDuringFirstHalfBurnin; 
		boolean accept = acceptanceDecision(totalLogLike,newLogLike,logProposalRatio,acceptAllMoves);
		if (accept) m.lowCounts = 0;		
		return accept;
	}	
	
	public boolean acceptanceDecision(double oldLogLikelihood, double newLogLikelihood, double logProposalRatio,
			boolean acceptMoveIfPossible) {		
		if (Utils.DEBUG) System.out.print("logLikelihoodRatio = "+(newLogLikelihood-oldLogLikelihood));
		if (logProposalRatio > Double.NEGATIVE_INFINITY) {
			cumulativeLogProposalRatio += logProposalRatio;
		}
		else {
			return false;
		}
		if (Utils.DEBUG) System.out.println("\tlogProposalRatio = "+logProposalRatio);
		if (acceptMoveIfPossible) {						
			return (newLogLikelihood > Double.NEGATIVE_INFINITY);
		}
		return (Math.log(Utils.generator.nextDouble()) < 
				(cumulativeLogProposalRatio + heat*(newLogLikelihood - oldLogLikelihood))
				+ (cumulativeLogProposalRatio=0));
	}	
	
	/**
	 * Returns a {@link State} object that describes the current state of the
	 * MCMC. This can then be passed on to other classes such as postprocessing
	 * plugins.
	 */
	public State getState() {
		return tree.getState();
	}

	protected void beginRandomisationPeriod() { }
	protected void endRandomisationPeriod() { }			
	
	/** 
	 * Additional initialisation routines before starting MCMC, if required.
	 */
	protected void beforeMCMC() {
		System.out.println("Starting MCMC...\n");
		Utils.generator = new Well19937c(mcmcpars.seed);
	}
	/** 
	 * @return Always <code>true</true> if running non-parallel version;
	 * when running parallel version, <code>true</code> if this is the
	 * master chain.
	 */
	protected boolean isMaster() { return true; }		
	/**
	 * Starts an MCMC run. 
	 * 
	 * If <tt>AutomateParameters.shouldAutomateProposalVariances() = true</tt>
	 * then the proposal distributions will be automatically adjusted during the 
	 * burnin.
	 * 
	 * If <tt>AutomateParameters.shouldAutomateNumberOfSamples() = true</tt>
	 * or <tt>AutomateParameters.shouldAutomateStepRate() = true</tt>
	 * or <tt>AutomateParameters.shouldAutomateBurnin() = true</tt>
	 * then these parameters will be adjusted automatically, although this  
	 * approach may affect the theoretical convergence properties of the MCMC
	 * chain, so this type of automation should be regarded more as a quick way
	 * of getting some initial results without tweaking the parameters.
	 * 
	 * This function also calls the appropriate functions of the PostpocessManager
	 * <tt>postprocMan</tt> to trigger data transfer to postprocessing modules 
	 * when necessary
	 */
	public int doMCMC() {
		
		beforeMCMC();		
		
		MainFrame frame = postprocMan.mainManager.frame;
		
		//edgeWeight *= tree.vertex.length;
		
		coreModel = new CoreMcmcModule(this,modelExtMan);
		initCoreModel(tree);

		// notifies MCMC modules (including plugins) of start of MCMC sampling		
		coreModel.beforeSampling(tree);

		// Triggers a /before first sample/ of the plugins.
		postprocMan.beforeFirstSample();		

		long currentTime, start = System.currentTimeMillis();

		// calculates initial log-likelihood (includes coreModel likelihood)
		totalLogLike = modelExtMan.totalLogLike(tree);
		
		ArrayList<Double> logLikeList = new ArrayList<Double>();

		int errorCode = 0;
		
		try {
			stoppable();
			
			// Recompute progressive alignment now that everything is initialised.
			//TreeAlgo treeAlgo = new TreeAlgo();
			//treeAlgo.alignSeqsRec(tree.root);
			
			//only to use if AutomateParameters.shouldAutomate() == true
//			final int SAMPLE_RATE_WHEN_DETERMINING_THE_SPACE = 100;
//			final int BURNIN_TO_CALCULATE_THE_SPACE = 25000;
//			ArrayList<String[]> alignmentsFromSamples = new ArrayList<String[]>(); 
			int burnIn = mcmcpars.burnIn;
//			boolean stopBurnIn = false;
//			if(AutomateParameters.shouldAutomateBurnIn()){
//				burnIn = 10000000;
//			} 
//			if(AutomateParameters.shouldAutomateStepRate()){
//				burnIn += BURNIN_TO_CALCULATE_THE_SPACE;
//			}

			// Randomise the initial starting configuration 
			// by accepting all moves for a period.
			tree.root.recomputeLogLike(); // For testing
			totalLogLike = modelExtMan.totalLogLike(tree);
			if (randomisationPeriod > 0) {
				System.out.println("Randomising initial configuration for "+randomisationPeriod+" steps.");
				acceptAllMoves = true;
				beginRandomisationPeriod();		
				for (int i = 0; i < randomisationPeriod; i++) {						
					sample(0);
				}
				endRandomisationPeriod();				
				acceptAllMoves = false;
			}
			
			burnin = true;
			firstHalfBurnin = true;			
			
			tree.root.recomputeLogLike(); 
			totalLogLike = modelExtMan.totalLogLike(tree);
			
			for (int i = 0; i < burnIn; i++) {
								
				if (firstHalfBurnin && i > burnIn / 2) {
					firstHalfBurnin = false;					
					coreModel.afterFirstHalfBurnin();
					modelExtMan.afterFirstHalfBurnin();
					coreModel.incrementWeights();
					modelExtMan.incrementWeights();
					if (simulatedAnnealing) {
						heat = 1;						
					}
					if (isParallel) {
						heat = heatWhenHot;
					}
				}
				else {
					if (simulatedAnnealing) {
						heat = Math.log(i) / Math.log(burnIn / 2); 
					}
				}
				// Perform an MCMC move
				sample(0);
				

				// Triggers a /new step/ and a /new peek/ (if appropriate) of
				// the plugins.
				//if (isMaster()) {
					// TODO do above inside sample() and add more info
					mcmcStep.newLogLike = modelExtMan.totalLogLike(tree);
					mcmcStep.burnIn = burnin;
					postprocMan.newStep(mcmcStep);
					if (i % mcmcpars.sampRate == 0) {
						postprocMan.newPeek();
					}
				//}
				if (i>0 && mcmcpars.doReportDuringBurnin && (i % mcmcpars.sampRate == 0)) {
					report(i, mcmcpars.cycles / mcmcpars.sampRate);
				}
//				if(AutomateParameters.shouldAutomateBurnIn() && i % 50 == 0){
//					// every 50 steps, add the current loglikelihood to a list
//					// and check if we find a major decline in that list 
//					logLikeList.add(getState().logLike);
//					if(!stopBurnIn){
//						stopBurnIn = AutomateParameters.shouldStopBurnIn(logLikeList);
//						if(AutomateParameters.shouldAutomateStepRate() && stopBurnIn){
//							burnIn = i + BURNIN_TO_CALCULATE_THE_SPACE;
//						}else if (stopBurnIn){
//							burnIn = i;
//						}
//					}
//				}
				currentTime = System.currentTimeMillis();
//				int realBurnIn = burnIn - BURNIN_TO_CALCULATE_THE_SPACE;
				if (frame != null) {
					String text = "";
//					if((i > realBurnIn ) && AutomateParameters.shouldAutomateStepRate()){
//						text = "Burn-in to aid automation of MCMC parameters: " + (i-realBurnIn + 1) ;
//					}else{
						text = "Burn-in: " + (i + 1);
//					}
					frame.statusText.setText(text);
				} else if (i % 1000 == 999) {
					if (Utils.DEBUG) {
						System.err.println("Burn in: " + (i + 1));
					}
					else {
						System.out.println("Burn in: " + (i + 1));
					}
				}
//				if( AutomateParameters.shouldAutomateStepRate() && (i >= realBurnIn) && i % SAMPLE_RATE_WHEN_DETERMINING_THE_SPACE == 0)   {
//					String[] align = getState().getLeafAlign();
//					alignmentsFromSamples.add(align);
//				}	
				if (AutomateParameters.shouldAutomateProposalVariances() && i % mcmcpars.sampRate == 0) {
					coreModel.modifyProposalWidths();
					modelExtMan.modifyProposalWidths();
				}
			}
			
			//both real burn-in and the one to determine the sampling rate have now been completed.
			burnin = false;
			coreModel.afterBurnin();
			modelExtMan.afterBurnin();
			coreModel.zeroAllMoveCounts();
			modelExtMan.zeroAllMoveCounts();
			
			//Utils.DEBUG = true;
			
			int period;
//			if(AutomateParameters.shouldAutomateNumberOfSamples()){
//				period = 1000000;
//			}else{
				period = mcmcpars.cycles / mcmcpars.sampRate;
//			}

			int sampRate;
//			if(AutomateParameters.shouldAutomateStepRate()){
//				if(frame != null)
//				{
//					frame.statusText.setText("Calculating the sample rate");
//				}
//				else
//				{
//					System.out.println("Calculating the sample rate");
//				}
//				ArrayList<Double> theSpace = Distance.spaceAMA(alignmentsFromSamples);
//				sampRate = AutomateParameters.getSampleRateOfTheSpace(theSpace,SAMPLE_RATE_WHEN_DETERMINING_THE_SPACE);
//
//			}else{
				sampRate = mcmcpars.sampRate;
//			}
		
			
//			AlignmentData alignment = new AlignmentData(getState().getLeafAlign());
//			ArrayList<AlignmentData> allAlignments = new ArrayList<AlignmentData>();
//			ArrayList<Double> distances = new ArrayList<Double>();

			boolean shouldStop = false;
//			double currScore = 0;
			for (int i = 0; i < period && !shouldStop; i++) {
				if (i > 0 && (i % LOG_INTERVAL == 0)) {
					postprocMan.flushAll();
					if (coreModel.printExtraInfo) printMcmcInfo();
				}
				for (int j = 0; j < sampRate; j++) {
					
					// Perform an MCMC move
					sample(0);				
					
					// Proposes a swap.
					if (isParallel && ((i*sampRate + j) % mcmcpars.swapRate == 0)) {
						doSwap();						
					}

					// Triggers a /new step/ and a /new peek/ (if appropriate)
					// of the plugins.
					//if (isMaster()) {
						mcmcStep.newLogLike = totalLogLike;
						mcmcStep.burnIn = burnin;
						postprocMan.newStep(mcmcStep);
						if (burnIn + i * period + j % mcmcpars.sampRate == 0) {
							postprocMan.newPeek();
						}
					//}

					currentTime = System.currentTimeMillis();
					if (frame != null) {
						String text = "Samples taken: " + Integer.toString(i);
//						//remainingTime((currentTime - start)
//						//		* ((period - i - 1) * sampRate
//						//				+ sampRate - j - 1)
//						//				/ (burnIn + i * sampRate + j + 1))
//
						text += "   The sampling rate: " + sampRate;
//						if(AutomateParameters.shouldAutomateNumberOfSamples()){
//							text +=  ",  Similarity(alignment n-1, alignment n): " + df.format(currScore) + " < " + df.format(AutomateParameters.PERCENT_CONST);
//						}
						frame.statusText.setText(text );
					}
				}
				if (frame == null && !isParallel) {
					if (Utils.DEBUG) {
						System.err.println("Sample: " + (i + 1));	
					}
					else {
						System.out.println("Sample: " + (i + 1));
					}
				}
//				if(AutomateParameters.shouldAutomateNumberOfSamples()){
//					alignment = new AlignmentData(getState().getLeafAlign());
//					allAlignments.add(alignment);
//					if (allAlignments.size() >1){
//						FuzzyAlignment Fa = FuzzyAlignment.getFuzzyAlignmentAndProject(allAlignments.subList(0, allAlignments.size()-1), 0);
//						FuzzyAlignment Fb = FuzzyAlignment.getFuzzyAlignmentAndProject(allAlignments, 0);
//						currScore = FuzzyAlignment.AMA(Fa, Fb);
//						System.out.println(currScore);
//						distances.add(currScore);
//						if (allAlignments.size() >5){
//							shouldStop = AutomateParameters.shouldStopSampling(distances);
//						}
//
//					}
//				}
				// Report the results of the sample.
				report(i, period);
			}
		} catch (StoppedException ex) {
			errorCode = 1;
			// stopped: report and save state
		}

		// Triggers a /after first sample/ of the plugins.
		if (isMaster()) {
			printMcmcInfo();
			postprocMan.afterLastSample();
		}
		
		// notifies model extension plugins of the end of sampling
		modelExtMan.afterSampling();		
		System.out.println(coreModel.getSummaryInfo());
		coreModel.afterSampling();
		
		if (frame != null) {
			frame.statusText.setText(MainFrame.IDLE_STATUS_MESSAGE);
		}
		return errorCode;
	}

	public String getInfoString() {
		return coreModel.getSummaryInfo();
	}
	private void printMcmcInfo() {
		String info = "\n"+Utils.repeatedString("#",64)+"\n";
		info += String.format("%-24s","Move name")+
		String.format("%8s","t")+
		String.format("%8s","nMoves")+
		String.format("%8s","t/move")+
		String.format("%8s", "acc")+
		String.format("%8s\n", "propVar");
		info += coreModel.getMcmcInfo();
		info += modelExtMan.getMcmcInfo();
		info += Utils.repeatedString("#",64)+"\n";
		System.out.println(info);
	}
	/** 
	 * Triggers <tt>postProcMan</tt> to print out a report of the current
	 * state of the chain.
	 * 
	 * @param no
	 * @param total
	 */
	private void report(int no, int total) {
		report(no,total,true);
	}
	private void reportDuringBurnin(int no, int total) {
		report(no,total,false);
	}
	protected void report(int no, int total, boolean useSample) {

		if (useSample) postprocSample(no,total);		
		// Log the accept ratios/params to the (.log) file. TODO: move to a plugin.
		try {			
			postprocMan.logFile.write(coreModel.getSummaryInfo() + "\n");
			coreModel.printParameters();			
		} catch (IOException e) {
			if (postprocMan.mainManager.frame != null) {
				ErrorMessage.showPane(postprocMan.mainManager.frame, e, true);
			} else {
				e.printStackTrace(System.out);
			}
		}
	}
	
	protected void postprocSample(int no, int total) {
		postprocMan.newSample(coreModel,getState(), no, total);
	}
	
	// This function is used (and defined) only by the parallel version		
	protected void doSwap() { }	
}
