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
public class Mcmc extends Stoppable {

	// Parallelization

	/** Is this a parallel chain? By-default false. */
	private boolean isParallel = false;
	
	private boolean simulatedAnnealing = false;

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

	/** Total log-likelihood of the current state, cached for speed */
	double totalLogLike;
	
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
	
	// TODO Move the parameters below into MCMCPars. Would be nice to have
	// a set of sliding bars in a menu in the GUI that go from 0 to 100, 
	// to select the relative frequency of the different moves etc.
	
	// Which indel parameter move scheme(s) to use
	private boolean lambdaMuMove = false;
	private boolean lambdaMove = true;
	private boolean lambdaPhiMove = false;
	private boolean rhoThetaMove = true;
	
	// Weights for coreModel McmcMoves
	private int rWeight = 8;
	private int lambdaWeight = 4;
	private int muWeight = 6;
	private int lambdaMuWeight = 6;
	private int phiWeight = 4;
	private int rhoWeight = 6;
	private int thetaWeight = 6;
	
	private int substWeight = 10;
	private int edgeWeight = 1; // per edge
	private int allEdgeWeight = 6;
	private int edgeWeightIncrement = 0; // Added after half of burnin
	
	//private int alignWeight = 25;
	private double[] Pvals = {0.9,0.99,0.999,0.9999,0.99999};
	private double[] rateMult = {3.0,2.5,2.0,1.5,1.0};
	//private double[] rateMult = {1,1,1,1,1};

	private int[] alignWeights = {0,0,0,0,12}; // ORIGINAL
	private int[] alignWeights2 = {0,0,0,0,12}; // ORIGINAL
	//private int[] alignWeights2 = {0,0,0,0,0};
	//private int[] alignWeights = {0,0,0,0,25};
	private int[] alignWeightIncrements = {5,2,2,2,-4}; // ORIGINAL
	private int[] alignWeightIncrements2 = {2,2,2,2,-8}; // ORIGINAL
	//private int[] alignWeightIncrements2 = {5,5,5,5,5};
	//private int[] alignWeightIncrements = {5,1,1,1,-10};
	//private int[] alignWeightIncrements = {11,1,1,1,-10};
	//private int[] alignWeightIncrements = {0,0,0,0,0};
	//private int[] alignWeightIncrements2 = {0,0,0,0,0};
	
	private int silentIndelWeight = 10;
	private int silentIndelWeightIncrement = -8; // Added after half of burnin
	private int topologyWeight = 8;
	private int localTopologyWeight = 8;
	private int topologyWeightIncrement = 12; // Added after half of burnin
	private int localTopologyWeightIncrement = 12; // Added after half of burnin
	private int topologyWeightDuringRandomisationPeriod = 100; // To use while we're randomising initial config
	
	LOCALTopologyMove localTopologyMove;
	TopologyMove topologyMove;
	
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
		this.tree.heat = 1.0d;
		randomisationPeriod = mcmcpars.randomisationPeriod;
	}
	public Mcmc(Tree tree, MCMCPars pars, PostprocessManager ppm, ModelExtManager modelExtMan,
			int noOfProcesses, int rank, double heat) {
		this(tree, pars, ppm, modelExtMan);
		this.noOfProcesses = noOfProcesses;
		this.rank = rank;
		this.tree.heat = heat;

		// Is parallel!
		isParallel = true;
	}
	
	private static final DecimalFormat df = new DecimalFormat("0.0000");

	/**
	 * Initialises the coreModel, adding the various MCMC moves
	 * with the specified priors and proposal distributions where 
	 * appropriate. Currently the coreModel cannot be modified from the 
	 * command line nor from within the GUI. 
	 */
	private void initCoreModel() {
		
		coreModel.printExtraInfo = Utils.VERBOSE;
		
		double[] lambdaMuPriorParams = {1,1};
		
		IndelMove rMove = new RMove(coreModel,new BetaPrior(1,1),new LogisticProposal(),"R");
		rMove.proposalWidthControlVariable = 1.0;
		coreModel.addMcmcMove(rMove,rWeight);				
		
		if (lambdaMuMove) {
			IndelMove lambdaMove = new LambdaMove(coreModel,new GammaPrior(lambdaMuPriorParams[0],lambdaMuPriorParams[1]),new GaussianProposal(),"Lambda");
			lambdaMove.proposalWidthControlVariable = 0.01;
			coreModel.addMcmcMove(lambdaMove,lambdaWeight);
			IndelMove muMove = new MuMove(coreModel,new GammaPrior(lambdaMuPriorParams[0],lambdaMuPriorParams[0]),new GaussianProposal(),"Mu");
			muMove.proposalWidthControlVariable = 0.01;
			coreModel.addMcmcMove(muMove,muWeight);
			ArrayList<McmcMove> lambdaMu = new ArrayList<McmcMove>();
			lambdaMu.add(lambdaMove);
			lambdaMu.add(muMove);
			coreModel.addMcmcMove(new McmcCombinationMove(lambdaMu),lambdaMuWeight);
		}
		if (lambdaMove || lambdaPhiMove) {
			IndelMove lambdaMove = new LambdaMove(coreModel,new GammaPrior(lambdaMuPriorParams[0],lambdaMuPriorParams[1]),new GaussianProposal(),"Lambda");
			coreModel.addMcmcMove(lambdaMove,lambdaWeight);
			lambdaMove.proposalWidthControlVariable = 0.01;
		}
		if (lambdaPhiMove) {
			// phi = lambda/mu 
			// Must be in the range (0,1)
			// The prior on phi below [Beta(1,1)] is not the prior implied by the priors on lambda and mu
			// (which would be of F-distribution type), so effectively it means we have two different 
			// priors on lambda and mu if we use this move, which is undesirable. 
			// Hence, this move should probably be avoided. 
			IndelMove phiMove = new PhiMove(coreModel,new BetaPrior(1,1),new LogisticProposal(),"Phi"); 
			phiMove.proposalWidthControlVariable = 0.5;
			coreModel.addMcmcMove(phiMove,phiWeight);
		}
		if (rhoThetaMove) {
			// rho = lambda + mu
			// Since lambda and mu are gamma distributed in the prior, this ratio
			// also follows a gamma distribution in the prior.
			IndelMove rhoMove = new RhoMove(coreModel,new GammaPrior(2*lambdaMuPriorParams[0],lambdaMuPriorParams[1]),new GaussianProposal(),"Rho");
			coreModel.addMcmcMove(rhoMove,rhoWeight);
			rhoMove.proposalWidthControlVariable = 0.02;

			// theta = lambda / (lambda + mu)
			// Must be in the range (0,0.5)
			// Since lambda and mu are gamma distributed in the prior, this ratio
			// follows a beta distribution in the prior.
			IndelMove thetaMove = new ThetaMove(coreModel,new BetaPrior(lambdaMuPriorParams[0],lambdaMuPriorParams[0]),new GaussianProposal(),"Theta");
			thetaMove.setMaxValue(0.5);
			thetaMove.proposalWidthControlVariable = 0.02;
			coreModel.addMcmcMove(thetaMove,thetaWeight);
		}
		if (!lambdaMuMove && !lambdaPhiMove && !rhoThetaMove) {
			throw new IllegalArgumentException("Invalid proposal scheme selected for indel parameters.");
		}
		
		SubstMove substMove = new SubstMove(coreModel,"Subst");
		coreModel.addMcmcMove(substMove, substWeight);

		if(!mcmcpars.fixAlign) {
			for (int i=0; i<Pvals.length; i++) {
				AlignmentMove alignMove = new AlignmentMove(coreModel,Pvals[i],"Alignment_"+i+"_"+Pvals[i]);
				if (i==(Pvals.length-1)) { // Keep one of them with longer windows
					alignMove.autoTunable = false;
					//alignMove.useModelExtInProposal(); // ORIGINAL
				}				
				alignMove.useModelExtInProposal(); // NOT ORIGINAL
				coreModel.addMcmcMove(alignMove, alignWeights[i],alignWeightIncrements[i]);
			}
			for (int i=0; i<rateMult.length; i++) {
				AlignmentMove alignMove = new AlignmentMove(coreModel,Pvals[4],"Alignment_"+(i+Pvals.length)+"_"+Pvals[4]+"_"+rateMult[i]);
				alignMove.setProposalParamMultiplier(rateMult[i]);
				if (i==(rateMult.length-1)) { // Keep one of them with longer windows
					alignMove.autoTunable = false;
					//alignMove.useModelExtInProposal(); // ORIGINAL
				}
				alignMove.useModelExtInProposal(); // NOT ORIGINAL
				coreModel.addMcmcMove(alignMove, alignWeights2[i],alignWeightIncrements2[i]);
			}
			
			SilentIndelMove silentIndelMove = new SilentIndelMove(coreModel,"SilentIndel");
			coreModel.addMcmcMove(silentIndelMove, silentIndelWeight,silentIndelWeightIncrement);
		}

		GammaPrior edgePrior = new GammaPrior(1,0.01);
		//GammaPrior edgePrior = new GammaPrior(1,2);
		double uniformProposalWidthControlVariable = 0.25;
		double multiplicativeProposalWidthControlVariable = 0.5;
		
		topologyMove = null;
		double fastSwapProb = 0.05;
		if (mcmcpars.fixAlign) fastSwapProb = 0.0;
		
		if(!mcmcpars.fixTopology && !mcmcpars.fixEdge && tree.vertex.length > 3) {
			topologyMove = new TopologyMove(coreModel,edgePrior,
					//0.5*multiplicativeProposalWidthControlVariable,"Topology"); // works ok with glob_25
					0.75*uniformProposalWidthControlVariable,fastSwapProb,"Topology"); // experimental
			coreModel.addMcmcMove(topologyMove, topologyWeight,topologyWeightIncrement);
						
//			LOCALTopologyMove localTopologyMove = new LOCALTopologyMove(coreModel,edgePrior,
//					0.5*multiplicativeProposalWidthControlVariable,"LOCALTopology");
			localTopologyMove = new LOCALTopologyMove(coreModel,edgePrior,
					1*uniformProposalWidthControlVariable,fastSwapProb,"LOCALTopology");
			coreModel.addMcmcMove(localTopologyMove, localTopologyWeight,localTopologyWeightIncrement);
		}
		if(!mcmcpars.fixEdge) {
			//HyperbolicPrior edgePrior = new HyperbolicPrior();
			for (int i=0; i<tree.vertex.length-1; i++) {
				EdgeMove edgeMove = new EdgeMove(coreModel,i,
						edgePrior,
						//new GaussianProposal(),
						new UniformProposal(),
						"Edge"+i);
				edgeMove.proposalWidthControlVariable = uniformProposalWidthControlVariable;
				// Default minimum edge length is 0.01
				coreModel.addMcmcMove(edgeMove, edgeWeight,edgeWeightIncrement);
			}		
			AllEdgeMove allEdgeMove = new AllEdgeMove(coreModel,edgePrior,
					new MultiplicativeProposal(),"AllEdge");
			allEdgeMove.proposalWidthControlVariable = multiplicativeProposalWidthControlVariable;			
			coreModel.addMcmcMove(allEdgeMove, allEdgeWeight);
						
		}
	}
	
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
				(cumulativeLogProposalRatio + tree.heat*(newLogLikelihood - oldLogLikelihood))
				+ (cumulativeLogProposalRatio=0));
	}	
	/**
	 * Returns a string representation describing the acceptance ratios of the current MCMC run.
	 * @return a string describing the acceptance ratios.
	 */
	public String getInfoString() {
		String info = "Acceptance rates: ";
		for (McmcMove m : coreModel.getMcmcMoves()) {
			info += m.name+": "+String.format(Locale.US, "%f ", m.acceptanceRate());
		}
		return info;
	}
	
	/**
	 * Returns a {@link State} object that describes the current state of the
	 * MCMC. This can then be passed on to other classes such as postprocessing
	 * plugins.
	 */
	public State getState() {
		return tree.getState();
	}

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
		Utils.generator = new Well19937c(mcmcpars.seed + rank);

		if(tree.substitutionModel.params == null || tree.substitutionModel.params.length == 0) {
			substWeight = 0;
		}
		//edgeWeight *= tree.vertex.length;
		
		coreModel = new CoreMcmcModule(this,modelExtMan);
		initCoreModel(); 
		// initCoreModel() uses the weights, so they need to be defined
		// and updated before this point.

		// notifies model extension plugins of start of MCMC sampling
		modelExtMan.beforeSampling(tree);

		// Triggers a /before first sample/ of the plugins.
		if ((isParallel && MPIUtils.isMaster(rank)) || !isParallel) {
			postprocMan.beforeFirstSample();
		}

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
				coreModel.setWeight("Topology",topologyWeightDuringRandomisationPeriod);
			}
			for (int i = 0; i < randomisationPeriod; i++) {						
				sample(0);
			}
			coreModel.setWeight("Topology",topologyWeight);
			coreModel.zeroAllMoveCounts();
			modelExtMan.zeroAllMoveCounts();
			acceptAllMoves = false;
			
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
						tree.heat = 1;
					}
				}
				else {
					if (simulatedAnnealing) {
						tree.heat = Math.log(i) / Math.log(burnIn / 2); 
					}
				}
				// Perform an MCMC move
				sample(0);
				

				// Triggers a /new step/ and a /new peek/ (if appropriate) of
				// the plugins.
				if ((isParallel && MPIUtils.isMaster(rank)) || !isParallel) {
					// TODO do above inside sample() and add more info
					mcmcStep.newLogLike = modelExtMan.totalLogLike(tree);
					mcmcStep.burnIn = burnin;
					postprocMan.newStep(mcmcStep);
					if (i % mcmcpars.sampRate == 0) {
						postprocMan.newPeek();
					}
				}
				if (mcmcpars.doReportDuringBurnin && (i % mcmcpars.sampRate == 0)) {
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
					System.out.println("Burn in: " + (i + 1));
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


			int swapNo = 0; // TODO: delete?
			swapCounter = mcmcpars.swapRate;
			
			
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
						mcmcStep.newLogLike = totalLogLike;
						mcmcStep.burnIn = burnin;
						postprocMan.newStep(mcmcStep);
						if (burnIn + i * period + j % mcmcpars.sampRate == 0) {
							postprocMan.newPeek();
						}
					}

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
					System.out.println("Sample: " + (i + 1));
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

		//if(Utils.DEBUG) {
			printMcmcInfo();
		//}

		// Triggers a /after first sample/ of the plugins.
		if ((isParallel && MPIUtils.isMaster(rank)) || !isParallel) {
			postprocMan.afterLastSample();
		}
		
		// notifies model extension plugins of the end of sampling
		modelExtMan.afterSampling();
		
		System.out.println(getInfoString());
		System.out.println("Acceptance rate for LOCAL moves resulting in topology change = "+
				String.format(Locale.US, "%.4f",(double)localTopologyMove.nTopologyChanges/(double)localTopologyMove.proposalCount));
		System.out.println("Acceptance rate for Topology moves resulting in topology change = "+
				String.format(Locale.US, "%.4f",(double)topologyMove.topologyChangeAcceptanceCount/(double)topologyMove.proposalCount));
		
		if (frame != null) {
			frame.statusText.setText(MainFrame.IDLE_STATUS_MESSAGE);
		}
		return errorCode;
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
	private void report(int no, int total, boolean useSample) {

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
				if (useSample) postprocMan.newSample(coreModel,getState(), no, total);
			} else if (isColdChain() && !MPIUtils.isMaster(rank)) {
				// Send state.
				State state = getState();
				MPIStateSend(state);
			} else if (!isColdChain() && MPIUtils.isMaster(rank)) {
				// Receive state.
				State state = MPIStateReceieve(coldChainLocation);
				if (useSample) postprocMan.newSample(coreModel,state, no, total);
			}

		} else {
			if (useSample) postprocMan.newSample(coreModel,getState(), no, total);
		}

		// Log the accept ratios/params to the (.log) file. TODO: move to a plugin.
		try {
			if ((isParallel && MPIUtils.isMaster(rank)) || !isParallel) {
				postprocMan.logFile.write(getInfoString() + "\n");
				postprocMan.logFile.write("Report\tLogLikelihood\t"
						+ (modelExtMan.totalLogLike(tree))
						+ "\tR\t" + tree.hmm2.params[0] + "\tLambda\t"
						+ tree.hmm2.params[1] + "\tMu\t" + tree.hmm2.params[2]
								+ "\t" + tree.substitutionModel.print() + "\n");
				if (isParallel) {
					postprocMan.logFile.write("Cold chain location: " + coldChainLocation + "\n");
				}

			}
		} catch (IOException e) {
			if (postprocMan.mainManager.frame != null) {
				ErrorMessage.showPane(postprocMan.mainManager.frame, e, true);
			} else {
				e.printStackTrace(System.out);
			}
		}
	}
	
	// The functions below are used only by the parallel version
	
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
			myStateInfo[0] = totalLogLike;
			myStateInfo[1] = modelExtMan.totalLogPrior(tree);
			//myStateInfo[1] = coreModel.totalLogPrior(tree) + modelExtMan.totalLogPrior(tree);
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
}
