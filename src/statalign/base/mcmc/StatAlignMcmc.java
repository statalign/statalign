package statalign.base.mcmc;

import java.util.ArrayList;

import statalign.base.MCMCPars;
import statalign.base.Mcmc;
import statalign.base.Tree;
import statalign.base.Utils;
import statalign.mcmc.BetaPrior;
import statalign.mcmc.GammaPrior;
import statalign.mcmc.GaussianProposal;
import statalign.mcmc.LogisticProposal;
import statalign.mcmc.McmcCombinationMove;
import statalign.mcmc.McmcMove;
import statalign.mcmc.MultiplicativeProposal;
import statalign.mcmc.UniformProposal;
import statalign.model.ext.ModelExtManager;
import statalign.postprocess.PostprocessManager;

/**
 * 
 * Contains the specifics of the MCMC scheme for StatAlign. Currently this 
 * includes various hard-coded parameters, and choice of move types.
 * 
 * @author herman
 *
 */
public class StatAlignMcmc extends Mcmc {
	
	public StatAlignMcmc(Tree tree, MCMCPars pars, PostprocessManager ppm, ModelExtManager modelExtMan) {
		super(tree,pars,ppm,modelExtMan);		
	}
	public StatAlignMcmc(Tree tree, MCMCPars pars, PostprocessManager ppm, ModelExtManager modelExtMan,
			int noOfProcesses, int rank, double heat) {
		super(tree,pars,ppm,modelExtMan,noOfProcesses,rank,heat);
	}
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
	
	protected void beginRandomisationPeriod() {
		coreModel.setWeight("Topology",topologyWeightDuringRandomisationPeriod);
	}
	protected void endRandomisationPeriod() {
		coreModel.setWeight("Topology",topologyWeight);
		coreModel.zeroAllMoveCounts();
		modelExtMan.zeroAllMoveCounts();
	}
	/**
	 * Initialises the coreModel, adding the various MCMC moves
	 * with the specified priors and proposal distributions where 
	 * appropriate. Currently the coreModel cannot be modified from the 
	 * command line nor from within the GUI. 
	 */
	@Override
	protected void initCoreModel(Tree tree) {
				
		coreModel.printExtraInfo = Utils.VERBOSE;
		
		if(tree.substitutionModel.params == null || tree.substitutionModel.params.length == 0) {
			substWeight = 0;
		}
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
}
