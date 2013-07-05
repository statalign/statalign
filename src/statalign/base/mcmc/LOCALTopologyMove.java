package statalign.base.mcmc;

import statalign.base.Tree;
import statalign.base.Utils;
import statalign.base.Vertex;
import statalign.mcmc.McmcModule;
import statalign.mcmc.McmcMove;
import statalign.mcmc.MultiplicativeProposal;
import statalign.mcmc.PriorDistribution;

public class LOCALTopologyMove extends McmcMove {

	Tree tree = null;

	Vertex nephew, parent, uncle;
	
	double w_ij, w_ai, w_aj, w_ac;
			
	double minEdgeLength;
	
	boolean topologyChange;
	boolean invalidProposal;
	/** Probability of selecting the fast nephew-uncle switch */
	double fastSwapProb;
	/** <tt>true</tt> if we used <tt>fastSwapWithUncle</tt>*/
	boolean didFastSwap;
	public int nTopologyChanges = 0;
	
	PriorDistribution<Double> edgePrior;
	
	public LOCALTopologyMove (McmcModule m, PriorDistribution<Double> pr, double propVar, double _fastSwapProb, String n) {
		owner = m;
		edgePrior = pr;
		proposalWidthControlVariable = propVar;
		name = n;
		autoTune = false; 
		// autoTune = true by default
		minEdgeLength = Utils.MIN_EDGE_LENGTH;
		fastSwapProb = _fastSwapProb;
		//fastSwapProb = 0.0;
	}
	
	/**
	 * Sets the three edges affected by the LOCAL move (nephew, parent and uncle) 
	 * to the specified values. If any of the values are less than <tt>minEdgeLength</tt>
	 * then the change is not made, <tt>invalidProposal</tt> is set to <tt>true</tt>, 
	 * and negative infinity is returned.
	 * @param a The edge length for the nephew
	 * @param b The edge length for the parent
	 * @param c The edge length for the uncle
	 * @return The log ratio of new versus old prior densities.
	 */
	public double setEdges(double a, double b, double c) {
		
		if (invalidProposal || a < minEdgeLength || b < minEdgeLength || c < minEdgeLength) {
			invalidProposal = true;
			return Double.NEGATIVE_INFINITY;
		}
		double logPrior = -edgePrior.logDensity(nephew.edgeLength);
		nephew.setEdgeLength(a);
		logPrior += edgePrior.logDensity(nephew.edgeLength);
		nephew.edgeChangeUpdate();
		nephew.calcAllUp(); //
		
		logPrior -= edgePrior.logDensity(parent.edgeLength);
		parent.setEdgeLength(b);
		logPrior += edgePrior.logDensity(parent.edgeLength);
		parent.edgeChangeUpdate();
		parent.calcAllUp();// 
		
		logPrior -= edgePrior.logDensity(uncle.edgeLength);
		uncle.setEdgeLength(c);
		logPrior += edgePrior.logDensity(uncle.edgeLength);
		uncle.edgeChangeUpdate();
		uncle.calcAllUp();//
//		tree.root.calcFelsRecursively();
//		tree.root.calcIndelLikeRecursively();
		return logPrior;
	}

	@Override
	public void copyState(Object externalState) {
		if (externalState instanceof Tree) {
			if (tree == null) {
				tree = (Tree) externalState;
			}
		}
		else {
			throw new IllegalArgumentException("TopologyMove.copyState must take an argument of type Tree.");
		}
		int vnum = tree.vertex.length;
		if (vnum <= 3) {
			return;
		}
		
		int rnd = Utils.generator.nextInt(vnum);
		Vertex vertex_a, vertex_j, vertex_i = tree.vertex[rnd];
		for (;;) {
			rnd = Utils.generator.nextInt(vnum);
			vertex_i = tree.vertex[rnd];
			if (vertex_i.right == null && vertex_i.left == null) {
				continue;
			}
			else {
				if (Utils.generator.nextDouble() < 0.5) {
					vertex_j = vertex_i.left;
					vertex_a = vertex_i.right;
				}	
				else {
					vertex_j = vertex_i.right;
					vertex_a = vertex_i.left;
				}
				if (vertex_j.right == null && vertex_j.left == null) {
					continue;
				}
				break;
			}
		}
		
		w_ij = vertex_j.edgeLength;
				
		Vertex vertex_c = (Utils.generator.nextDouble() < 0.5) ? vertex_j.left : vertex_j.right; 
		
		w_ai = vertex_a.edgeLength;
		w_aj = w_ai + w_ij;
		w_ac = w_aj + vertex_c.edgeLength;
		
		nephew = vertex_c;
		parent = vertex_j;
		uncle = vertex_a;
		// grandpa = vertex_i;
		
		((CoreMcmcModule) owner).getModelExtMan().beforeTreeChange(tree, vertex_c);		
		// Should also do a beforeAlignChange here, but not obvious what to pass
		// as the selectedRoot argument.
		// The rest of the state copying is handled inside the Vertex
	}
	@Override
	public double proposal(Object externalState) {
		
		double u1 = Utils.generator.nextDouble();
		double u2 = Utils.generator.nextDouble();
		double r = Math.exp(proposalWidthControlVariable*(u1-0.5));
		double w_ac_new = r * w_ac; 
		double w_ai_new = r * w_ai;
		double w_aj_new = u2 * w_ac_new;
		double logProposalRatio = 3 * Math.log(r);
		if (Utils.DEBUG) System.out.println("Before LOCAL: "+tree.printedTreeWithNumbers());
		invalidProposal = false;
		if (w_aj_new < w_ai_new) { // Then we have a topology switch
	        if (Utils.DEBUG) tree.root.printToScreenAlignment(0,0,true);
			topologyChange = true;
			double w_ij_new = w_ai_new - w_aj_new;
			double w_ic_new = w_ac_new - w_ai_new;
			// Add in priorDensity to logProposalRatio
			logProposalRatio += setEdges(w_ic_new, w_ij_new, w_aj_new); 
		}
		else {
			//System.out.println("Only changing edges.");
			topologyChange = false;
			double w_jc_new = w_ac_new - w_aj_new;
			double w_ij_new = w_aj_new - w_ai_new;
			logProposalRatio += setEdges(w_jc_new, w_ij_new, w_ai_new);
		}
		if (topologyChange && !invalidProposal) {
			// Do the topology switch:
	    	if (Utils.USE_MODEXT_EM) fastSwapProb = 0;
			if (Utils.generator.nextDouble() < fastSwapProb) {
				logProposalRatio += nephew.fastSwapWithUncle();
				didFastSwap = true;
		        if (Utils.DEBUG) tree.root.printToScreenAlignment(0,0,true);
			}
			else {
				//logProposalRatio += nephew.swapWithUncleAlignToParent();
				//logProposalRatio += nephew.nephewUncleSwapFixedColumns();
				logProposalRatio += nephew.nephewUncleSwapFixedColumns3();
				didFastSwap = false;
			}
		}
		if (Utils.DEBUG) System.out.println("After  LOCAL: "+tree.printedTreeWithNumbers());
		return logProposalRatio;
	}
	@Override
	public double logPriorDensity(Object externalState) {
		return 0.0;
	}
	@Override
	public void updateLikelihood(Object externalState) {
		owner.setLogLike(((CoreMcmcModule) owner).getModelExtMan().logLikeTreeChange(tree, nephew));
	}
	@Override
	public void restoreState(Object externalState) {
		if (topologyChange && !invalidProposal) {
			if (didFastSwap) {
				uncle.fastSwapBackUncle();
			}
			else {
				uncle.restoreFiveWay();
		        //uncle.swapBackUncleAlignToParent();
			}
		}
		setEdges(w_ac-w_aj, w_ij, w_ai);
		if (didFastSwap) tree.root.recomputeLogLike();
	}
	
	@Override
	public void afterMove(Object externalState) {
		if (lastMoveAccepted && topologyChange && nephew.parent != tree.root) {
			System.out.println("Topology change (LOCAL).");
			nTopologyChanges++;
		}
		((CoreMcmcModule) owner).getModelExtMan().afterTreeChange(tree,lastMoveAccepted ? uncle : nephew,lastMoveAccepted);
		// Should also do an afterAlignChange here, but not obvious what to pass
		// as the selectedRoot argument.
		
		if (Utils.DEBUG) {
			tree.checkPointers();
		}
		if (Utils.USE_MODEXT_EM && lastMoveAccepted) {
			uncle.parent.updateAligned();
			nephew.updateAlignedParent();
			if (Utils.DEBUG) tree.root.updateAlignedRecursivelyWithCheck();
		}
	}
	
	 
	
}
