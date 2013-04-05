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
		
	double edgeProposalWidthControlVariable;
	
	boolean topologyChange;
	
	PriorDistribution<Double> edgePrior;
	
	public LOCALTopologyMove (McmcModule m, PriorDistribution<Double> pr, double propVar, String n) {
		owner = m;
		edgePrior = pr;
		edgeProposalWidthControlVariable = propVar;
		name = n;
		autoTune = false; 
		// Already true by default, but let's write it here to be explicit
	}
	
	public double setEdges(double a, double b, double c) {
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
		//System.out.println("Before LOCAL: "+tree.printedTree());
		if (w_aj_new < w_ai_new) { // Then we have a topology switch
			topologyChange = true;
			double w_ij_new = w_ai_new - w_aj_new;
			double w_ic_new = w_ac_new - w_ai_new;
			// Add in priorDensity to logProposalRatio
			logProposalRatio += setEdges(w_ic_new, w_ij_new, w_aj_new);

			// Do the topology switch:
			logProposalRatio += nephew.fastSwapWithUncle();
			// Below is another version, slow and slightly better mixing
			// logProposalRatio += nephew.swapWithUncleAlignToParent();
		}
		else {
			//System.out.println("Only changing edges.");
			topologyChange = false;
			double w_jc_new = w_ac_new - w_aj_new;
			double w_ij_new = w_aj_new - w_ai_new;
			logProposalRatio += setEdges(w_jc_new, w_ij_new, w_ai_new);
		}
		//System.out.println("After LOCAL ("+logProposalRatio+"): "+tree.printedTree());
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
		if (topologyChange) {
			uncle.fastSwapBackUncle();
		}
		// If using the alternative move:
        // uncle.swapBackUncleAlignToParent();
		
		setEdges(w_ac-w_aj, w_ij, w_ai);
	}
	
	@Override
	public void afterMove(Object externalState) {
		((CoreMcmcModule) owner).getModelExtMan().afterTreeChange(tree,lastMoveAccepted ? uncle : nephew,lastMoveAccepted);
		// Should also do an afterAlignChange here, but not obvious what to pass
		// as the selectedRoot argument.
		
		if (Utils.DEBUG) {
			tree.checkPointers();
		}
	}
	
	 
	
}
