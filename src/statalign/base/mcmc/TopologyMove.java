package statalign.base.mcmc;

import statalign.base.Tree;
import statalign.base.Utils;
import statalign.base.Vertex;
import statalign.mcmc.McmcModule;
import statalign.mcmc.McmcMove;
import statalign.mcmc.MultiplicativeProposal;
import statalign.mcmc.PriorDistribution;
import statalign.mcmc.UniformProposal;

public class TopologyMove extends McmcMove {

	Tree tree = null;
	Vertex nephew;
	Vertex uncle;
	
	double edgeProposalWidthControlVariable;
	EdgeMove nephewEdgeMove;
	EdgeMove parentEdgeMove;
	EdgeMove uncleEdgeMove;
	
	PriorDistribution<Double> edgePrior;
	
	public TopologyMove (McmcModule m, PriorDistribution<Double> pr, double propVar, String n) {
		owner = m;
		edgePrior = pr;
		edgeProposalWidthControlVariable = propVar;
		name = n;
		autoTune = false;
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
		if (vnum <= 3) { // then we can't do the nephew-uncle swap
			return;
		}
		
		int rnd = Utils.generator.nextInt(vnum);
		int vertId = tree.getTopVertexId(rnd);
		// Need to choose a valid nephew, i.e. not the root or either of its children
		while (vertId >= 0) { // while nephew is invalid
			rnd = Utils.generator.nextInt(vnum); // select a new nephew
			vertId = tree.getTopVertexId(rnd);
		}
		
		nephew = tree.vertex[rnd];
		uncle = nephew.parent.brother();
		
		if (nephewEdgeMove == null) {
			nephewEdgeMove = new EdgeMove(
					owner,nephew.index,edgePrior,new MultiplicativeProposal(),
					edgeProposalWidthControlVariable,"nephewEdge");
		}
		else {
			nephewEdgeMove.setEdgeIndex(nephew.index);
		}
		if (parentEdgeMove == null) {
			parentEdgeMove = new EdgeMove(
					owner,nephew.parent.index,edgePrior,new MultiplicativeProposal(),
					edgeProposalWidthControlVariable,"parentEdge");
		}
		else {
			parentEdgeMove.setEdgeIndex(nephew.parent.index);
		}
		if (uncleEdgeMove == null) {
			uncleEdgeMove = new EdgeMove(
					owner,uncle.index,edgePrior,new MultiplicativeProposal(),
					edgeProposalWidthControlVariable,"uncleEdge");
		}
		else {
			uncleEdgeMove.setEdgeIndex(uncle.index);
		}
		
		nephewEdgeMove.copyState(externalState);
		parentEdgeMove.copyState(externalState);
		uncleEdgeMove.copyState(externalState);
		
		((CoreMcmcModule) owner).getModelExtMan().beforeTreeChange(tree, nephew);		
		// Should also do a beforeAlignChange here, but not obvious what to pass
		// as the selectedRoot argument.
		// The rest of the state copying is handled inside the Vertex
	}
	@Override
	public double proposal(Object externalState) {
		System.out.println("Before: "+tree.printedTree());
		double logProposalRatio = nephewEdgeMove.proposal(externalState);
		logProposalRatio += parentEdgeMove.proposal(externalState);
		logProposalRatio += uncleEdgeMove.proposal(externalState);
//		nephew.calcAllUp();
//		nephew.parent.calcAllUp();
//		uncle.calcAllUp();
		
		logProposalRatio += nephew.fastSwapWithUncle();
		// Below is another version, slow but slightly better mixing
		// double logProposalRatio = nephew.swapWithUncleAlignToParent();
		
		return logProposalRatio;
	}
	
	public void printChildren(Vertex v){
		if(v.left != null){
			System.out.println(v.index + " " + v.left.index);
			System.out.println(v.index + " " + v.right.index);
			printChildren(v.left);
			printChildren(v.right);
		}
	}

	public void printEdges(Vertex v){
		System.out.println(v.index + " " + v.edgeLength);
		if(v.left != null){
			printEdges(v.left);
			printEdges(v.right);
		}
	}
	
	
	@Override
	public double logPriorDensity(Object externalState) {
		return 0.0;
	}
	@Override
	public void updateLikelihood(Object externalState) {
		System.out.println("After: "+tree.printedTree());
		System.out.print("LogLikelihood before: "+owner.curLogLike+"\t");
		owner.setLogLike(((CoreMcmcModule) owner).getModelExtMan().logLikeTreeChange(tree, nephew));
		System.out.println("LogLikelihood after: "+owner.curLogLike);
	}
	@Override
	public void restoreState(Object externalState) {
		
		// Note these are restored in the reverse order to the proposal
		uncle.fastSwapBackUncle();
		// If using the alternative move:
        // uncle.swapBackUncleAlignToParent();
		
		uncleEdgeMove.restoreState(externalState);
		parentEdgeMove.restoreState(externalState);
		nephewEdgeMove.restoreState(externalState);
		
		
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
