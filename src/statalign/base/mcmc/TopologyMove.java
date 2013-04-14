package statalign.base.mcmc;

import statalign.base.Tree;
import statalign.base.Utils;
import statalign.base.Vertex;
import statalign.mcmc.McmcModule;
import statalign.mcmc.McmcMove;
import java.io.FileWriter;
import java.io.IOException;
import statalign.mcmc.MultiplicativeProposal;
import statalign.mcmc.PriorDistribution;
import statalign.mcmc.UniformProposal;

public class TopologyMove extends McmcMove {
	
	FileWriter topMoves;
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
		if(Utils.DEBUG){
			try{
				topMoves = new FileWriter("topMoves.txt");
				topMoves.write("origInd \t origSeq \t propInd \t propSeq \t Hastings \t logAccept \n");
			} catch (IOException e){}
		}
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
			nephewEdgeMove.setMinValue(Utils.MIN_EDGE_LENGTH);
		}
		else {
			nephewEdgeMove.setEdgeIndex(nephew.index);
		}
		if (parentEdgeMove == null) {
			parentEdgeMove = new EdgeMove(
					owner,nephew.parent.index,edgePrior,new MultiplicativeProposal(),
					edgeProposalWidthControlVariable,"parentEdge");
			parentEdgeMove.setMinValue(Utils.MIN_EDGE_LENGTH);
		}
		else {
			parentEdgeMove.setEdgeIndex(nephew.parent.index);
		}
		if (uncleEdgeMove == null) {
			uncleEdgeMove = new EdgeMove(
					owner,uncle.index,edgePrior,new MultiplicativeProposal(),
					edgeProposalWidthControlVariable,"uncleEdge");
			uncleEdgeMove.setMinValue(Utils.MIN_EDGE_LENGTH);
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
		System.out.println("Before Topology: "+tree.printedTree());

//		if(Utils.DEBUG){
//			System.out.println("Before:");
//			System.out.println("Likelihood:");
//			System.out.println(owner.curLogLike);
//			double[] params = tree.getState().indelParams;
//			for(int i = 0; i < params.length; i++)
//				System.out.println(params[i]);
//			String[] fullAlign = tree.getState().getFullAlign();
//			for(int i = 0; i < fullAlign.length; i++)
//				System.out.println(fullAlign[i]);
//			String printTree = tree.printedTree();
//			Vertex.printChildren(tree.root);
//			Vertex.printEdges(tree.root);
//			System.out.println(printTree);
//			System.out.println(tree.root.indelLogLike);
//			System.out.println(tree.root.orphanLogLike);
//			try{
//				topMoves.write(tree.root.indelLogLike + "\t" + tree.root.orphanLogLike + "\t");
//			} catch(IOException e){}
//		}
		
//		double logProposalRatio = nephewEdgeMove.proposal(externalState);
//		System.out.println("logProposalRatio after nephew = "+logProposalRatio);
//		logProposalRatio += parentEdgeMove.proposal(externalState);
//		System.out.println("logProposalRatio after parent = "+logProposalRatio);
//		logProposalRatio += uncleEdgeMove.proposal(externalState);
//		System.out.println("logProposalRatio after uncle = "+logProposalRatio);
		
//		logProposalRatio += nephew.fastSwapWithUncle();
		//double logProposalRatio = nephew.fastSwapWithUncle();
		String[] s = nephew.printedAlignment();
		System.out.println(nephew.index+" "+s[1]+"\n"+nephew.parent.index+" "+s[0]);
		s = uncle.printedAlignment();
		System.out.println(uncle.index+" "+s[1]+"\n"+uncle.parent.index+" "+s[0]);
		//double logProposalRatio = nephew.fastSwapWithUncle();
		double logProposalRatio = nephew.swapWithUncleAlignToParent();
		System.out.println("logProposalRatio after swap = "+logProposalRatio);
		System.out.println("After  Topology: "+tree.printedTree());
		s = nephew.printedAlignment();
		System.out.println(s[0]+"\n"+s[1]);

		if (logProposalRatio == Double.POSITIVE_INFINITY) {
			System.out.println("Likelihood: "+tree.getLogLike());
			throw new RuntimeException("Let's stop now and have a rest.");
		}
//		if(Utils.DEBUG){
//			System.out.println("Proposed:");
//			String[] fullAlign = tree.getState().getFullAlign();
//			for(int i = 0; i < fullAlign.length; i++)
//				System.out.println(fullAlign[i]);
//			String printTree = tree.printedTree();
//			System.out.println(printTree);
//			Vertex.printChildren(tree.root);
//			Vertex.printEdges(tree.root);
//			System.out.println("logProposalRatio:");
//			System.out.println(logProposalRatio);
//			try{
//				topMoves.write(tree.root.indelLogLike + "\t" + tree.root.orphanLogLike + "\t" + logProposalRatio + "\n");
//			} catch(IOException e){}
//			
//		}
		
		return logProposalRatio;
	}
	
	
	@Override
	public double logPriorDensity(Object externalState) {
		return 0.0;
	}
	@Override
	public void updateLikelihood(Object externalState) {
//		System.out.println("After: "+tree.printedTree());
//		System.out.print("LogLikelihood before: "+owner.curLogLike+"\t");
		owner.setLogLike(((CoreMcmcModule) owner).getModelExtMan().logLikeTreeChange(tree, nephew));
//		System.out.println("LogLikelihood after: "+owner.curLogLike);
	}
	@Override
	public void restoreState(Object externalState) {
		
		// Note these are restored in the reverse order to the proposal
	//	uncle.fastSwapBackUncle();
		// If using the alternative move:
       uncle.swapBackUncleAlignToParent();
		
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
		//throw new RuntimeException("Let's stop now and have a rest.");
	}
	
	 
	
}
