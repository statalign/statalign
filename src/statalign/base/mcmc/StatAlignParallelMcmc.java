package statalign.base.mcmc;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Random;

import org.apache.commons.math3.random.Well19937c;

import mpi.MPI;
import statalign.MPIUtils;
import statalign.base.MCMCPars;
import statalign.base.Mcmc;
import statalign.base.State;
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
public class StatAlignParallelMcmc extends StatAlignMcmc {
	
	/** The number of processes. */
	protected int noOfProcesses;

	/** The rank of the process. */
	protected int rank;
	
	protected boolean verbose = false;
	
	/** The random number generator used for swapping. */
	protected Random swapGenerator;

	public StatAlignParallelMcmc(Tree tree, MCMCPars pars, PostprocessManager ppm, ModelExtManager modelExtMan,
			int noOfProcesses, int rank, double heat) {
		super(tree,pars,ppm,modelExtMan);
		this.noOfProcesses = noOfProcesses;
		this.rank = rank;
		heatWhenHot = heat;
		isParallel = true;
	}
	protected boolean isMaster() {
		return MPIUtils.isMaster(rank);
	}
	protected boolean isColdChain() {
		return heat == 1.0d;
	}	
	protected void postprocSample(int no, int total) {
		int[] ranks = new int[] { (isColdChain() ? rank : 0) };
		int[] coldChainLoc = new int[1];
		int coldChainLocation = -1;
		MPI.COMM_WORLD.Reduce(ranks, 0, coldChainLoc, 0, 1, MPI.INT, MPI.SUM, 0);
		coldChainLocation = coldChainLoc[0];

		// TODO: Remove - for debugging purposes
		if (MPIUtils.isMaster(rank)) {
			MPIUtils.println(rank, "Cold chain is at: " + coldChainLocation);
		}
		
		postprocMan.newSample(coreModel,getState(), no, total);			
	}
	protected void beforeMCMC() {
		String str = String.format(
				"Starting MCMC chain no. %d/%d (heat: %.2f)\n\n", 
				rank + 1, noOfProcesses, heat);
		MPIUtils.println(rank, str);
		swapGenerator = new Random(mcmcpars.swapSeed);
		// the above ensures that all chains are using the same
		// uniform random numbers to decide whether to switch temp

		//Utils.generator = new Well19937c(mcmcpars.seed + rank);
		Utils.generator = new Well19937c(mcmcpars.seed); 
		// seed is already changed on a per-chain basis inside MainThread.java 
	}
	protected void doSwap() {
		int chainA, chainB;
		chainA = swapGenerator.nextInt(noOfProcesses);
		do {
			chainB = swapGenerator.nextInt(noOfProcesses);
		} while (chainA == chainB);
		
		// TODO: Should only attempt to swap chains with similar heat?		
		if (verbose) {
			System.out.printf("chainA = %d, chainB = %d\n",chainA, chainB);
		}

		double logU = Math.log(swapGenerator.nextDouble());
		// this has to be outside the if condition below,
		// to ensure that all chains advance their generators
		// at the same rate.

		if (rank == chainA || rank == chainB) {
						
			double[] myStateInfo = new double[3];
			myStateInfo[0] = totalLogLike;
			myStateInfo[1] = modelExtMan.totalLogPrior(tree);
			//myStateInfo[1] = coreModel.totalLogPrior(tree) + modelExtMan.totalLogPrior(tree);
			myStateInfo[2] = heat;

			double[] partnerStateInfo = MPI_exchange(chainA,chainB,myStateInfo); 
			
			if (verbose) {
				System.out
				.printf("[Worker %d] Heat: [%f] - Sent: [%f,%f,%f] - Recv: [%f,%f,%f]\n",
						rank, heat, myStateInfo[0], myStateInfo[1],
						myStateInfo[2], partnerStateInfo[0],
						partnerStateInfo[1], partnerStateInfo[2]);
			}

			double myLogLike = myStateInfo[0];
			double myLogPrior = myStateInfo[1];
			double myTemp = myStateInfo[2];
			double hisLogLike = partnerStateInfo[0];
			double hisLogPrior = partnerStateInfo[1];
			double hisTemp = partnerStateInfo[2];			
			
			if (heat == hisTemp) return;
						
			double logRatio = myTemp * (hisLogLike + hisLogPrior) + hisTemp
					* (myLogLike + myLogPrior);
			logRatio -= hisTemp * (hisLogLike + hisLogPrior) + myTemp
					* (myLogLike + myLogPrior);			
			
			if (verbose) {
				MPIUtils.println(rank,
						"logU: " + logU);
				MPIUtils.println(rank, "acceptance:           "
						+ logRatio);				
			}
			if (logU < logRatio) {				
				if (verbose) {
					MPIUtils.println(rank,
							"Just swapped heat with my partner. New heat: "
									+ hisTemp);
				}
				
				// Swap the heats
				heat = hisTemp;
								
				// Swap the proposal width variables
				double[] myProposalWidths, partnerProposalWidths; 
				
				myProposalWidths = modelExtMan.getProposalWidths();
				partnerProposalWidths = MPI_exchange(chainA,chainB,myProposalWidths);							
				modelExtMan.setProposalWidths(partnerProposalWidths);
				
				myProposalWidths = coreModel.getProposalWidths();
				partnerProposalWidths = MPI_exchange(chainA,chainB,myProposalWidths);							
				coreModel.setProposalWidths(partnerProposalWidths);
			}			
		}
	}	
	
	public double[] MPI_exchange(int rankA, int rankB, double[] x) {
		mpi.Request send, recieve;
		double[] y = new double[x.length];
		if (rank == rankA) {
			send = MPI.COMM_WORLD.Isend(x,0,x.length, MPI.DOUBLE,rankB,0);
			recieve = MPI.COMM_WORLD.Irecv(y,0,x.length,MPI.DOUBLE,rankB, 1);
		} else {
			send = MPI.COMM_WORLD.Isend(x,0,x.length,MPI.DOUBLE,rankA, 1);
			recieve = MPI.COMM_WORLD.Irecv(y,0,x.length,MPI.DOUBLE,rankA, 0);
		}		
		mpi.Request.Waitall(new mpi.Request[] { send, recieve });
		return y;
	}
}
