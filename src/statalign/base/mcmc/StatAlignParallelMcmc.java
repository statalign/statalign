package statalign.base.mcmc;

import java.util.Arrays;
import java.util.Random;

import org.apache.commons.math3.random.Well19937c;

import mpi.MPI;
import statalign.base.MCMCPars;
import statalign.base.Tree;
import statalign.base.Utils;
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
			int noOfProcesses, int rank, double beta) {
		super(tree,pars,ppm,modelExtMan);
		this.noOfProcesses = noOfProcesses;
		this.rank = rank;
		betaWhenHot = beta;
		isParallel = true;
	}
	protected boolean isMaster() {
		return rank == 0;
	}
	protected boolean isColdChain() {
		return beta == 1.0d;
	}	
	protected void postprocSample(int no, int total) {
//		int[] ranks = new int[] { (isColdChain() ? rank : 0) };
//		int[] coldChainLoc = new int[1];
//		int coldChainLocation = -1;
//		MPI.COMM_WORLD.Reduce(ranks, 0, coldChainLoc, 0, 1, MPI.INT, MPI.SUM, 0);
//		coldChainLocation = coldChainLoc[0];
		
		postprocMan.newSample(coreModel,getState(), no, total);			
	}
	protected void beforeMCMC() {
		String str = String.format(
				"Starting MCMC chain no. %d/%d (beta: %.2f)\n\n", 
				rank + 1, noOfProcesses, beta);
		System.out.println(str);
		swapGenerator = new Random(mcmcpars.swapSeed);
		// the above ensures that all chains are using the same
		// uniform random numbers 

		//Utils.generator = new Well19937c(mcmcpars.seed + rank);
		Utils.generator = new Well19937c(mcmcpars.seed); 
		// seed is already changed on a per-chain basis inside MainThread.java 
	}

	protected void doSwap() {
		int chainA, chainB;
		
		double logU = Math.log(swapGenerator.nextDouble());

		// each chain starts with the same seed, so they all
		// generate the same sequence of integers
		chainA = swapGenerator.nextInt(noOfProcesses-1);
		do {
			chainB = swapGenerator.nextInt(noOfProcesses);
		} while(chainB==chainA);

		if (verbose) {
			System.out.printf("chainA = %d, chainB = %d\n",chainA, chainB);
		}

		if (rank == chainA || rank == chainB) {
						
			double[] myStateInfo = new double[3];
			myStateInfo[0] = totalLogLike;
			myStateInfo[1] = modelExtMan.totalLogPrior(tree);
			//myStateInfo[1] = coreModel.totalLogPrior(tree) + modelExtMan.totalLogPrior(tree);
			myStateInfo[2] = beta;

			double[] partnerStateInfo = MPI_exchange(chainA,chainB,myStateInfo); 
			
			if (verbose) {
				System.out
				.printf("[Worker %d] Beta: [%f] - Sent: [%f,%f,%f] - Recv: [%f,%f,%f]\n",
						rank, beta, myStateInfo[0], myStateInfo[1],
						myStateInfo[2], partnerStateInfo[0],
						partnerStateInfo[1], partnerStateInfo[2]);
			}

			double myLogLike = myStateInfo[0];
			double myLogPrior = myStateInfo[1];
			double myBeta = myStateInfo[2];
			double hisLogLike = partnerStateInfo[0];
			double hisLogPrior = partnerStateInfo[1];
			double hisBeta = partnerStateInfo[2];
			
			if (beta == hisBeta) return;
						
			double logRatio = myBeta * (hisLogLike + hisLogPrior) + hisBeta
					* (myLogLike + myLogPrior);
			logRatio -= hisBeta * (hisLogLike + hisLogPrior) + myBeta
					* (myLogLike + myLogPrior);			
			
			if (verbose) {
				System.out.println("logU: " + logU);
				System.out.println("acceptance: " + logRatio);				
			}
			if (logU < logRatio) {				
				if (verbose) {
					System.out.println("Just swapped beta with my partner. New beta: "+ hisBeta);
				}
				
				// Swap the betas
				beta = hisBeta;
								
				// Swap the proposal width variables
				double[] myProposalWidths, partnerProposalWidths; 
				
				double[] modelExtProposalWidths = modelExtMan.getProposalWidths();				
				double[] coreModelProposalWidths = coreModel.getProposalWidths();
				
				// Combine the two arrays into one, to reduce the number of communication
				// operations required.
				myProposalWidths = new double[modelExtProposalWidths.length + 
				                              coreModelProposalWidths.length];
				System.arraycopy(modelExtProposalWidths, 0, myProposalWidths, 
						0, modelExtProposalWidths.length);
				System.arraycopy(coreModelProposalWidths, 0, myProposalWidths, 
						modelExtProposalWidths.length, coreModelProposalWidths.length);

				partnerProposalWidths = MPI_exchange(chainA,chainB,myProposalWidths);							
				modelExtMan.setProposalWidths(Arrays.copyOfRange(partnerProposalWidths,0,modelExtProposalWidths.length));											
				coreModel.setProposalWidths(Arrays.copyOfRange(partnerProposalWidths,modelExtProposalWidths.length,
						modelExtProposalWidths.length+coreModelProposalWidths.length));
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
