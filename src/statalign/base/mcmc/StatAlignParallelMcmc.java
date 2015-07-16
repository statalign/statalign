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

	/** The random number generator used for swapping. */
	protected Random swapGenerator;

	public StatAlignParallelMcmc(Tree tree, MCMCPars pars, PostprocessManager ppm, ModelExtManager modelExtMan,
			int noOfProcesses, int rank, double heat) {
		super(tree,pars,ppm,modelExtMan);
		this.noOfProcesses = noOfProcesses;
		this.rank = rank;
		this.heat = heat;
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

		if (isColdChain() && MPIUtils.isMaster(rank)) {
			// Sample normally.
			postprocMan.newSample(coreModel,getState(), no, total);			
			
		} else if (isColdChain() && !MPIUtils.isMaster(rank)) {
			/*
			 * TODO: We don't want to be passing around the full state for each postprocessing
			 * plugin, in the case where the cold chain is in a non-master process, 
			 * because we won't know in general what needs to be passed, and it'll be very slow.
			 * 
			 * Perhaps easiest to just write to separate files for each chain, although
			 * need to figure out how this can be achieved via MPI (can slave processes
			 * do IO without using MPI-IO?)
			 */
			postprocMan.newSample(coreModel,getState(), no, total);
		}
		if (MPIUtils.isMaster(rank)) {
			try {
				postprocMan.logFile.write("Cold chain location: " + coldChainLocation + "\n");
			}
			catch (IOException e) {
				e.printStackTrace();
			}
		}
	}
	protected void beforeMCMC() {
		String str = String.format(
				"Starting MCMC chain no. %d/%d (heat: %.2f)\n\n", 
				rank + 1, noOfProcesses, heat);
		MPIUtils.println(rank, str);
		swapGenerator = new Random(mcmcpars.swapSeed);
		Utils.generator = new Well19937c(mcmcpars.seed + rank);
	}
	protected void doSwap() {
		int swapA, swapB;
		swapA = swapGenerator.nextInt(noOfProcesses);
		do {
			swapB = swapGenerator.nextInt(noOfProcesses);
		} while (swapA == swapB);

		System.out.printf("SwapNo: %d - SwapA: %d - SwapB: %d\n",swapA, swapB);

		double swapAccept = swapGenerator.nextDouble();

		if (rank == swapA || rank == swapB) {
			double[] myStateInfo = new double[3];
			myStateInfo[0] = totalLogLike;
			myStateInfo[1] = modelExtMan.totalLogPrior(tree);
			//myStateInfo[1] = coreModel.totalLogPrior(tree) + modelExtMan.totalLogPrior(tree);
			myStateInfo[2] = heat;

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
}
