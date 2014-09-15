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
		this.tree.heat = heat;
		isParallel = true;
	}
	protected boolean isMaster() {
		return MPIUtils.isMaster(rank);
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
			// Send state.
			State state = getState();
			MPIStateSend(state);
		} else if (!isColdChain() && MPIUtils.isMaster(rank)) {
			// Receive state.
			State state = MPIStateReceieve(coldChainLocation);
			postprocMan.newSample(coreModel,state, no, total);
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
				rank + 1, noOfProcesses, tree.heat);
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
	protected boolean isColdChain() {
		return tree.heat == 1.0d;
	}
		
	protected State MPIStateReceieve(int peer) {
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

	protected void MPIStateSend(State state) {

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
