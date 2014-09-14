package tests;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;

import statalign.base.Utils;
import statalign.mcmc.MultiplicativeProposal;
import statalign.utils.BetaDistribution;

public class McmcTest {

	
	public static void main5(String[] args) {
		
		HashMap<Integer,String> test = new HashMap<Integer,String>();
		test.put(1, "1");
		test.put(2, "2");
		test.put(4, "4");
		
		for (int i=1; i<=4; i++) {
			String s = test.get(i);
			if (s==null) {
				test.put(i, ""+i);
			}
			System.out.println(s+" "+test.get(i));
		}
		
	}
	
	public static void main4(String[] args) {
		
		int nSamples = 10000;
		RandomGenerator master = new Well19937c(1);
		ArrayList<RandomGenerator> rands = new ArrayList<RandomGenerator>();
		//for (int i=0; i<10; i++) rands.add(new Well19937c(master.nextInt()));
		for (int i=0; i<10; i++) rands.add(new Well19937c(i));
		try {
			FileWriter output = new FileWriter("rngTest.txt");
			for (int i=0; i<nSamples; i++) {
				for (int j=0; j<10; j++) {
					output.write(rands.get(j).nextDouble()+"\t");
				}
				output.write("\n");
			}
			output.flush();
			output.close();
		} catch (IOException e) { throw new RuntimeException(e.toString());}
	}
	/**
	 * Code to test that the moves are being sampled according to the
	 * correct weights (i.e. that there is not some bias arising from
	 * the choose function and/or random number generator).
	 */
	public static void main3(String[] args) {
		final int rWeight = 8;
		final int lambdaWeight = 4;
		final int muWeight = 6;
		final int lambdaMuWeight = 6;
		final int phiWeight = 4;
		final int rhoWeight = 6;
		final int thetaWeight = 6;
		
		final int substWeight = 10;
		final int edgeWeight = 1; // per edge
		final int allEdgeWeight = 6;
		final int edgeWeightIncrement = 0; // Added after half of burnin
		final int alignWeight = 25;
		final int topologyWeight = 8;
		final int localTopologyWeight = 8;
		
		int[] weights = {
				rWeight,
				lambdaWeight, 
				rhoWeight,
				thetaWeight,				
				edgeWeight,
				edgeWeight,
				edgeWeight,
				edgeWeight,
				edgeWeight,
				edgeWeight,
				edgeWeight,
				allEdgeWeight,				
				alignWeight,
				topologyWeight,
				localTopologyWeight
		};
		
		final int nSamples = 10000;
		int[] samples = new int[nSamples];
		for (int i=0; i<nSamples; i++) {
			samples[i] = Utils.weightedChoose(weights);
		}
		try {
			FileWriter output = new FileWriter("sampleTest.txt");
			for (int i=0; i<nSamples; i++) {
				output.write(samples[i]+"\n");
			}
			output.flush();
			output.close();
		} catch (IOException e) { throw new RuntimeException(e.toString());}		
		System.out.println("Done.");
	}
		
	/**
	 * Code to test that various samplers are working correctly.
	 */
	public static void main2(String[] args) {
		
		MultiplicativeProposal proposalDistribution = new MultiplicativeProposal();
		//LogisticProposal proposalDistribution = new LogisticProposal();
		//GaussianProposal proposalDistribution = new GaussianProposal();
		//GammaDistribution n = new GammaDistribution(2,0.5);
		BetaDistribution n = new BetaDistribution(2,5);
		
		int ITERS = 100000;
		//double var = 3.5; // LogisticProposal
		//double var = 0.6; // GaussianProposal
		double var = 2.5; // MultiplicativeProposal
		double oldParam = 0.5; 
		double[] params = new double[ITERS]; 
		double MIN_PARAM = 0;
		double MAX_PARAM = 1; // Beta
		//double MAX_PARAM = Double.POSITIVE_INFINITY; // Gamma
		int acc = 0;
		
		for (int i=0; i<ITERS; i++) {						
			
			if (i>0) oldParam = params[i-1];
			
			proposalDistribution.updateProposal(var,oldParam);
			
		    double newParam = proposalDistribution.sample();						
		    
		    if (newParam <= MIN_PARAM || newParam > MAX_PARAM) {
		    	params[i] = oldParam;
		    	continue;
		    }
		    
		    double logProposalRatio = -proposalDistribution.logDensity(newParam);
		    
		    proposalDistribution.updateProposal(var,newParam);
		
		    logProposalRatio += proposalDistribution.logDensity(oldParam);
				    
		    double logLikeRatio = Math.log(n.density(newParam)) - Math.log(n.density(oldParam));
		    
		    //System.out.println(oldParam+"\t"+newParam+"\t"+logProposalRatio+"\t"+logLikeRatio);

		    if (Math.log(Utils.generator.nextDouble())< (logLikeRatio + logProposalRatio)) {
		    	params[i] = newParam;
		    	acc++;
		    }
		    else {
		    	params[i] = oldParam; 
		    }
		}
		System.out.println("Acceptance rate = "+((double)acc/(double)ITERS));
		try {
			FileWriter output = new FileWriter("multiplicativeTest.txt");
			for (int i=0; i<ITERS; i++) {
				output.write(params[i]+"\n");
			}
			output.close();
		} catch (IOException e) {}	

	}
}
