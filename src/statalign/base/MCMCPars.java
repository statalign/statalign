package statalign.base;

/**
 * 
 * This is a container class containing the MCMC parameters described below.
 * 
 * @author miklos, novak
 *
 */
public class MCMCPars {
	/**
	 * 
	 * Number of burn-in steps
	 * 
	 */
	public int burnIn;
	/**
	 * Number of cycles total after burn-in period
	 * 
	 */
	public int cycles;       
	
	/**
	 * One sampling after each sampRate cycles
	 */
	public int sampRate;

	/**
	 * The seed for the random number generator
	 */
	public long seed;
	
	/**
	 * The swap seed for the swap random number generator.
	 */
	public long swapSeed;
	
	/**
	 * How often to propose swaps between chains (per cycle). 
	 */
	public int swapRate;
	
	/**
	 * Difference in temperature between adjacent parallel chains.
	 */
	public double tempDiff;
	
	/**
	 * Minimum chain temperature.
	 */
	public double minTemp;
	
	/**
	 * Length of the period in which the MCMC chain should be allowed to 
	 * drift randomly, with all moves accepted, before starting the 
	 * burnin. This can be useful for generating different starting 
	 * configurations to test convergence.
	 */
	public int randomisationPeriod;
	
	public boolean fixAlign = false;
	public boolean fixTopology = false;
	public boolean fixEdge = false;
	/** If <code>true</true> then logging information is also printed during the burnin. */
	public boolean doReportDuringBurnin = false;
	
	/** If <code>true</true> then output is printed for each parallel chain only when its temperature is 1.0d. */
	public boolean reportOnlyColdChain = false;

	/**
	 * This constructor sets the values in the class
	 * 
	 * @param burnIn this.burnIn is set to this value.
	 * @param cycles this.cycles is set to this value.
	 * @param sampRate this.sampRate is set to this value.
	 */
	public MCMCPars(int burnIn, int cycles, int sampRate, long seed, long swapSeed, int swapRate) {
		this.burnIn = burnIn;
		this.cycles = cycles;
		this.sampRate = sampRate;
		this.swapRate = swapRate;
		this.seed = seed;
		this.swapSeed = swapSeed;
	}
	/**
	 * This constructor sets the values in the class
	 * 
	 * @param burnIn this.burnIn is set to this value.
	 * @param cycles this.cycles is set to this value.
	 * @param sampRate this.sampRate is set to this value.
	 * @param randomisationPeriod this.randomisationPeriod is set to this value.	 
	 */
	public MCMCPars(int burnIn, int cycles, int sampRate, long seed, long swapSeed, int swapRate, int randomisationPeriod, double tempDiff, double minTemp) {
		this.burnIn = burnIn;
		this.cycles = cycles;
		this.sampRate = sampRate;
		this.swapRate = swapRate;
		this.seed = seed;
		this.swapSeed = swapSeed;
		this.randomisationPeriod = randomisationPeriod;
		this.tempDiff = tempDiff;
		this.minTemp = minTemp;
	}
	
	/*public MCMCPars(int burnIn, int cycles, int sampRate, long seed) {
		this.burnIn = burnIn;
		this.cycles = cycles;
		this.sampRate = sampRate;
		this.seed = seed;
	}*/
	
}
