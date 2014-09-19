package statalign.mcmc;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;

import statalign.base.Mcmc;
import statalign.base.Tree;
import statalign.base.Utils;

/**
 * 
 * Generic class for a group of McmcMove objects. Each ModelExtension
 * is an instance of this class, as is the CoreMcmcModule
 * 
 * @author herman
 *
 */

public abstract class McmcModule {
	
	protected Mcmc mcmc;
	protected String moduleName = "";
	public boolean logParametersToFile = false;
	FileWriter parameterLog;
	
	public McmcModule() { }	
	public McmcModule(String name) {
		moduleName = name;
		logParametersToFile = true;		
	}
	public void setOutputFile(String baseFileName) {
		try {
			parameterLog = new FileWriter(baseFileName+moduleName+".params");
			//System.out.println(baseFileName+moduleName+".params");
		}
		catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public void setMcmc(Mcmc m) {
		mcmc = m;
	}
	public boolean printExtraInfo = false;	
	
	public boolean isFirstHalfBurnin() {
		return mcmc.firstHalfBurnin;
	}
	public boolean isBurnin() {
		return mcmc.burnin;
	}
	/** Current log-likelihood contribution */
	public double curLogLike = 0;
	
	protected List<McmcMove> mcmcMoves = new ArrayList<McmcMove>();
	protected List<Integer> mcmcMoveWeights = new ArrayList<Integer>();
	protected List<Integer> mcmcMoveWeightIncrements = new ArrayList<Integer>();
	
	public int getParamChangeWeight() {
		int w = 0;
		for (int i=0; i<mcmcMoveWeights.size(); i++) {
			w += mcmcMoveWeights.get(i);
		}
		return w;
	}
	public void setWeight(String name, int weight) {
		for (int i=0; i<mcmcMoves.size(); i++) {
			if (mcmcMoves.get(i).name.contains(name)) {
				mcmcMoveWeights.set(i, weight);
				if (printExtraInfo) System.out.println("Move \""+mcmcMoves.get(i).name+"\" now has weight "+weight);
			}
		}
	}
	public void addMcmcMove(McmcMove m, int weight) {
		mcmcMoves.add(m);
		mcmcMoveWeights.add(weight);
		mcmcMoveWeightIncrements.add(0);
	}
	public void addMcmcMove(McmcMove m, int weight, int increment) {
		mcmcMoves.add(m);
		mcmcMoveWeights.add(weight);
		mcmcMoveWeightIncrements.add(increment);
	}
	public List<McmcMove> getMcmcMoves() {
		return mcmcMoves;
	}
	public void setAllMovesNotProposed() {
		for (McmcMove mcmcMove : mcmcMoves) {
			mcmcMove.moveProposed = false;
		}
	}
	public void zeroAllMoveCounts() {
		for (McmcMove mcmcMove : mcmcMoves) {
			mcmcMove.proposalCount = 0;
			mcmcMove.acceptanceCount = 0;
			mcmcMove.lowCounts = 0;
		}
	}
	public McmcMove getMcmcMove(String name) {
		for (McmcMove mcmcMove : mcmcMoves) {
			if (mcmcMove.name.equals(name)) {
				return mcmcMove;
			}
		}
		throw new RuntimeException("McmcMove "+name+" not found.");
	}
	public String getMcmcInfo() {
		String info = "";
		for (McmcMove mcmcMove : mcmcMoves) {
			String infoFormat = "%-24s%8s%8d%8d%8.4f%8.4f\n";
			info += String.format(Locale.US, infoFormat,
					mcmcMove.name,
					Utils.convertTime(mcmcMove.getTime()),
					mcmcMove.proposalCount,
					mcmcMove.getTime()/(mcmcMove.proposalCount>0 ? mcmcMove.proposalCount : 1),
					mcmcMove.acceptanceRate(),
					mcmcMove.proposalWidthControlVariable);
		}
		return info;
	}
	public String getSummaryInfo() {
		String info = "Acceptance rates: ";
		for (McmcMove m : mcmcMoves) {
			info += m.name+": "+String.format(Locale.US, "%f ", m.acceptanceRate());
		}
		return info;
	}
	public void printParameters() {		
		if (logParametersToFile) {
			String params = "";
			for (McmcMove m : mcmcMoves) {			
				if (m.printableParam) {
					if (params != "") params += ", ";
					params += m.getParameterString();			
				}
			}
			if (params != null) {
				try {
					parameterLog.write(params+"\n");
				}
				catch (IOException e) {
					e.printStackTrace();
				}
			}
		}
	}
	/**
	 * Called before the start of MCMC sampling, but after the initial tree, alignment etc. have been
	 * generated. Override to initialise data structures etc.
	 * @param tree the starting tree
	 */
	public void beforeSampling(Tree tree) {
		if (logParametersToFile) {
			String paramNames = "";
			for (McmcMove m : mcmcMoves) {			
				if (m.printableParam) {
					if (paramNames != "") paramNames += ", ";
					paramNames += m.getNameString();
				}
			}
			if (paramNames != null) {
				try {
					if (parameterLog == null) System.out.println("null log");
					parameterLog.write(paramNames+"\n");
				}
				catch (IOException e) {
					e.printStackTrace();
				}
			}
		}
	}
	
	public void afterSampling() {
		if (parameterLog != null) {
			try {
				parameterLog.close();
			}
			catch (IOException e) {
				e.printStackTrace();
			}
		}
		for (McmcMove m : mcmcMoves) {
			m.printInfo();
		}
	}
	
	/**
	 * This should return the log of the model's contribution to the likelihood, it will be added on to
	 * the log-likelihood of the current point in the MCMC state space. Normally it will be called once at the
	 * initialisation of the MCMC process and from then on once in each MCMC step, when proposing any change.
	 * In debug mode, will be called more often (including after proposed changes) to ensure consistency.
	 * @param tree current tree
	 * @return log of model extension likelihood, conditional on current tree, alignment and params
	 */
	public abstract double logLikeFactor(Tree tree);
	
	public double getLogLike() {
		return curLogLike;
	}
	public void setLogLike(double ll) {
		curLogLike = ll;
	}

	/**
	 * This should return the log of the total prior calculated for the model parameters. It is only used
	 * in parallel mode when proposing swaps between chains. By default returns 0.
	 */
	public double logPrior(Tree tree) {
		return 0;
	}
	
	public boolean proposeParamChange(Tree tree) {
		int selectedMoveIndex = Utils.weightedChoose(mcmcMoveWeights);
		McmcMove selectedMove = mcmcMoves.get(selectedMoveIndex); 
		selectedMove.move(tree);
		return selectedMove.lastMoveAccepted;
	}
	
	public void modifyProposalWidths() {
		for (McmcMove m : mcmcMoves) {
			if (!m.autoTune) { continue; }
			if (m.proposalCount > Utils.MIN_SAMPLES_FOR_ACC_ESTIMATE) {
				if (m.acceptanceRate() < m.minAcceptance) {
					m.proposalWidthControlVariable *= m.spanMultiplier;
					m.proposalCount = 0;
					m.acceptanceCount = 0;
					m.lowCounts++;
				}
				else if (m.acceptanceRate() > m.maxAcceptance && 
						m.proposalWidthControlVariable <= m.maxProposalWidthControlVariable) {
					m.proposalWidthControlVariable /= m.spanMultiplier;
					m.proposalCount = 0;
					m.acceptanceCount = 0;
					m.lowCounts = 0;
				}
				else m.lowCounts = 0;
			}
		}
	}
	public boolean isParamChangeAccepted(double logProposalRatio,McmcMove m) {
		return mcmc.isParamChangeAccepted(logProposalRatio,m);
	}
	public void incrementWeights() {
		for (int i=0; i<mcmcMoves.size(); i++) {
			if (mcmcMoveWeightIncrements.get(i) != 0) {
				mcmcMoveWeights.set(i,mcmcMoveWeights.get(i)+
						mcmcMoveWeightIncrements.get(i));				
				if (printExtraInfo) System.out.println("Move \""+mcmcMoves.get(i).name+"\" now has weight "+mcmcMoveWeights.get(i));
			}
		}
	}

	public void afterFirstHalfBurnin() { 
		for (McmcMove m : mcmcMoves) m.afterFirstHalfBurnin();
	}
	public void afterBurnin() { 
		for (McmcMove m : mcmcMoves) m.afterBurnin();
	}
}
