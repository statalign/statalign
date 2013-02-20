package statalign.base.thread;

import statalign.base.MainManager;
import statalign.base.Mcmc;
import statalign.base.Tree;
import statalign.io.RawSequences;

/**
 * The main (suspendable) thread for background MCMC calculation.
 *
 * @author novak
 */
public class MainThread extends StoppableThread {
	
	/**
	 * Reference to the (singleton) MainManager object encapsulating all settings and data
	 * that an MCMC run depends on.
	 */
	public MainManager owner;
	
	/**
	 * Constructs a new MainThread that can be used to fire a background MCMC calculation.
	 * 
	 * @param owner Reference to the MainManager object.
	 */
	public MainThread(MainManager owner) {
		this.owner = owner;
	}
	/**
	 * Start background MCMC calculation.
	 */
	@Override
	public synchronized void run() {
		try {

			owner.postProcMan.initRun(owner.inputData);
			
			RawSequences seqs = owner.inputData.seqs;
			
			if(owner.frame != null) {
				owner.frame.statusText.setText(" Generating initial tree and alignment...");
			}

			System.out.println("\nPreparing initial tree and alignment...\n");

			// remove gaps and whitespace
			String[] nongapped = new String[seqs.size()];
			StringBuilder builder = new StringBuilder();
			int i, j;
			char ch;
			for(i = 0; i < nongapped.length; i++) {
				builder.setLength(0);
				String seq = seqs.getSequence(i);
				for(j = 0; j < seq.length(); j++) {
					ch = seq.charAt(j);
					if(Character.isWhitespace(ch) || ch == '-')
						continue;
					builder.append(ch);
				}
				nongapped[i] = builder.toString();
			}
			
			Tree tree = new Tree(nongapped, seqs.getSeqnames().toArray(new String[seqs.size()]), 	
					owner.inputData.model,
					owner.inputData.model.attachedScoringScheme);
			Mcmc mcmc = new Mcmc(tree, owner.inputData.pars, owner.postProcMan);
			int errorCode = mcmc.doMCMC();
			owner.postProcMan.finalizeRun();
			owner.finished(errorCode, null);
			System.out.println(errorCode == 0 ? "Ready." : "Stopped.");
		} catch (StoppedException e) {
			// stopped during initial alignment
			owner.postProcMan.finalizeRun();
			owner.finished(2, null);
			System.out.println("Stopped.");
		} catch(Exception e) {
			owner.postProcMan.finalizeRun();
			e.printStackTrace();
			owner.finished(-1, e);
		}
		
	}
}
