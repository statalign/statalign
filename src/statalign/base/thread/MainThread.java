package statalign.base.thread;

import java.io.File;

import statalign.base.MainManager;
import statalign.base.Mcmc;
import statalign.base.Tree;
import statalign.base.TreeAlgo;
import statalign.io.RawSequences;
import statalign.ui.MainFrame;

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
			
			owner.inputData.pars.fixAlign = owner.inputData.useAlign >= 2;
			owner.inputData.pars.fixTopology = owner.inputData.useTree >= 2;
			owner.inputData.pars.fixEdge = owner.inputData.useTree >= 3;
			owner.inputData.pars.doReportDuringBurnin = owner.inputData.doReportDuringBurnin;

			// initialise model extension plugins for run
			owner.modelExtMan.initRun(owner.inputData);
			
			RawSequences seqs = owner.inputData.seqs;
			
			if(owner.frame != null) {
				owner.frame.statusText.setText(" Generating initial tree and alignment...");
			}

			System.out.println("\nPreparing initial tree and alignment...\n");

			// remove gaps and whitespace
			//owner.inputData.title = new File(owner.fullPath).getName();
			owner.inputData.setBaseFile(new File(owner.fullPath));			
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
			String[] names = seqs.getSeqnames().toArray(new String[seqs.size()]);

			TreeAlgo treeAlgo = new TreeAlgo();
			Tree tree;
			if(owner.inputData.useTree > 0) {				
				tree = treeAlgo.rearrangeTree(owner.inputData.tree, names);
				treeAlgo.addAlignSeqsToTree(tree, nongapped, names,
						owner.inputData.model, new File(owner.fullPath).getName());
			} else {
				tree = treeAlgo.buildNJTree(nongapped, seqs.getSeqnames().toArray(new String[seqs.size()]), 	
						owner.inputData.model, new File(owner.fullPath).getName());
			}
			System.out.println("Initial tree: "+tree.printedTree()+"\n");
//			Tree tree = new Tree(nongapped, seqs.seqNames.toArray(new String[seqs.seqNames.size()]), 	
//					owner.inputData.model,
//					owner.inputData.model.attachedScoringScheme,
//					new File(owner.fullPath).getName());
			
			Mcmc mcmc = new Mcmc(tree, owner.inputData.pars, owner.postProcMan, owner.modelExtMan);
			int errorCode = mcmc.doMCMC();
			owner.finished(errorCode, null);
			System.out.println(errorCode == 0 ? "Ready." : "Stopped.");
		} catch(StoppedException e) {
			if (owner.frame != null) {
				owner.frame.statusText.setText(MainFrame.IDLE_STATUS_MESSAGE);
			}
			owner.finished(2, null);
			System.out.println("Stopped.");
		} catch (IllegalArgumentException e) {
			e.printStackTrace();
			printStateInfo();

			String msg = e.getMessage();
			if(msg != null)
				System.err.println("Plugin error: "+msg);
		} catch(Exception e) {
			e.printStackTrace();
			printStateInfo();

			if(owner.frame != null) {				
				System.out.println("Here is the error: " + e.getClass());
				owner.frame.statusText.setText(MainFrame.IDLE_STATUS_MESSAGE);
				//ErrorMessage.showPane(owner.frame,e.getLocalizedMessage(),true);
			}			
			owner.finished(-1, e);
		}
	}
	
	private void printStateInfo() {
		System.err.println(owner.modelExtMan.getMcmc().tree.printedTree());	
		owner.modelExtMan.getMcmc().tree.checkPointers();
		owner.modelExtMan.getMcmc().tree.root.recomputeCheckLogLike();
		System.err.println(owner.fullPath);		
		String align[] = owner.modelExtMan.getMcmc().getState().getFullAlign();
		for (int i=0; i<align.length; i++) {
			System.err.println(align[i]);				
		}
	}
}
