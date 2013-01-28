package statalign.model.ext.plugins.structalign;

import statalign.base.Utils;
import statalign.model.ext.plugins.StructAlign;

public class TranslationMove extends RotationOrTranslationMove {
	
	public TranslationMove (StructAlign s, String n) {
		owner = s;
		name = n;
	}

	public double proposal(Object externalState) {
		double[] shift = new double[3];
		for(int i = 0; i < 3; i++)
			shift[i] = Utils.generator.nextGaussian() * owner.xlatP;
		for(int l = 0; l < subtreeLeaves.size(); l++){
			int j = subtreeLeaves.get(l);
			for(int i = 0; i < 3; i++)
				owner.xlats[j][i] += shift[i];  
		}	
		return 0;
		// logProposalRatio is 0 because prior is uniform and proposal is symmetric
	}

}
