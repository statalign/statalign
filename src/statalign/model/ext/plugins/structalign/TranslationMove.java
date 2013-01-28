package statalign.model.ext.plugins.structalign;

import statalign.base.Utils;
import statalign.model.ext.plugins.StructAlign;

public class TranslationMove extends RotationOrTranslationMove {

	StructAlign owner;
	
	public TranslationMove (StructAlign s) {
		owner = s;
		omit = 0;
	}

	public double proposal(Object externalState) {
		for(int i = 0; i < 3; i++)
			owner.xlats[ind][i] = Utils.generator.nextGaussian() * owner.xlatP + owner.xlats[ind][i];  
				
		return 0;
		// logProposalRatio is 0 because prior is uniform and proposal is symmetric
	}

}
