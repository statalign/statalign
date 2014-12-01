package statalign.model.ext.plugins.structalign;

import statalign.base.Utils;
import statalign.model.ext.plugins.StructAlign;

public class TranslationMove extends RotationOrTranslationMove {
	
	public TranslationMove (StructAlign s, String n) {
		owner = s;
		structAlign = s;
		name = n;
		proposalWidthControlVariable = structAlign.xlatP;
	}

	public double proposal(Object externalState) {
		double[] shift = new double[3];
		for(int i = 0; i < 3; i++)
			shift[i] = Utils.generator.nextGaussian() * proposalWidthControlVariable;
		for(int l = 0; l < subtreeLeaves.size(); l++){
			int j = subtreeLeaves.get(l);
			if (structAlign.coords[j] == null) continue;
			for(int i = 0; i < 3; i++)
				structAlign.xlats[j][i] += shift[i];  
		}	
		return 0;
		// logProposalRatio is 0 because prior is uniform and proposal is symmetric
	}

}
