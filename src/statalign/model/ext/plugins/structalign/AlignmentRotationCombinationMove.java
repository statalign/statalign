package statalign.model.ext.plugins.structalign;

import java.util.List;
import statalign.model.ext.McmcMove;
import statalign.model.ext.McmcCombinationMove;
import statalign.model.ext.plugins.StructAlign;

public class AlignmentRotationCombinationMove extends McmcCombinationMove {

	StructAlign owner;
		
	public AlignmentRotationCombinationMove (StructAlign s, 
			List<McmcMove> mcmcMoves) {
		if (mcmcMoves.size() < 2) {
			throw new IllegalArgumentException("McmcCombinationMove must contain at least two McmcMove objects");
		}
		owner = s;
		name = mcmcMoves.get(0).name;
		for (int i=1; i<mcmcMoves.size(); i++) {
			name += "_"+mcmcMoves.get(i).name;
		}
	}
	public void updateLikelihood(Object externalState) {
		oldll = owner.curLogLike;
		for (McmcMove mcmcMove : mcmcMoves) {
			mcmcMove.updateLikelihood(externalState);
		}
	}
}
